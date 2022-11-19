use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use bio::alignment::pairwise::*;
use crossbeam::channel::bounded;
use pyo3::prelude::*;
use rand::distributions::{Distribution, Uniform};
use rand_distr::weighted_alias::WeightedAliasIndex;
use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

use std::thread;

// TODO: Next step is to speed it up by splitting up the generation of kmer pairs to another thread
// And keeping a buffer always full

static ALPHABET: [u8; 5] = *b"ACGTN";

fn convert_to_onehot(seq: &Vec<u8>) -> Vec<u8> {
    let mut onehot = vec![0; seq.len() * 5];
    for (i, c) in seq.iter().enumerate() {
        let idx = ALPHABET.iter().position(|&x| x == *c).unwrap();
        onehot[i * 5 + idx] = 1;
    }
    onehot
}

#[pyclass]
struct KmerGenerator {
    k: usize,
    dist: WeightedAliasIndex<u32>,
    scoring: Scoring<MatchParams>,
    aligner: Aligner<MatchParams>,
    scoring_dist: Uniform<u8>,
    substitution_dist: Uniform<usize>,
    threads: usize,
    buffer_size: usize,
    // rx: Option<crossbeam::channel::Receiver<(Vec<u8>, Vec<u8>, String, String, i32, i32)>>,
    rx: Option<crossbeam::channel::Receiver<(Vec<u8>, Vec<u8>, i32)>>,
    handles: Vec<std::thread::JoinHandle<()>>,
    seed: u64,
}

#[pymethods]
impl KmerGenerator {
    #[new]
    fn new() -> Self {
        println!("Version 1");
        let scoring = Scoring::from_scores(-5, -1, 1, -1);
        KmerGenerator {
            k: 21,
            dist: WeightedAliasIndex::new(vec![25, 25, 25, 25, 5]).unwrap(),
            aligner: Aligner::with_capacity_and_scoring(21, 21, scoring.clone()),
            scoring,
            scoring_dist: Uniform::new_inclusive(0, 21),
            substitution_dist: Uniform::new_inclusive(0, 21),
            threads: 1,
            buffer_size: 8 * 1024 * 4,
            rx: None,
            handles: Vec::with_capacity(64),
            seed: 42,
        }
    }

    fn set_seed(&mut self, seed: u64) {
        self.seed = seed;
    }

    fn set_k(&mut self, k: usize) {
        self.k = k;
        let aligner = Aligner::with_capacity_and_scoring(k, k, self.scoring.clone());
        self.aligner = aligner;
        self.scoring_dist = Uniform::new(0, k as u8);
        self.substitution_dist = Uniform::new(0, k);
    }

    fn set_scoring(
        &mut self,
        match_score: i32,
        mismatch_score: i32,
        gap_open_score: i32,
        gap_extend_score: i32,
    ) {
        let scoring = Scoring::from_scores(
            match_score,
            mismatch_score,
            gap_open_score,
            gap_extend_score,
        );
        let aligner = Aligner::with_capacity_and_scoring(self.k, self.k, scoring.clone());
        self.scoring = scoring;
        self.aligner = aligner;
    }

    fn set_threads(&mut self, threads: usize) {
        self.threads = threads;
    }

    fn set_weights(&mut self, weights: Vec<u32>) {
        assert!(weights.len() == 5);
        self.dist = WeightedAliasIndex::new(weights).unwrap();
    }

    fn set_capacity(&mut self, buffer_size: usize) {
        self.buffer_size = buffer_size;
    }

    fn start(&mut self, py: Python<'_>) {
        py.allow_threads(|| {
            let (tx, rx) = bounded(self.buffer_size);
            self.rx = Some(rx);

            for threadnum in 0..self.threads {
                let tx = tx.clone();
                let k = self.k;
                let scoring = self.scoring.clone();
                let scoring_dist = Uniform::new_inclusive(0, self.k);
                let substitution_dist = Uniform::new(0, self.k);
                let dist = self.dist.clone();
                let mut aligner = Aligner::with_capacity_and_scoring(21, 21, scoring.clone());

                let mut rng = Xoshiro256PlusPlus::seed_from_u64(self.seed);

                for _ in 0..threadnum {
                    rng.long_jump();
                }

                let handle = std::thread::spawn(move || {
                    let mut score;
                    let backoff = crossbeam::utils::Backoff::new();
                    let mut iter;

                    let mut k1: Vec<u8> = vec![0; k];
                    let mut k2;

                    let arch = pulp::Arch::new();

                    'main: loop {
                        backoff.reset();

                        iter = 0;

                        for x in 0..k {
                            k1[x] = ALPHABET[dist.sample(&mut rng) as usize];
                        }

                        k2 = k1.clone();

                        let target_score = scoring_dist.sample(&mut rng) as i32;

                        score = 0;

                        if target_score - score > 6 {
                            k2.remove(substitution_dist.sample(&mut rng));
                            k2.push(ALPHABET[dist.sample(&mut rng) as usize]);
                        }

                        score = k as i32 - aligner.local(&k1, &k2).score;

                        // If score diff is big enough, create a new random kmer
                        if target_score - score > 14 {
                            for x in 0..k {
                                k2[x] = ALPHABET[dist.sample(&mut rng) as usize];
                            }
                        }

                        while score != target_score {
                            k2[substitution_dist.sample(&mut rng)] =
                                ALPHABET[dist.sample(&mut rng) as usize];

                            score = k as i32 - aligner.local(&k1, &k2).score;

                            iter += 1;

                            if iter > 128 {
                                break;
                            }

                            if score > target_score {
                                k1 = k2.clone();
                                score = 0;
                            }
                        }

                        let mut msg = arch.dispatch(|| {
                            (
                                convert_to_onehot(&k1),
                                convert_to_onehot(&k2),
                                score,
                            )
                        });

                        while let Err(e) = tx.try_send(msg) {
                            match e {
                                crossbeam::channel::TrySendError::Full(msg_e) => {
                                    msg = msg_e;
                                    backoff.snooze();

                                    if backoff.is_completed() {
                                        thread::park();
                                        backoff.reset();
                                    }

                                }
                                crossbeam::channel::TrySendError::Disconnected(_) => {
                                    break 'main;
                                }
                            };
                        }
                    }
                });

                self.handles.push(handle);
            }
        });
    }

    fn shutdown(&mut self) {
        self.rx = None;
        for handle in self.handles.drain(..) {
            handle.join().unwrap();
        }
    }

    fn generate_pairs<'py>(&mut self, py: Python<'py>) -> &'py pyo3::types::PyList {
        let data = pyo3::types::PyList::empty(py);

        let mut count: usize = 0;

        let rx = self.rx.as_ref().unwrap();

        py.allow_threads(|| {
            for handle in self.handles.iter() {
                handle.thread().unpark();
            }
        });

        while count < 1024 {
            count = count.saturating_add(1);
            if let Ok(msg) = rx.try_recv() {
                data.append(msg).expect("Unable to append to PyList");
            } else {
                py.allow_threads(|| {
                    println!("Waiting for data");
                    std::thread::sleep(std::time::Duration::from_millis(64));
                });
            }
        }
        data
    }

}


#[pymodule]
fn beaker_kmer_generator(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<KmerGenerator>()?;
    Ok(())
}
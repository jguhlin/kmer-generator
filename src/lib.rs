use bio::alignment::pairwise::*;
use pyo3::prelude::*;
use rand::distributions::{Distribution, Uniform};
use rand_distr::weighted_alias::WeightedAliasIndex;
use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

static ALPHABET: [u8; 5] = *b"ACGTN";

fn convert_to_onehot(seq: &str) -> Vec<u8> {
    let mut onehot = vec![0; seq.len() * 5];
    for (i, c) in seq.chars().enumerate() {
        let idx = ALPHABET.iter().position(|&x| x == c as u8).unwrap();
        onehot[i * 5 + idx] = 1;
    }
    onehot
}

#[pyclass]
struct KmerGenerator {
    rng: Xoshiro256PlusPlus,
    k: usize,
    dist: WeightedAliasIndex<u32>,
    scoring: Scoring<MatchParams>,
    aligner: Aligner<MatchParams>,
    scoring_dist: Uniform<u8>,
    substitution_dist: Uniform<usize>,
}

#[pymethods]
impl KmerGenerator {
    #[new]
    fn new() -> Self {
        let scoring = Scoring::from_scores(-5, -1, 1, -1);
        KmerGenerator {
            rng: Xoshiro256PlusPlus::seed_from_u64(42),
            k: 21,
            dist: WeightedAliasIndex::new(vec![25, 25, 25, 25, 5]).unwrap(),
            aligner: Aligner::with_capacity_and_scoring(21, 21, scoring.clone()),
            scoring,
            scoring_dist: Uniform::new(0, 21),
            substitution_dist: Uniform::new(0, 21),
        }
    }

    fn set_seed(&mut self, seed: u64) {
        self.rng = Xoshiro256PlusPlus::seed_from_u64(seed);
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

    fn set_weights(&mut self, weights: Vec<u32>) {
        assert!(weights.len() == 5);
        self.dist = WeightedAliasIndex::new(weights).unwrap();
    }

    fn generate_kmer(&mut self) -> String {
        let mut kmer = String::with_capacity(self.k);
        for _ in 0..self.k {
            kmer.push(ALPHABET[self.dist.sample(&mut self.rng) as usize] as char);
        }
        kmer
    }

    fn generate_pair(&mut self) -> (Vec<u8>, Vec<u8>, String, String, i32) {

        let mut k1 = String::with_capacity(self.k);
        let mut k2; // = String::with_capacity(self.k);

        let target_score = self.scoring_dist.sample(&mut self.rng) as i32;
        let mut score = 0;

        for _ in 0..self.k {
            k1.push(ALPHABET[self.dist.sample(&mut self.rng) as usize] as char);
            // k2.push(ALPHABET[self.dist.sample(&mut self.rng) as usize] as char);
        }

        k2 = k1.clone();

        let mut iter = 0;

        while score < target_score {
            let to_replace = self.substitution_dist.sample(&mut self.rng);
            k2.replace_range(
                to_replace..to_replace+1,
                std::str::from_utf8(&[ALPHABET[self.dist.sample(&mut self.rng) as usize]]).unwrap(),
            );
            score = self.k as i32 - self.aligner.local(&k1.as_bytes(), &k2.as_bytes()).score;

            iter += 1;
            if iter > self.k * 2 {
                break;
            }
        }

        

        (convert_to_onehot(&k1), convert_to_onehot(&k2), k1, k2, score)
    }
}

#[pymodule]
fn beaker_kmer_generator(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<KmerGenerator>()?;
    Ok(())
}

/*
fn main() {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
    let weights = vec![25, 25, 25, 25, 5];
    let dist = WeightedAliasIndex::new(weights).unwrap();

    let k = 21;

    let scoring = Scoring::from_scores(-5, -1, 1, -1);
    let max_score = k as i32;

    let mut aligner = Aligner::with_capacity_and_scoring(21, 21, scoring);

    for _ in 0..10 {
        let k1 = generate_kmer(&dist, &mut rng, k);
        let k2 = generate_kmer(&dist, &mut rng, k);
        let mut k3 = k1.clone();
        k3.replace_range(20..21, "C");

        let alignment = aligner.local(&k1.as_bytes(), &k2.as_bytes());

        println!("{} {} Score: {}", k1, k2, max_score - alignment.score);

        let alignment = aligner.local(&k1.as_bytes(), &k3.as_bytes());

        println!("{} {} Score: {}", k1, k3, max_score - alignment.score);
    }
} */

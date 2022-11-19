from beaker_kmer_generator import KmerGenerator as kmer_generator
import plotille

kg = kmer_generator()
kg.set_threads(32)
kg.set_seed(42)
kg.start()

j = list()
x = list()
for _ in range(0, 100000):
    if len(j) >= 1000000:
        break
    js = kg.generate_pairs()
    for i in js:
        j.append(i[2])
        x.append(i[2])

fig = plotille.hist(
    j,
    bins=19,
    width=80,
    log_scale=False,
    linesep="\n",
    lc=None,
    bg=None,
    color_mode="names",
)

print(fig)

fig = plotille.hist(
    x,
    bins=19,
    width=80,
    log_scale=False,
    linesep="\n",
    lc=None,
    bg=None,
    color_mode="names",
)

print(fig)

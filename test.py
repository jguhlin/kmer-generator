from beaker_kmer_generator import KmerGenerator as kmer_generator
# mport beaker_kmer_generator as bkg

kg = kmer_generator()

for _ in range(0, 1000):
  print(kg.generate_pair()[4])
#  print(kg.generate_pair()[4])

from mcomp_project.algorithms.markov_chain import DNAMarkovChain
from Bio import SeqIO
import numpy as np
import time

def test_markov_chain(fasta_file_name, pickle_file_name=None, test_num=10000, order=6, region_length=10000, read_length=100, substitution_rate=0.02, prior=0.001):
    mc = DNAMarkovChain(order, region_length, read_length, substitution_rate, prior)
    start_time = time.time()
    if pickle_file_name is None:
        mc.read(fasta_file_name)
    else:
        mc.load(pickle_file_name)
    print("Reading time:", time.time() - start_time, "s")
    fasta_sequences = SeqIO.parse(open(fasta_file_name), 'fasta')
    for fasta in fasta_sequences:
        sequence_length = len(str(fasta.seq))
        correct = 0
        start_time = time.time()
        for _ in range(test_num):
            index = np.random.randint(0, sequence_length - read_length - 1)
            region = index // region_length
            #print(fasta.seq[index:index + read_length])
            predict_region = mc.query(fasta.seq[index:index + read_length])
            correct += (region == predict_region)
            #print(region, predict_region)
    
        print("Query time:", time.time() - start_time, "s")
    
    print("Correct rate:", correct / test_num)

if __name__ == "__main__":
    test_markov_chain("/home/zhenhao/mcomp-dissertation/sequence_sample.fasta")#, "/home/zhenhao/mcomp-dissertation/sequence_sample.fasta_markov_chain.pickle")
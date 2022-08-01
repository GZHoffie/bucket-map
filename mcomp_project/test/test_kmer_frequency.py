from mcomp_project.algorithms.kmer_frequency import KMerFrequency, GappedKMerFrequency
from Bio import SeqIO
import numpy as np
import time
import math
from sklearn.decomposition import PCA
from mcomp_project.utils import *

from matplotlib import pyplot as plt


def test_kmer_frequency(fasta_file_name, pickle_file_name=None, test_num=1000, order=7, bucket_length=10000, read_length=100, substitution_rate=0.02, prior=0.001):
    mc = KMerFrequency(order, bucket_length, read_length, substitution_rate, prior)
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
            bucket = index // bucket_length
            predict_bucket = mc.query(fasta.seq[index:index + read_length])
            correct += (bucket == predict_bucket)
            #print(bucket, predict_bucket)
    
        print("Query time:", time.time() - start_time, "s")
    
    print("Correct rate:", correct / test_num)


def test_gapped_kmer_frequency(fasta_file_name, pickle_file_name=None, test_num=100000, order=8, bucket_length=10000, read_length=100, sample_num=30, substitution_rate=0, prior=0.001, test_reads=False):
    mc = GappedKMerFrequency(order, bucket_length, read_length, sample_num, substitution_rate, prior, 4)
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

        hit_probability = np.zeros(test_num)
        other_largest_probability = np.zeros(test_num)

        for i in range(test_num):
            index = np.random.randint(0, sequence_length - read_length - 1)
            bucket = index // bucket_length
            test_sequence = fasta.seq[index:index + read_length]
            if test_reads:
                test_sequence = modify_sequence(test_sequence, math.floor(substitution_rate * len(test_sequence)))
            
            predict_bucket, log_probability = mc.query(test_sequence)

            hit_probability[i] = log_probability[bucket]
            other_largest_probability[i] = max(list(log_probability[:bucket]) + list(log_probability[bucket+1:]))
            correct += (bucket == predict_bucket)
            #print(bucket, predict_bucket)
        total_time = time.time() - start_time
        print("Query time:", total_time, "s")
        print("Throughput:", read_length * test_num / total_time, "bp/s")

    print("Expectation for hit bucket:", np.mean(hit_probability))
    print("Expectation for other buckets:", np.mean(other_largest_probability))
    plt.hist(hit_probability, bins='sturges')
    plt.hist(other_largest_probability, bins='sturges')
    plt.show()
    
    print("Correct rate:", correct / test_num)

def try_PCA(fasta_file_name, pickle_file_name, test_num=1000, order=7, bucket_length=20000, read_length=100, sample_num=50, substitution_rate=0.02, prior=0.001):
    mc = KMerFrequency(order, bucket_length, read_length, sample_num, substitution_rate, prior)
    mc.load(pickle_file_name)

    matrix = mc.M
    print(matrix.shape)
    pca = PCA(n_components='mle')
    pca.fit_transform(matrix.T)
    print(pca.explained_variance_ratio_.sum())





if __name__ == "__main__":
    #test_kmer_frequency("/home/zhenhao/mcomp-dissertation/sequence_sample.fasta")#, "/home/zhenhao/mcomp-dissertation/sequence_sample.fasta_markov_chain.pickle")
    test_gapped_kmer_frequency("/home/zhenhao/mcomp-dissertation/sequence_sample.fasta", "/home/zhenhao/mcomp-dissertation/sequence_sample.fasta_frequency_list_gapped.pickle", test_reads=False)
    
    #try_PCA("/home/zhenhao/mcomp-dissertation/sequence_sample.fasta", "/home/zhenhao/mcomp-dissertation/sequence_sample.fasta_markov_chain.pickle")
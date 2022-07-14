from mcomp_project.algorithms.kmer_frequency import KMerFrequency, GappedKMerFrequency
from Bio import SeqIO
import numpy as np
import time
from sklearn.decomposition import PCA

def test_kmer_frequency(fasta_file_name, pickle_file_name=None, test_num=1000, order=7, region_length=10000, read_length=100, substitution_rate=0.02, prior=0.001):
    mc = KMerFrequency(order, region_length, read_length, substitution_rate, prior)
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
            predict_region = mc.query(fasta.seq[index:index + read_length])
            correct += (region == predict_region)
            #print(region, predict_region)
    
        print("Query time:", time.time() - start_time, "s")
    
    print("Correct rate:", correct / test_num)


def test_gapped_kmer_frequency(fasta_file_name, pickle_file_name=None, test_num=1000, order=7, region_length=10000, read_length=100, substitution_rate=0.02, prior=0.001):
    mc = GappedKMerFrequency(order, region_length, read_length, substitution_rate, prior)
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
            predict_region = mc.query(fasta.seq[index:index + read_length])
            correct += (region == predict_region)
            #print(region, predict_region)
    
        print("Query time:", time.time() - start_time, "s")
    
    print("Correct rate:", correct / test_num)

def try_PCA(fasta_file_name, pickle_file_name, test_num=1000, order=7, region_length=10000, read_length=100, substitution_rate=0.02, prior=0.001):
    mc = KMerFrequency(order, region_length, read_length, substitution_rate, prior)
    mc.load(pickle_file_name)

    matrix = mc.M
    print(matrix.shape)
    pca = PCA(n_components='mle')
    pca.fit_transform(matrix.T)
    print(pca.explained_variance_ratio_.sum())





if __name__ == "__main__":
    #test_kmer_frequency("/home/zhenhao/mcomp-dissertation/sequence_sample.fasta")#, "/home/zhenhao/mcomp-dissertation/sequence_sample.fasta_markov_chain.pickle")
    test_gapped_kmer_frequency("/home/zhenhao/mcomp-dissertation/sequence_sample.fasta")#, "/home/zhenhao/mcomp-dissertation/sequence_sample.fasta_markov_chain.pickle")
    
    #try_PCA("/home/zhenhao/mcomp-dissertation/sequence_sample.fasta", "/home/zhenhao/mcomp-dissertation/sequence_sample.fasta_markov_chain.pickle")
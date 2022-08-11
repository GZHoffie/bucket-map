from mcomp_project.utils import *
import numpy as np

class KMerExistence(BucketMapper):
    def __init__(self, order=7, bucket_length=100000, read_length=100, sample_num=10, tolerant_num=2) -> None:
        super().__init__(bucket_length, read_length)
        self.order = order
        self.sample_num = sample_num
        self.tolerant_num = tolerant_num

        # Storage data structure
        self.kmer_record = None

        self.compressed_kmer_record = []
    
    def read_bucket(self, bucket):
        kmer_record = np.zeros(4 ** self.order, dtype=bool)
        for i in range(len(bucket) - self.order + 1):
            kmer_record[kmer_to_index(bucket[i:i+self.order])] = True
        
        self.kmer_record = kmer_record if self.kmer_record is None else np.vstack([self.Mkmer_record, kmer_record])
    
    def end_read(self):
        for i in range(4 ** self.order):

            kmer_string = self.kmer_record[:, i] * i
            kmer_string = ''.join(kmer_string.astype(str))
            self.compressed_kmer_record.append(int(kmer_string, 2))
        
        self.kmer_record = None
    
    def query(self, sequence):
        tolerate_vectors = [int('1' * (4 ** self.order), 2)] * self.tolerant_num
        sample_indices = np.linspace(0, len(sequence) - self.order)
        for i in sample_indices:
            tolerate_vectors[0] = tolerate_vectors[0] | self.compressed_kmer_record[]
            for j in range(self.tolerant_num):

    
    def store(self, pickle_file_name):
        """
        Store the self.M matrix into a pickle file.
        """
        with open(pickle_file_name, "wb") as f:
            pickle.dump((self.compressed_kmer_record, self.bucket_id), f)
    
    def load(self, pickle_file_name):
        """
        load the pickle file into self.M.
        """
        with open(pickle_file_name, 'rb') as f:
            self.compressed_kmer_record, self.bucket_id = pickle.load(f)

            
import random
import math
import pickle
import time
from Bio import SeqIO
import numpy as np

char_to_index_map = {
    "A": 0, "a": 0,
    "C": 1, "c": 1,
    "G": 2, "g": 2,
    "T": 3, "t": 3
}

def kmer_to_index(sequence):
    """
    Turn a `self.order`-mer into an index between 0 and
    4 ** `self.order` - 1. Used to store the frequency of the k-mer.

    Args:
        sequence: A k-mer.
        
    Returns:
        The index of the sequence or None if illegal characters appear
        in the sequence.
    """
    #assert len(sequence) == self.order, f"The length of input sequence {sequence} should be of length {self.order}"
    try:
        index = [char_to_index_map[c] for c in sequence]
    except:
        return None
    return np.ravel_multi_index(index, [4] * len(sequence))


def modify_sequence(sequence, num_errors, substitution_percentage=0.95):
    """
    Randomly modify a sequence: add several mismatches or indels, according
    to the specified `num_errors`.

    Args:
        sequence: A DNA sequence.
        num_errors: The total number of errors to be introduced.
        substitution_percentage: the percentage of mismatches in all the errors
            (mismatch / (mismatch + indel)).
    """
    sequence_copy = sequence
    #print("--", sequence_copy)
    for _ in range(num_errors):
        random_number = random.uniform(0, 1)
        if random_number < substitution_percentage:
            # add a mismatch in the sequence
            i = random.randint(0, len(sequence_copy) - 1)
            current_char = sequence_copy[i]
            candidate_char = ['A', 'C', 'G', 'T']
            candidate_char.remove(current_char)
            sequence_copy = sequence_copy[:i] + random.choice(candidate_char) + sequence_copy[i+1:]
        elif random_number < substitution_percentage + (1-substitution_percentage)/2:
            # add an insert
            i = random.randint(0, len(sequence_copy))
            sequence_copy = sequence_copy[:i] + random.choice(['A', 'C', 'G', 'T']) + sequence_copy[i:]
        else:
            # add a delete
            i = random.randint(0, len(sequence_copy) - 1)
            sequence_copy = sequence_copy[:i] + sequence_copy[i+1:]
        
        #print("->", sequence_copy)
        
    return sequence_copy


class BucketMapper:
    """
    Base class for the bucket mapper algorithm. This class includes
    basic functionalities such as dividing the reference genome
    into buckets, finding an approximate match of a read inside a given
    bucket, etc.
    """
    def __init__(self, bucket_length=100000, read_length=100) -> None:
        self.bucket_length = bucket_length
        self.read_length = read_length

        # A list storing where does the bucket come from
        self.bucket_id = []

    def read_bucket(self, bucket):
        """
        Read and store the sequence in the bucket in the indexing
        data structure.

        Args:
            bucket: A sequence or a string of length self.bucket_length + self.read_length
                each bucket represent all reads that originate from this area.
        """
        pass

    def read_fasta(self, fasta_file_name):
        """
        Read the fasta file and divide the reference genome into
        buckets. Store the bucket id in self.bucket_id.
        """
        print("Opening FASTA file:", fasta_file_name)
        fasta_sequences = SeqIO.parse(open(fasta_file_name), 'fasta')
        for fasta in fasta_sequences:
            sequence_length = len(str(fasta.seq))
            print("... processing sequence", fasta.id, "with length", sequence_length)
            
            num_buckets = math.ceil(float(sequence_length) / self.bucket_length)


            for bucket in range(num_buckets):
                start_index = bucket * self.bucket_length
                end_index = (bucket + 1) * self.bucket_length + self.read_length + 1
                bucket_sequence = fasta.seq[start_index:end_index]
                self.read_bucket(bucket_sequence)
                self.bucket_id.append((fasta.id, bucket))
        
        self.end_read()
    

    def query(self, sequence):
        pass

    def end_read(self):
        pass

    def store(self, pickle_file_name):
        pass
    
    def load(self, pickle_file_name):
        pass
    

    def test(self, fasta_file_name, test_num=100000, test_reads=False, error_rate=0.02):
        fasta_sequences = SeqIO.parse(open(fasta_file_name), 'fasta')
        for fasta in fasta_sequences:
            sequence_length = len(str(fasta.seq))
            correct = 0
            start_time = time.time()

            hit_probability = np.zeros(test_num)
            other_largest_probability = np.zeros(test_num)

            for i in range(test_num):
                index = np.random.randint(0, sequence_length - self.read_length - 1)
                bucket = index // self.bucket_length
                test_sequence = fasta.seq[index:index + self.read_length]
                if test_reads:
                    test_sequence = modify_sequence(test_sequence, math.floor(error_rate * len(test_sequence)))
            
                predict_bucket, log_probability = self.query(test_sequence)

                hit_probability[i] = log_probability[bucket]
                other_largest_probability[i] = max(list(log_probability[:bucket]) + list(log_probability[bucket+1:]))
                correct += (bucket == predict_bucket)
                #print(bucket, predict_bucket)
            total_time = time.time() - start_time
            print("Query time:", total_time, "s")
            print("Throughput:", self.read_length * test_num / total_time, "bp/s")

        print("Correct rate:", correct / test_num)


if __name__ == "__main__":
    print(modify_sequence("ACGCT", 2))
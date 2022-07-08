from Bio import SeqIO
from mcomp_project.utils import *
import numpy as np
import pickle
import math

class DNAMarkovChain:
    def __init__(self, order=7, region_length=100000, read_length=100, substitution_rate=0.02, prior=0.001) -> None:
        self.order = order
        self.region_length = region_length
        self.read_length = read_length
        self.substitution_rate = substitution_rate
        self.prior = prior
        self.pi = None
        self.A = None
    

    def _new_markov_chain(self):
        """
        Initialize a new vector and fill everything with `self.prior`
        """
        pi = np.zeros([4] * self.order)
        A = np.zeros([4] * (self.order + 1))
        pi.fill(self.prior)
        A.fill(self.prior)
        return pi, A
    

    def _insert_into_A(self, sequence, A, allow_substitution=True):
        """
        Count the sequence into the markov chain. In particular, we want
        markov_chain[sequence] += 1. If we allow substitution, we also do
        markov_chain[variant_of_sequence] += self.substitution_rate.

        The variant of a sequence is basically all the sequences that can be 
        obtained by substituting a character in the original sequence.

        Args:
            sequence: A sequence of length `self.order`
            markov_chain: A matrix storing markov chain parameters
            allow_substitution: whether we allow substitution.
        """
        try:
            index = [char_to_index_map[c] for c in sequence]
        except:
            return
        if allow_substitution:
            for i in range(len(sequence)):
                index_i = index[i]
                index[i] = ...
                A[tuple(index)] += self.substitution_rate
                index[i] = index_i
            
            A[tuple(index)] += 1 - self.order * self.substitution_rate
        else:
            A[tuple(index)] += 1
    

    def _insert_into_pi(self, sequence, pi, allow_substitution=True):

        try:
            index = [char_to_index_map[c] for c in sequence]
        except:
            return
        if allow_substitution:
            for i in range(len(sequence)):
                index_i = index[i]
                index[i] = ...
                pi[tuple(index)] += self.substitution_rate
                index[i] = index_i
            
            pi[tuple(index)] += 1 - self.order * self.substitution_rate
        else:
            pi[tuple(index)] += 1



        
    
    def _sequence_to_index(self, sequence):
        """
        Turn a `self.order`-mer into an index between 0 and
        4 ** `self.order` - 1. Used to store the frequency of the k-mer.

        Args:
            sequence: A sequence of length `self.order`
        
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
    

    def _store_markov_chain(self, pi, A):
        """
        Store the current markov chain parameters in self.pi and self.A. We do a
        normalization and store the log of each parameters.
        """
        normalized_pi = pi / np.sum(pi)
        normalized_pi = np.log(normalized_pi)
        self.pi = normalized_pi.flatten() if self.pi is None else np.vstack([self.pi, normalized_pi.flatten()])

        # Normalize A with respect to the last axis
        #print(np.sum(A, axis=-1, keepdims=True).shape)
        normalized_A = A / np.sum(A, axis=0, keepdims=True)
        #print(normalized_A.sum(axis=-1))
        normalized_A = np.log(normalized_A)
        self.A = normalized_A.flatten() if self.A is None else np.vstack([self.A, normalized_A.flatten()])


            

    def read(self, fasta_file_name, output_file=True):
        """
        Read the fasta file and store the Markov Chain frequencies in
        self.A and self.pi, and print the related information.

        Args:
            fasta_file_name: the string consisting of the path to the
                fasta file.
        """
        print("Opening FASTA file:", fasta_file_name)
        fasta_sequences = SeqIO.parse(open(fasta_file_name), 'fasta')
        for fasta in fasta_sequences:
            sequence_length = len(str(fasta.seq))
            print("... processing sequence", fasta.id, "with length", sequence_length)
            print("Estimated space usage: ", sequence_length / self.region_length * (4 ** self.order) * 4 / (1024 ** 2), "MB")
            
            # start counting the frequency the k-mers
            num_chunks = math.ceil(float(sequence_length) / self.region_length)


            for chunk in range(num_chunks):
                # Initialize pi and A
                pi, A = self._new_markov_chain()

                # First step: store the first read in the chunk to pi and A
                index = chunk * self.region_length
                self._insert_into_pi(fasta.seq[index:index + self.order], pi)
                for i in range(1, self.read_length + 1):
                    self._insert_into_A(fasta.seq[index + i:index + i + self.order + 1], A)
                
                # Repeat until we read all k-mers in the chunk
                for i in range(1, self.region_length):
                    self._insert_into_pi(fasta.seq[index + i:index + i + self.order], pi)
                    self._insert_into_A(fasta.seq[index + i + self.read_length:index + i + self.read_length + self.order], A)
                
                # Complete this region, store the markov chain parameters
                self._store_markov_chain(pi, A)

        #print(self.M)
        if output_file:
            with open(f"{fasta_file_name}_markov_chain.pickle", "wb") as f:
                pickle.dump((self.pi, self.A), f, pickle.HIGHEST_PROTOCOL)
    
    def store(self, pickle_file_name):
        """
        Store the self.M matrix into a pickle file.
        """
        with open(pickle_file_name, "wb") as f:
            pickle.dump((self.pi, self.A), f, pickle.HIGHEST_PROTOCOL)
    
    def load(self, pickle_file_name):
        """
        load the pickle file into self.M.
        """
        with open(pickle_file_name, 'rb') as f:
            self.pi, self.A = pickle.load(f)
    
    def query(self, sequence):
        # Find the first k-mer's corresponding probability
        first_k_mer = sequence[:self.order]
        init_probability = self.pi[:, self._sequence_to_index(first_k_mer)]
        #print(init_probability)

        k_mers = np.zeros((4 ** (self.order + 1), 1))
        for i in range(1, len(sequence) - self.order - 1):
            current_k_mer = sequence[i:i + self.order + 1]
            index = self._sequence_to_index(current_k_mer)
            if index is not None:
                k_mers[index] += 1
        
        #print(np.dot(self.A, k_mers).flatten())
        #print(init_probability)
        log_probability = np.dot(self.A, k_mers).flatten() + init_probability
        return np.argmax(log_probability)


        

            
        



if __name__ == "__main__":
    mc = DNAMarkovChain(order=6, region_length=20000)
    #mc._insert_into_markov_chain("ACNCN", mc._new_markov_chain())
    #mc.read("/home/zhenhao/data/SRR611076/sequence.fasta")
    #print(mc.query("GCTCTTTCCCCGGAAACCATTGAAATCGGACGGTTTAGTGAAAATGGAGGATCAAGTTGGGTTTGGGTTC"))
        
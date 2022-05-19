from Bio import SeqIO
from mcomp_project.utils import *
import numpy as np
import pickle

class DNAMarkovChain:
    def __init__(self, order=7, region_length=100000, read_length=100, substitution_rate=0.02, prior=0.001) -> None:
        self.order = order
        self.region_length = region_length
        self.read_length = read_length
        self.substitution_rate = substitution_rate
        self.prior = prior
        self.M = None
    

    def _new_markov_chain(self):
        """
        Initialize a new vector and fill everything with `self.prior`
        """
        markov_chain = np.zeros([4] * self.order)
        markov_chain.fill(self.prior)
        return markov_chain
    

    def _insert_into_markov_chain(self, sequence, markov_chain, allow_substitution=True):
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
            for i in range(self.order):
                index_i = index[i]
                index[i] = ...
                markov_chain[tuple(index)] += self.substitution_rate
                index[i] = index_i
            
            markov_chain[tuple(index)] += 1 - self.order * self.substitution_rate
        else:
            markov_chain[tuple(index)] += 1


        
    
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
        assert len(sequence) == self.order, f"The length of input sequence {sequence} should be of length {self.order}"
        try:
            index = [char_to_index_map[c] for c in sequence]
        except:
            return None
        return np.ravel_multi_index(index, [4] * self.order)
    

    def _store_markov_chain(self, markov_chain):
        """
        Store the current markov chain parameters in self.M. We do a
        normalization and store the log of each parameters.

        Args:
            markov_chain: an array of markov chain parameters
        """
        markov_chain /= np.sum(markov_chain)
        markov_chain = np.log(markov_chain)
        self.M = markov_chain.flatten() if self.M is None else np.vstack([self.M, markov_chain.flatten()])

            

    def read(self, fasta_file_name, output_file=True):
        """
        Read the fasta file and store the Markov Chain frequencies in
        self.M, and print the related information.

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
            index = 0
            markov_chain = self._new_markov_chain()
            while index + self.order < sequence_length:

                current_k_mer = fasta.seq[index:index + self.order]
                self._insert_into_markov_chain(current_k_mer, markov_chain)
                #print(current_k_mer, self._sequence_to_index(current_k_mer))
                index += 1
                if index % self.region_length == 0:
                    # Complete this region, store the markov chain parameters
                    self._store_markov_chain(markov_chain)
                    
                    # Create new Markov chain
                    markov_chain = self._new_markov_chain()
            
            # Complete this region, store the markov chain parameters
            self._store_markov_chain(markov_chain)
                    
        
        #print(self.M)
        if output_file:
            with open(f"{fasta_file_name}_markov_chain.pickle", "wb") as f:
                pickle.dump(self.M, f, pickle.HIGHEST_PROTOCOL)
    
    def store(self, pickle_file_name):
        """
        Store the self.M matrix into a pickle file.
        """
        with open(pickle_file_name, "wb") as f:
            pickle.dump(self.M, f)
    
    def load(self, pickle_file_name):
        """
        load the pickle file into self.M.
        """
        with open(pickle_file_name, 'rb') as f:
            self.M = pickle.load(f)
    
    def query(self, sequence):
        k_mers = np.zeros((4 ** self.order, 1))
        for i in range(len(sequence) - self.order):
            current_k_mer = sequence[i:i + self.order]
            index = self._sequence_to_index(current_k_mer)
            if index is not None:
                k_mers[index] += 1
        
        log_probability = np.dot(self.M, k_mers)
        #print(log_probability)
        return np.argmax(log_probability)

        

            
        



if __name__ == "__main__":
    mc = DNAMarkovChain(order=7, region_length=20000)
    mc._insert_into_markov_chain("ACNCN", mc._new_markov_chain())
    #mc.read("/home/zhenhao/data/SRR611076/sequence.fasta")
    #print(mc.query("GCTCTTTCCCCGGAAACCATTGAAATCGGACGGTTTAGTGAAAATGGAGGATCAAGTTGGGTTTGGGTTC"))
        
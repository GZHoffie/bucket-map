from typing import overload
from Bio import SeqIO
from mcomp_project.utils import *
import numpy as np
import pickle
from typing import List

class KMerFrequency:
    def __init__(self, order=7, region_length=100000, read_length=100, substitution_rate=0.02, prior=0.001) -> None:
        self.order = order
        self.region_length = region_length
        self.read_length = read_length
        self.substitution_rate = substitution_rate
        self.prior = prior
        self.M = None
    

    def _new_frequency_list(self):
        """
        Initialize a new vector and fill everything with `self.prior`
        """
        frequency_list = np.zeros([4] * self.order)
        frequency_list.fill(self.prior)
        return frequency_list
    

    def _insert_into_frequency_list(self, sequence, frequency_list, allow_substitution=True):
        """
        Count the sequence into the markov chain. In particular, we want
        frequency_list[sequence] += 1. If we allow substitution, we also do
        frequency_list[variant_of_sequence] += self.substitution_rate.

        The variant of a sequence is basically all the sequences that can be 
        obtained by substituting a character in the original sequence.

        Args:
            sequence: A sequence of length `self.order`
            frequency_list: A matrix storing markov chain parameters
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
                frequency_list[tuple(index)] += self.substitution_rate
                index[i] = index_i
            
            frequency_list[tuple(index)] += 1 - self.order * self.substitution_rate
        else:
            frequency_list[tuple(index)] += 1


        
    
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
    

    def _store_frequency_list(self, frequency_list):
        """
        Store the current markov chain parameters in self.M. We do a
        normalization and store the log of each parameters.

        Args:
            frequency_list: an array of markov chain parameters
        """
        frequency_list /= np.sum(frequency_list)
        frequency_list = np.log(frequency_list)
        self.M = frequency_list.flatten() if self.M is None else np.vstack([self.M, frequency_list.flatten()])

            

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
            frequency_list = self._new_frequency_list()
            while index + self.order < sequence_length:

                current_k_mer = fasta.seq[index:index + self.order]
                self._insert_into_frequency_list(current_k_mer, frequency_list)
                #print(current_k_mer, self._sequence_to_index(current_k_mer))
                index += 1
                if index % self.region_length == 0:
                    # Complete this region, store the markov chain parameters
                    self._store_frequency_list(frequency_list)
                    
                    # Create new Markov chain
                    frequency_list = self._new_frequency_list()
            
            # Complete this region, store the markov chain parameters
            self._store_frequency_list(frequency_list)
                    
        
        #print(self.M)
        if output_file:
            with open(f"{fasta_file_name}_frequency_list.pickle", "wb") as f:
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
            self.M = pickle.load(f)[0]
    
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

        
class GappedKMerFrequency(KMerFrequency):
    def __init__(self, order=7, region_length=100000, read_length=100, substitution_rate=0.02, prior=0.001, gapped_k_mer_sequence=5) -> None:
        super().__init__(order, region_length, read_length, substitution_rate, prior)
        if isinstance(gapped_k_mer_sequence, int):
            # Randomly generate a gapped k-mer sequence between 0-(order+gapped_k_mer_sequence)
            self.gapped_k_mer = random.sample(range(order + gapped_k_mer_sequence), k=order)
            self.gapped_k_mer.sort()
        elif isinstance(gapped_k_mer_sequence, List):
            self.gapped_k_mer = gapped_k_mer_sequence
        
        print("Using gapped k-mer sequence:", self.gapped_k_mer)
    

    def _get_current_k_mer(self, sequence, starting_index):
        """
        Return the current gapped k-mer.
        """
        return "".join([sequence[i + starting_index] for i in self.gapped_k_mer])

            

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
            frequency_list = self._new_frequency_list()
            while index + np.max(self.gapped_k_mer) < sequence_length:

                current_k_mer = self._get_current_k_mer(fasta.seq, index)
                self._insert_into_frequency_list(current_k_mer, frequency_list)
                #print(current_k_mer, self._sequence_to_index(current_k_mer))
                index += 1
                if index % self.region_length == 0:
                    # Complete this region, store the markov chain parameters
                    self._store_frequency_list(frequency_list)
                    
                    # Create new Markov chain
                    frequency_list = self._new_frequency_list()
            
            # Complete this region, store the markov chain parameters
            self._store_frequency_list(frequency_list)
                    
        
        #print(self.M)
        if output_file:
            with open(f"{fasta_file_name}_frequency_list_gapped.pickle", "wb") as f:
                pickle.dump(self.M, f, pickle.HIGHEST_PROTOCOL)
    

    def query(self, sequence):
        k_mers = np.zeros((4 ** self.order, 1))
        for i in range(len(sequence) - np.max(self.gapped_k_mer)):
            current_k_mer = self._get_current_k_mer(sequence, i)
            index = self._sequence_to_index(current_k_mer)
            if index is not None:
                k_mers[index] += 1
        
        log_probability = np.dot(self.M, k_mers)
        #print(log_probability)
        return np.argmax(log_probability)
    
    def store(self, pickle_file_name):
        """
        Store the self.M matrix into a pickle file.
        """
        with open(pickle_file_name, "wb") as f:
            pickle.dump((self.M, self.gapped_k_mer), f)
    
    def load(self, pickle_file_name):
        """
        load the pickle file into self.M.
        """
        with open(pickle_file_name, 'rb') as f:
            self.M, self.gapped_k_mer = pickle.load(f)
    


        



if __name__ == "__main__":
    mc = KMerFrequency(order=7, region_length=20000)
    mc._insert_into_frequency_list("ACNCN", mc._new_frequency_list())
    #mc.read("/home/zhenhao/data/SRR611076/sequence.fasta")
    #print(mc.query("GCTCTTTCCCCGGAAACCATTGAAATCGGACGGTTTAGTGAAAATGGAGGATCAAGTTGGGTTTGGGTTC"))
        
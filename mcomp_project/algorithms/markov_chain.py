from Bio import SeqIO
from mcomp_project.utils import *
import numpy as np

class DNAMarkovChain:
    def __init__(self, order=5, segment_length=100000, prior=0.1) -> None:
        self.order = order
        self.region_length = segment_length
        self.prior = prior
        self.M = None
    

    def _new_markov_chain(self):
        """
        Initialize a new vector and fill everything with `self.prior`
        """
        markov_chain = np.zeros((1, 4 ** self.order))
        markov_chain.fill(self.prior)
        return markov_chain
        
    
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
        index = 0
        for index, char in enumerate(sequence):
            try:
                index += (4 ** index) * char_to_index_map[char]
            except:
                print(f"[WARNING]\tillegal char {char} occurred, ignored.")
                return None
        
        return index
            
    
    
    def read(self, fasta_file_name):
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
            
            # start putting the sequence into
            index = 0

            current_k_mer = "ACNCT"#fasta.seq[index:index + self.order]
            print(current_k_mer, self._sequence_to_index(current_k_mer))



if __name__ == "__main__":
    mc = DNAMarkovChain()
    mc.read("/home/zhenhao/data/SRR611076/sequence.fasta")
        
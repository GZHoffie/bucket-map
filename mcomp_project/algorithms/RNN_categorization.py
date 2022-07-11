from torch.utils.data import Dataset
from Bio import SeqIO
import random


class ShortReadsDataset(Dataset):
    def __init__(self, fasta_file_name, substitution_rate=0.02, indel_rate=0.001) -> None:
        self.substitution_rate = substitution_rate
        self.indel_rate = self.indel_rate
        # Read the sequence and store
        fasta_sequences = SeqIO.parse(open(fasta_file_name), 'fasta')
        for fasta in fasta_sequences:
            self.sequence = fasta.seq
    
    def __len__(self):
        return len(str(self.sequence))
    
    def __getitem__(self, index):
        return None

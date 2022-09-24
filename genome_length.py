from Bio import SeqIO
import math

def read(fasta_file_name, bucket_length=50000):
    print("Opening FASTA file:", fasta_file_name)
    fasta_sequences = SeqIO.parse(open(fasta_file_name), 'fasta')
    total_buckets = 0
    for fasta in fasta_sequences:
        sequence_length = len(str(fasta.seq))
        print("... processing sequence", fasta.id, "with length", sequence_length)
        total_buckets += math.ceil(sequence_length / bucket_length)
    
    print("total number of buckets", total_buckets)

if __name__ == "__main__":
    read("/mnt/d/genome/Egu.v3.genome_f.fasta")
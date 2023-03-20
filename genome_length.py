from Bio import SeqIO
import math

def read(fasta_file_name, bucket_length=65536):
    print("Opening FASTA file:", fasta_file_name)
    fasta_sequences = SeqIO.parse(open(fasta_file_name), 'fasta')
    total_buckets = 0
    total_length = 0
    num_references = 0
    for fasta in fasta_sequences:
        sequence_length = len(str(fasta.seq))
        print("... processing sequence", fasta.id, "with length", sequence_length)
        total_buckets += math.ceil(sequence_length / bucket_length)
        total_length += sequence_length
        num_references += 1
    
    print("total length of sequence", total_length)
    print("total number of buckets", total_buckets)
    print("total number of references", num_references)

if __name__ == "__main__":
    read("/mnt/d/genome/GCA_000005845.2.fasta")
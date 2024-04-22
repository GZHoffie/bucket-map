
import torch
import json
import math
import numpy as np
from tqdm import tqdm
import random
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from Bio import SeqIO


# Load dictionary that maps k-mer to their corresponding index.
# A k-mer and its reverse complement are mapped to the same index.


with open("/home/zhenhao/bucket-map/seed_selection/index/9-mers.json", 'r') as dict_file:
    canonical_kmer_dict = json.load(dict_file)

# We define a utility function here that turns sequences to their 9-mer profiles.

def sequence_to_kmer_profile(sequence : str, k : int = 9):
    """
    Return the k-mer profile of the input sequence (string)
    """
    res = np.zeros(len(set(canonical_kmer_dict.values())), dtype=np.float32)
    for i in range(len(sequence) - k + 1):
        k_mer = sequence[i:i + k]
        if k_mer in canonical_kmer_dict:
            res[canonical_kmer_dict[k_mer]] = 1

    return res

def read_buckets_from_file(sequence_file_name : str, bucket_len : int, overlap_len : int):
    """
    Read the DNA sequence and store the sequences in buckets.
    """
    res = []

    for record in SeqIO.parse(sequence_file_name, "fasta"):
        record_sequence = str(record.seq)

        # Split the sequences into buckets
        num_buckets = math.ceil(len(record_sequence) / bucket_len)
        for i in range(num_buckets):
            bucket_sequence = record_sequence[i*bucket_len : (i+1)*bucket_len + overlap_len]
            res.append(bucket_sequence)

    return res

def sample_read_from_sequence(sequence : str, read_len : int):
    sample_range = max(0, len(sequence) - read_len)
    starting_pos = random.randint(0, sample_range)
    return sequence[starting_pos:starting_pos + read_len]


class DNAReadDataset(Dataset):
    def __init__(self, sequence_file, k=9, bucket_len=65536, read_len=150, num_reads_per_epoch=10000):
        """
        Dataset class to load large CS4220 sequence database. 

        Args:
            - data_file (`str`): Can either be a *.fasta file if the input is raw reads, or *.npy file
                                 if the input is k-mer profile.
            - label_df (`pd.DataFrame` or `None`): A dataframe with "labels" column indicating the label
                                                   of the data (must match with data_file), or `None` if there is
                                                   no label (in the case of test sets).
            - k (`int`): The lengt of k-mer. We use 6 in this project.
            - samples_index (`List` or `None`): list of indices of data we sample from the data file. You
                                                can use this if the dataset is very large and can't fit in memory.
                                                set this to `None` if you want to use all the data.
            - kmer_profile_on_the_fly (`bool`): If input data_file is raw reads and this set to `True`, 
                                                we will build k-mer profile on the fly. This is helpful if you want to
                                                alter the input sequences during training, or the k-mer profile can't fit in memory.
                                                Otherwise, we build k-mer profile in advance, which will speed up the
                                                training process.
            - dtype: type to store the k-mer profile. You may use, for example, `np.float32` for better precision,
                     or `np.float16` for smaller memory usage. If loaded from ".npy" file, it is always `np.float16`.
        """
        self.data_file = sequence_file
        self.bucket_sequences = read_buckets_from_file(sequence_file, bucket_len, read_len)

        # k-mer length, set to be 9.
        self.k = k
        self.read_len = read_len
        self.bucket_len = bucket_len

        self.num_reads_per_epoch = num_reads_per_epoch

    
    def __len__(self):
        return self.num_reads_per_epoch

    def __getitem__(self, idx):
        bucket = random.randrange(len(self.bucket_sequences))
        read = sample_read_from_sequence(self.bucket_sequences[bucket], self.read_len)
        read_tensor = torch.tensor(sequence_to_kmer_profile(read, self.k))
        return read_tensor, bucket


# Example usage
#input_file_path = './training_data/train_raw_reads.fasta'
input_file_path = '/home/zhenhao/mapping_data/GCA_000005845.2_ASM584v2_genomic.fna'
dataset = DNAReadDataset(input_file_path)
dataloader = DataLoader(dataset, batch_size=32, shuffle=False)


## Define the model
class BucketClassifier(nn.Module):
    # build the constructor
    def __init__(self, n_inputs, n_outputs):
        super(BucketClassifier, self).__init__()
        self.linear = torch.nn.Linear(n_inputs, n_outputs)
    # make predictions
    def forward(self, x):
        y_pred = torch.sigmoid(self.linear(x))
        return y_pred

n_inputs = len(set(canonical_kmer_dict.values()))
n_outputs = len(dataset.bucket_sequences)
classifier = BucketClassifier(n_inputs, n_outputs)

optimizer = torch.optim.Adam(classifier.parameters(), lr=0.001)
# defining Cross-Entropy loss
criterion = torch.nn.CrossEntropyLoss()
 
epochs = 50
Loss = []
acc = []
for epoch in range(epochs):
    train_loss = 0.0
    train_correct_predictions = 0
    for i, (reads, labels) in enumerate(dataloader):
        optimizer.zero_grad()
        outputs = classifier(reads)
        loss = criterion(outputs, labels)

        # Record the performance
        train_loss += loss.item()
        _, predicted = torch.max(outputs.data, 1)
        train_correct_predictions += (predicted == labels).sum()

        # Training
        loss.backward()
        optimizer.step()

    Loss.append(train_loss)
    accuracy = 100 * (train_correct_predictions.item()) / len(dataset)
    acc.append(accuracy)
    print('Epoch: {}. Loss: {}. Accuracy: {}'.format(epoch, train_loss, accuracy))
import json
import os
from tqdm import tqdm
import numpy as np
from pathlib import Path

class Seq2kmer:
    def __init__(self, k, dict_dir=None, dict_prefix=None) -> None:
        """
        Util class that convert DNA sequences to k-mer profiles.
        """
        # k-mer length
        self._k = k

        # Utils to find canonical k-mers
        self._nucleotide_to_hash = {'A': 0, 'a': 0,
                                    'C': 1, 'c': 1,
                                    'G': 2, 'g': 2,
                                    'T': 3, 't': 3}
        self._hash_to_nucleotide = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
        self._seed_mask = int("1" * 2 * self._k, 2)

        if dict_dir is None or dict_prefix is None:
            self.canonical_kmer_dict, self.canonical_kmer_list = self.init_kmer_dict()
        else:
            self.canonical_kmer_dict, self.canonical_kmer_list = self.read_dict(dict_dir, dict_prefix)


    def read_dict(self, load_directory, filename_prefix):
        """
        Read the dict file that maps kmer to their canonical k-mer.
        """
        dict_path = os.path.join(load_directory, filename_prefix + ".json")
        with open(dict_path, 'r') as json_file:
            kmer_dict = json.load(json_file)
        
        kmer_list = []
        kmer_list_path = os.path.join(load_directory, filename_prefix + ".txt")
        with open(kmer_list_path, 'r') as f:
            for line in f.readlines():
                kmer_list.append(line.rstrip())
        
        return kmer_dict, kmer_list
    

    def save_dict(self, save_directory, filename_prefix):
        """
        Save the dict to files in save_directory/filename_prefix.{json, txt}.
        """
        dict_path = os.path.join(save_directory, filename_prefix + ".json")
        with open(dict_path, "w") as json_file:
            json.dump(self.canonical_kmer_dict, json_file)


        kmer_list_path = os.path.join(save_directory, filename_prefix + ".txt")
        with open(kmer_list_path, "w") as f:
            for i in self.canonical_kmer_list:
                f.write(i + "\n")


    def _rev_comp(self, encoded_sequence):
        """
        Given a sequence of hash values returned from self._encode(sequence),
        return the reverse complement of that sequence.
        """
        rev_comp_hash = [(3 - h) for h in reversed(encoded_sequence)]
        return rev_comp_hash
    
    def _encode(self, sequence):
        """
        Convert a DNA sequence to a sequence of numbers in 0-3.
        We assume that there is no ambiguous characters such as 'N' or 'n'.
        """
        sequence_hash = [self._nucleotide_to_hash[n] for n in sequence if n in self._nucleotide_to_hash]
        return sequence_hash
    
    def _decode(self, encoded_sequence):
        """
        Convert a sequence of numbers in 0-3 to a DNA sequence.
        """
        #print(encoded_sequence)
        sequence = "".join([self._hash_to_nucleotide[n] for n in encoded_sequence if n in self._hash_to_nucleotide])
        return sequence

    
    def init_kmer_dict(self):
        """
        Initialize the canonical k-mers.

        Returns:
            - kmer_dict: A dict that maps the kmer hash to the index of the corresponding canonical k-mer.
            - kmer_list: list of canonical k-mers.
        """
        MASK = 3 # Mask to filter out single nucleotide in k-mer
        kmer_dict = {}
        kmer_list = []
        kmer_index = 0
        for i in tqdm(range(4 ** self._k), desc=f"[INFO]\t\tBuilding up kmer dict using {self._k}-mers"):
            kmer_sequence = [(i >> (j * 2)) & MASK for j in reversed(range(self._k))]
            rev_comp_sequence = self._rev_comp(kmer_sequence)

            kmer_hash = 0
            rev_comp_hash = 0
            for i in range(self._k):
                kmer_hash = kmer_hash << 2
                rev_comp_hash = rev_comp_hash << 2
                kmer_hash = kmer_hash | kmer_sequence[i]
                rev_comp_hash = rev_comp_hash | rev_comp_sequence[i]
        
            # Pick the canonical k-mer only
            canonical_kmer = kmer_sequence if kmer_hash < rev_comp_hash else rev_comp_sequence
            if self._decode(kmer_sequence) not in kmer_dict:
                # insert the kmer and its reverse complement to the dict
                kmer_dict[self._decode(kmer_sequence)] = kmer_index
                kmer_dict[self._decode(rev_comp_sequence)] = kmer_index
                # Store the canonical kmer also
                kmer_list.append(self._decode(canonical_kmer))
                kmer_index += 1
        
        # Add 'IGNORE' to indicate ambiguous k-mers
        kmer_dict['IGNORE'] = kmer_index
        kmer_list.append('IGNORE')

        print(f"[INFO]\t\tK-mer dict build complete. Dict size: {kmer_index + 1}.")
        return kmer_dict, kmer_list


    def sequence_to_kmer_profile(self, sequence):
        res = np.zeros(len(self.canonical_kmer_list))
        for i in range(len(sequence) - self._k + 1):
            k_mer = sequence[i:i + self._k]
            if k_mer in self.canonical_kmer_dict:
                res[self.canonical_kmer_dict[k_mer]] += 1
            else:
                res[-1] += 1

        res /= np.sum(res)
        return res
    
    def fasta_to_kmer_profile(self, fasta_file_path, output_file_path='default'):
        from Bio import SeqIO
        import subprocess
        
        if output_file_path == 'default':
            # output npy file that contain the k-mer profile 
            output_dir = os.path.dirname(fasta_file_path)
            output_file_name = Path(fasta_file_path).stem + f"_{self._k}mers.npy"
            output_file_path = os.path.join(output_dir, output_file_name)

        proc = subprocess.run(["wc", "-l", fasta_file_path], capture_output=True)
        num_reads = int(int(str(proc.stdout).split(" ")[0][2:]) / 2)

        res = np.zeros((num_reads, len(self.canonical_kmer_list)), dtype=np.float16)
        i = 0
        for record in tqdm(SeqIO.parse(fasta_file_path, "fasta"), desc=f"Parsing fasta file {fasta_file_path}", total=num_reads):
            sequence = str(record.seq).upper()
            kmer_profile = self.sequence_to_kmer_profile(sequence)
            res[i, :] = kmer_profile
            i += 1

        
        if output_file_path is None:
            return res
        
        elif output_file_path == 'default':
            # output npy file that contain the k-mer profile 
            output_dir = os.path.dirname(fasta_file_path)
            output_file_name = Path(fasta_file_path).stem + f"_{self._k}mers.npy"
            output_file_path = os.path.join(output_dir, output_file_name)

        np.save(output_file_path, res) 
        return res
        


if __name__ == "__main__":
    import glob

    tool = Seq2kmer(9)
    #tool.fasta_to_kmer_profile("/home/zhenhao/CS4220-project/training_data/train_raw_reads.fasta", 'default')
    #for f in glob.glob("/home/zhenhao/CS4220-project/test_data/*.fasta"):
    #    tool.fasta_to_kmer_profile(f, 'default')
    tool.save_dict("/home/zhenhao/bucket-map/seed_selection/index", "9-mers")
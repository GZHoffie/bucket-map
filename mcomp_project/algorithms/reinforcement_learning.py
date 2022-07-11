from gym import Env, spaces
from Bio import SeqIO
from mcomp_project.utils import *
import numpy as np
import math
from stable_baselines3 import DQN

# TODO: add random variance (indel/substitution) to the reference genome
class ReferenceGenome(Env):
    def __init__(self, fasta_file_name, region_length=100000, read_length=100, substitution_rate=0.02, prior=0.001):
        super(ReferenceGenome, self).__init__()

        # Parameters for genome dividing
        self.region_length = region_length
        self.read_length = read_length
        self.substitution_rate = substitution_rate

        # Read the fasta file
        fasta_sequences = SeqIO.parse(open(fasta_file_name), 'fasta')
        self.fasta_sequence = None
        for fasta in fasta_sequences:
            if self.fasta_sequence is None:
                self.fasta_sequence = fasta.seq
            else:
                self.fasta_sequence += fasta.seq

        self.sequence_length = len(self.fasta_sequence)
        self.num_chunks = math.ceil(float(self.sequence_length) / self.region_length) 

        # Define observation space
        self.observation_space = spaces.MultiDiscrete([4] * self.read_length)
        
        # Define the action space
        self.action_space = spaces.Discrete(self.num_chunks)

        # Store the answer to the last observation
        self.last_observation_region = None
    
    def step(self, action):
        reward = 1 if self.last_observation_region == action else 0

        index = np.random.randint(0, self.sequence_length - self.read_length - 1)
        self.last_observation_region = index // self.region_length
        observation = np.array([char_to_index_map[c] for c in self.fasta_sequence[index:index + self.read_length]])
        
        return observation, reward, True, {}
    
    def reset(self):
        index = np.random.randint(0, self.sequence_length - self.read_length - 1)
        self.last_observation_region = index // self.region_length
        observation = np.array([char_to_index_map[c] for c in self.fasta_sequence[index:index + self.read_length]])
        return observation


if __name__ == "__main__":
    env = ReferenceGenome("/home/zhenhao/mcomp-dissertation/sequence_sample.fasta")
    model = DQN("MlpPolicy", env, verbose=1)
    model.learn(total_timesteps=1000000, log_interval=4)
    model.save("mapper")
import math
import random
import numpy as np
from matplotlib import pyplot as plt

class SimulateKMerFrequency:
    """
    Calculate a theoretical time/space usage and performance analysis (correctness)
    given certain parameters set.

    TODO: Simulate the case for error in reads.
    TODO: Try negative binomial distribution instead of Poisson distribution.
    """
    def __init__(self, genome_size, bucket_num, sample_num, k, l, r=100, prior=0.01) -> None:
        """
        Args:
            genome_size: length of the reference genome.
            bucket_num: number of buckets we divide the reference genome into.
            k: the effective neucleotides of the seed we are using 
                (e.g. for seed 11011011000011011, we have k=10)
            l: the length of the seed
                (e.g. for seed 11011011000011011, we have k=17)
            r: length of each short read
        """
        self.genome_size = genome_size
        self.bucket_num = bucket_num
        self.sample_num = sample_num
        self.bucket_size = math.ceil(self.genome_size / self.bucket_num)
        print("Expected bucket size:", self.bucket_size)
        self.k = k
        self.l = l
        self.r = r
        self.prior = prior

        # Calculate parameters for the distribution
        self.lam = (self.bucket_size - self.l) / (4 ** self.k)
        self.hit_lam = 1 + (self.bucket_size - self.r) / (4 ** self.k)
        print("Expectation for hit bucket:", self.hit_lam)
        print("Expectation for other buckets:", self.lam)


    def simulate_frequency_distribution(self, simulate_num=1000, plot=False):
        """
        Simulate the calculation of product of probabilities and
        return the samples.
        """
        largest_samples = np.zeros(simulate_num)
        for i in range(simulate_num):
            largest_num = float('-inf')
            samples = np.random.poisson(lam=self.lam, size=self.sample_num * (self.bucket_num - 1)) + self.prior
            samples = np.reshape(samples, (self.bucket_num - 1, self.sample_num))
            samples = np.sum(np.log(samples), axis=1)
            largest_samples[i] = np.max(samples)
        #print(largest_samples)
        if plot:
            plt.hist(largest_samples, bins='sturges')
            plt.show()
        return largest_samples


    
    def simulate_correctness(self, simulate_num=1000, plot=False):
        # expected count from the correct bucket
        expectation = np.log(self.hit_lam + self.prior) * self.sample_num - 0.5
        print("Expected sum of log probability:", expectation)

        # calculate a lower bound
        lower_bound = np.log(1 + self.prior) * self.sample_num

        # start simulation
        samples = self.simulate_frequency_distribution(simulate_num, plot)
        print("Simulated correctness:", np.sum(samples < expectation) / simulate_num)
        print("Worst case correctness:", np.sum(samples < lower_bound) / simulate_num)

if __name__ == "__main__":
    sim = SimulateKMerFrequency(7000000, 700, 50, 7, 10)
    sim.simulate_correctness(1000, plot=True)

            

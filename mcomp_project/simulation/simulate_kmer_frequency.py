import math
import random
import numpy as np

class SimulateKMerFrequency:
    """
    Calculate a theoretical time/space usage and performance analysis (correctness)
    given certain parameters set.

    TODO: Simulate the case for error in reads.
    TODO: Try negative binomial distribution instead of Poisson distribution.
    """
    def __init__(self, genome_size, bucket_num, sample_num, k, l, prior=0.01) -> None:
        """
        Args:
            genome_size: length of the reference genome.
            bucket_num: number of buckets we divide the reference genome into.
            k: the effective neucleotides of the seed we are using 
                (e.g. for seed 11011011000011011, we have k=10)
            l: the length of the seed
                (e.g. for seed 11011011000011011, we have k=17)
        """
        self.genome_size = genome_size
        self.bucket_num = bucket_num
        self.sample_num = sample_num
        self.bucket_size = math.ceil(self.genome_size / self.bucket_num)
        print("Expected bucket size:", self.bucket_size)
        self.k = k
        self.l = l
        self.prior = prior

        # Calculate parameters for the distribution
        self.lam = (self.bucket_size - self.l) / (4 ** self.k)
        self.hit_lam = 1 + (self.bucket_size - self.l - self.sample_num) / (4 ** self.k)


    def simulate_frequency_distribution(self, simulate_num=1000):
        """
        Simulate the calculation of product of probabilities and
        return the samples.
        """
        largest_samples = np.zeros(simulate_num)
        for i in range(simulate_num):
            largest_num = float('-inf')
            for _ in range(self.bucket_num - 1):
                res = np.sum(np.log((np.random.poisson(lam=self.lam, size=self.sample_num) + self.prior)))
                if res > largest_num:
                    largest_num = res
            largest_samples[i] = largest_num
        print(largest_samples)
        return largest_samples


    
    def simulate_correctness(self, simulate_num=1000):
        # expected count from the correct bucket
        expectation = np.log(self.hit_lam) * self.sample_num - 0.5
        print("Expected sum of log probability:", expectation)
        # start simulation
        samples = self.simulate_frequency_distribution(simulate_num)
        print("Simulated correctness:", np.sum(samples < expectation) / simulate_num)

if __name__ == "__main__":
    sim = SimulateKMerFrequency(30000000, 3000, 50, 6, 10)
    sim.simulate_correctness(100)

            

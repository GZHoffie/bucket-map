import sys
sys.path.append("/home/zhenhao/bucket-map/")

import numpy as np

from seed_selection.dataset import sequence_to_kmer_profile, read_buckets_from_file

class RepetitiveRegionFilter:
    def __init__(self, bucket_len, overlap_len, k) -> None:
        self.bucket_len = bucket_len
        self.overlap_len = overlap_len
        self.k = k
    

    def read(self, sequence_file):
        buckets = read_buckets_from_file(sequence_file, self.bucket_len, self.overlap_len)
        res = []
        for i in range(len(buckets)):
            buckets_profile = sequence_to_kmer_profile(buckets[i], self.k)
            res.append(buckets_profile)
        
        return res
    
    def ji_matrix(self, bucket_profile):
        res = np.zeros((len(bucket_profile), len(bucket_profile)))
        for i in range(len(bucket_profile)):
            for j in range(i + 1, len(bucket_profile)):
                jaccard_index = np.sum(bucket_profile[i] * bucket_profile[j]) / np.sum(np.logical_or(bucket_profile[i], bucket_profile[j]))
                res[i][j] = res[j][i] = jaccard_index
        
        return res
    
if __name__ == "__main__":
    print("yes")
    filter = RepetitiveRegionFilter(65536, 150, 9)
    profile = filter.read('/home/zhenhao/mapping_data/GCA_000005845.2_ASM584v2_genomic.fna')
    matrix = filter.ji_matrix(profile)
    print(len(matrix))
    print(np.max(matrix))
    #print(filter.ji_matrix(profile))



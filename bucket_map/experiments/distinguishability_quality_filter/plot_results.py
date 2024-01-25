import pandas as pd
import matplotlib.pyplot as plt


# Benchmarking items
items = ["Average number of buckets an original read is mapped to", 
         "Average number of buckets a reverse complement of the read is mapped to", 
         "Number of Q-grams with distinguishability"]
data = {"Number of Candidate Buckets", "Remaining k-mers"}


# directory of the log files
log_directory = "./bucket_map/experiments/distinguishability_quality_filter/log/"

for distinguishability in ["0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"]:
    with open(log_directory + "bucketmap_" + distinguishability + "_0_map.log", "r") as index_log:
        orig_candidate_buckets = 0
        rev_candidate_buckets = 0
        remaining_k_mers = 0


        for line in index_log:
            if items[0] in line:
                orig_candidate_buckets = float((line.split(' ')[-1]).strip()[:-1])
            elif items[1] in line:
                rev_candidate_buckets = float((line.split(' ')[-1]).strip()[:-1])
            elif items[2] in line:
                remaining_k_mers = int((line.split(' ')[-2]).strip())
            
        print(orig_candidate_buckets, rev_candidate_buckets, remaining_k_mers)


        
print(pd.DataFrame(data))

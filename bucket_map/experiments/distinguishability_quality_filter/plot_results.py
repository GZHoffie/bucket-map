import pandas as pd
import matplotlib.pyplot as plt

# Named of the tools used for benchmarking
tools = ["bowtie2", "bwa", "minimap2", "subread", "bucketmap"]

# Benchmarking items
items = ["User time (seconds)", "Maximum resident set size (kbytes)", "Major (requiring I/O) page faults", "Minor (reclaiming a frame) page faults"]
data = {"tools": tools}
for item in items:
    data[item] = []

# directory of the log files
log_directory = "./bucket_map/benchmark/short_read/log/"

for tool in tools:
    with open(log_directory + tool + "_map.time", "r") as index_log:
        for line in index_log:
            for item in items:
                if item in line:
                    data[item].append(float((line.split(' ')[-1]).strip()))
                    break
        
print(pd.DataFrame(data))

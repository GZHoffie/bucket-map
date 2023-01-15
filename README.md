![](./bucketmap-logo.png)

# *BucketMap*: Hierarchical DNA Read Mapping

**This repository contains my work for the MComp Dissertation at School of Computing, National University of Singapore. This work is supervised by [Prof. Wong Limsoon](https://www.comp.nus.edu.sg/~wongls/).**

*BucketMap* is a novel DNA mapping tool that is extremely light-weight, memory-saving and fast. For a 1.7 billion base-pair long reference genome, and 1 million short reads of length 300, *BucketMap* is able to complete the mapping within 8 minutes, with 97% accuracy and a peak memory of only 0.8 GB.

Named after the classic *Bucket Sort* algorithm, the idea of *BucketMap* is simple.

- We first divide the reference genome into overlapping buckets,
- For each read, we use samping and bit-parallel data structures to quickly find the candidate buckets that might contain the read,
- Within each bucket, we go through all reads that are mapped to it and find their exact locations.

## Installation

This tool is built with [Seqan3](https://docs.seqan.de/seqan/3-master-user/index.html). To properly build the package, you need to have GCC >= 11.3, G++ and CMake installed.

```bash
git clone https://github.com/GZHoffie/bucket-map.git
cd bucket-map/bucket_map
```
## DNA read simulator

In `tools/short_read_simulator.h`, we implement a more realistic generator of DNA reads that is used for our benchmarking. The tool is easy to use,

```C++
#include "tools/short_read_simulator.h"

int main() {
    // path to the fasta reference genome file
    std::filesystem::path genome_file = "/mnt/d/genome/Egu.v3.genome_f.fasta";

    // path to store the generated fastq file
    std::filesystem::path output_path = "/mnt/d/genome/test";


    // length of buckets, used specifically for accuracy benchmark for BucketMap
    int bucket_length = 65536;

    // length of each read
    int read_length = 300;

    // total number of reads
    int num_reads = 1000000;

    // setting error rates 
    float substitution_rate = 0.002;
    float insertion_rate = 0.00025;
    float deletion_rate = 0.00025;

    short_read_simulator sim(bucket_length, read_length, 0.002, 0.00025, 0.00025);
    sim.read(genome_file);
    sim.generate_fastq_file(output_path, "sim_illumina_1M", num_reads);
}
```

The number of errors (substitutions and indels) are randomly generated with a Poisson distribution with mean being `error_rate * read_length`. Upon running the generator, 3 files will be generated in the `output_path`.

1. `sim_illumina_1M.fastq`: the query file containing all the simulated reads.
2. `sim_illumina_1M.bucket_ground_truth`: specific to benchmarking the accuracy of *BucketMap*, where each line contains the `bucket_id` and `offset` the read comes from, as well as the ground truth CIGAR string.
2. `sim_illumina_1M.bucket_ground_truth`: each line contains the `reference_id` and `offset` the read comes from, as well as the ground truth CIGAR string.
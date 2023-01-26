<p align="center">
  <img src="./bucketmap-logo.png"/>
</p>

# *BucketMap*: Hierarchical DNA Read Mapping

**This repository contains my experimental work for the MComp Dissertation at School of Computing, National University of Singapore. This work is supervised by [Prof. Wong Limsoon](https://www.comp.nus.edu.sg/~wongls/).**

*BucketMap* is a novel DNA mapping tool that is extremely light-weight, memory-saving and fast. For a 1.7 billion base-pair long reference genome, and 1 million short reads of length 300, *BucketMap* is able to complete the mapping within 8 minutes, with 97% accuracy and a peak memory of only 0.8 GB.

Named after the classic *Bucket Sort* algorithm, the idea of *BucketMap* is simple.

- We first divide the reference genome into overlapping buckets,
- For each read, we use samping and bit-parallel data structures to quickly find the candidate buckets that might contain the read,
- Within each bucket, we go through all reads that are mapped to it and find their exact locations.

## Installation

This tool is built with [Seqan3](https://docs.seqan.de/seqan/3-master-user/index.html). To properly build the package, you need to have GCC >= 11.3, G++ and CMake installed.

```bash
git clone https://github.com/GZHoffie/bucket-map.git
cd bucket-map
mkdir build
cd build

# Set the path to the fasta file and the desired size of bucket.
FASTA_FILE=/mnt/d/genome/Egu.v3.genome_f.fasta # Should be an absolute path to the genome file.
BUCKET_LEN=65536

# Build the project with the following command
cmake ../bucket_map/ -DBM_FASTA_FILE=${FASTA_FILE} -DBM_BUCKET_LEN=${BUCKET_LEN}
```
Note that if the reference genome and bucket length is to be changed, we need to rebuild the project by rerunning the above commands.


Now, we can build the binary file using CMake. For the alignment-free version of BucketMap, which doesn't output the CIGAR string in the output SAM file, we can use

```bash
cmake --build . --target bucketmap
```

Or we may build the full version of BucketMap, which further uses Smith-Waterman algorithm to verify the mapping,

```bash
cmake --build . --target bucketmap_align
```

and the binary files named `bucketmap` and `bucketmap_align` will be built under the `./build` directory.

## Usage

The usage of BucketMap is simple. To build the index files, you may use

```bash
<path_to_bucketmap_binary>/bucketmap -x -i <index_name>
```

which will output `<index_name>.qgram` and `<index_name>.bucket_id` under the current directory.

To map the reads in a fastq file, you may run the following command (under the directory of the index file if you have built the index file using the above command)

```bash
<path_to_bucketmap_binary>/bucketmap -i <index_name> -q <fastq_file_path> -o <output_sam_file_path>
```

which will output the `.sam` file in the desired location. You may also run the second command directly, which will build the index and do the mapping in one go. To perform pairwise alignment and output the CIGAR string in the sam files, you may simply replace `bucketmap` with `bucketmap_align`.



## Benchmarking

We benchmark our read mapper against several popular DNA read mappers available, including

- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), 
- [BWA-MEM](https://bio-bwa.sourceforge.net/), 
- [subread](https://github.com/ShiLab-Bioinformatics/subread).

We document a detailed guide on how we do the benchmarking, which metrics we are using, as well as some preliminary benchmark results, in [this page](./bucket_map/benchmark/README.md).

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

    short_read_simulator sim(bucket_length, 
                             read_length, 
                             substitution_rate, 
                             insertion_rate, 
                             deletion_rate
                             );
    sim.read(genome_file);
    sim.generate_fastq_file(output_path, "sim_illumina_1M", num_reads);
}
```

The number of errors (substitutions and indels) are randomly generated with a Poisson distribution with mean being `error_rate * read_length`. Upon running the generator, 3 files will be generated in the `output_path`.

1. `sim_illumina_1M.fastq`: the query file containing all the simulated reads.
2. `sim_illumina_1M.bucket_ground_truth`: specific to benchmarking the accuracy of *BucketMap*, where each line contains the `bucket_id` and `offset` the read comes from, as well as the ground truth CIGAR string.
3. `sim_illumina_1M.bucket_ground_truth`: each line contains the `reference_id` and `offset` the read comes from, as well as the ground truth CIGAR string.
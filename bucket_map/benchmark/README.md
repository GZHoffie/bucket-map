# Benchmarking of DNA Read Mapping Tools

In this page, we note down the benchmarking methods and results of some DNA read mapping tools. The tools we use for benchmarking include [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), [BWA-MEM](https://bio-bwa.sourceforge.net/), [subread](https://github.com/ShiLab-Bioinformatics/subread).

## Installation of Tools

We write down how we download the read mapping tools for future reference. We also add script of adding the binary applications to `PATH` for easier usage. All tools can be built from source

- **Bowtie2** can be downloaded directly using [bioconda](https://anaconda.org/bioconda).
  ```bash
  conda install -c bioconda bowtie2
  ```
  
  or from source
  ```bash
  git clone https://github.com/BenLangmead/bowtie2.git
  cd bowtie2; make
  echo "export PATH=\${PATH}:$(pwd)" >> ~/.bashrc
  ```

- **BWA-MEM**.
  
  ```bash
  git clone https://github.com/lh3/bwa.git
  cd bwa; make
  echo "export PATH=\${PATH}:$(pwd)" >> ~/.bashrc
  ```

- **Subread**.

  ```bash
  git clone https://github.com/ShiLab-Bioinformatics/subread.git
  cd subread/src
  make -f Makefile.Linux
  cd ..
  echo "export PATH=\${PATH}:$(pwd)/bin" >> ~/.bashrc
  ```

 - **Minimap2**.

   ```bash
   git clone https://github.com/lh3/minimap2
   cd minimap2 && make
   echo "export PATH=\${PATH}:$(pwd)" >> ~/.bashrc
   ```

After adding everything to `.bashrc`, we can run
```bash
source ~/.bashrc
```
to use the binaries by their names directly.

## Dataset generation

```bash
python3 gen_reads.py -r /home/guzh/data/mapping/GCA_022991725.1_PDT001286877.1_genomic.fna -R 150 -o simulated_data --bam --vcf

# Convert BAM file to SAM file
samtools index simulated_data_golden.bam
samtools view -h -o simulated_data_golden.sam simulated_data_golden.bam

# Decompress fq.gz file
zcat simulated_data_read1.fq.gz >> simulated_data_read1.fq
```

## Benchmarking Time & Memory Usage
To benchmark the time and memory usage, we mainly use Linux's [`time`](https://man7.org/linux/man-pages/man1/time.1.html) command. An example usage is

```bash
/usr/bin/time -v bowtie2 -x "${INDEX_INDICATOR}_bowtie" -U ${FASTQ_PATH} -S "${BENCHMARK_PATH}/output/bowtie2_map.sam" &> "${BENCHMARK_PATH}/log/bowtie2_map.log"
```

Note that we write `/usr/bin/time` instead of simply `time` to specify the binary we are using. An output is shown below.

```log
	Command being timed: "bowtie2 -x egu_bowtie -U /mnt/d/genome/test/sim_illumina_1M.fastq -S /home/zhenhao/bucket-map/bucket_map/benchmark/output/bowtie2_map.sam"
	User time (seconds): 852.35
	System time (seconds): 1.33
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 14:21.50
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 1968900
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 102
	Minor (reclaiming a frame) page faults: 8740
	Voluntary context switches: 97183
	Involuntary context switches: 636
	Swaps: 0
	File system inputs: 3866056
	File system outputs: 1382648
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

We mainly look at the following fields:

- `User time (seconds)`: the total time used by the program.
- `Maximum resident set size (kbytes)`: peak memory usage of the program.
- `Major (requiring I/O) page faults`, `Minor (reclaiming a frame) page faults`: number of page faults caused by the program.

We write down a script that gathers the output of the mappers and the returned log of `time` command in [`benchmark_index.sh`](./benchmark_index.sh) and [`benchmark_map.sh`](./benchmark_map.sh). The generated logs and outputs of programs are stored under [`./log`](./log) and the output sam files under [`./output`](./output).

## Benchmarking Accuracy

To benchmark accuracy: whether the mapper returns the true location of the read, we use the following metrics.

- Percentage of reads mapped,
- Average number of locations a read is mapped to,
- Percantage of reads mapped to the correct reference (chromosome or contig),
- Percentage of reads mapped to the correct position.

The metrics can be checked directly using the sam file. The benchmarking script is written in [`sam_file_analyzer.cpp`](./sam_file_analyzer.cpp), which checks all sam files under `./output`.

## Benchmarking Results

We primarily use the reference genome of *Elaeis guineensis* (to be checked) with length 1,701,312,507. The dataset consists of 16 chromosomes and 916 contigs. For all tools, we simply use the default settings.

### Indexing

We record the time and memory performance when producing the index files. Usually this is less important than performance of mapping, as it only requires a one-time run.

|Tool|User time (s)|Peak memory usage (Kb)|Major page faults|Minor page faults|
|----|-------------|----------------------|-----------------|-----------------|
|Bowtie2|2424.85|7648524|832483|7510117|
|BWA-MEM|1529.78|2495668|13|424721|
|Subread|633.29|3197744|8|1347347|

*Note*: The indexing of Bowtie2 is very memory consuming. We need to make sure that we have enough memory before running the indexing (reboot if possible before running), otherwise it will run into exit status `247`.

### Mapping

The time and memory usage of the tools with 1 million reads of length 300 (with sequencing error rates similar to Illumina) are shown below.

|Tool|User time (s)|Peak memory usage (Kb)|Major page faults|Minor page faults|
|----|-------------|----------------------|-----------------|-----------------|
|Bowtie2|852.35|1968900|102|8740|
|BWA-MEM|563.19|2995356|25|436960|
|Subread|573.93|4201588|13|734950|
|**BucketMap**|426.78|865568|0|138474|
|**Alignment-free BucketMap**|320.95|865344|0|135648|

And the accuracy of the tools when comparing with the ground truth is

|Tool|% reads mapped|Average mapped locations|% correct reference|% correct positions|
|----|--------------|------------------------|-------------------|-------------------|
|Bowtie2|99.8904%|1|98.4358%|97.6452%|
|BWA-MEM|100%|1.00142|98.9025%|98.0892%|
|Subread|96.7519%|1|95.735%|94.8024%|
|**BucketMap**|98.6348%|1.14538|97.8566%|97.3444%|
|**Alignment-free BucketMap**|98.6348%|1.14538|97.8566%|97.3444%|
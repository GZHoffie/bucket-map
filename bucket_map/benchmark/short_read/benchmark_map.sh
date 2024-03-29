#!/bin/bash


# the path to the fasta/fastq file is required to be an absolute path
FASTQ_PATH=$1
BENCHMARK_PATH=$2
INDEX_INDICATOR=$3


# go to the index path
cd "${BENCHMARK_PATH}/../index"

# run bowtie2 to map the reads
echo "Mapping using bowtie2"
#/usr/bin/time -o "${BENCHMARK_PATH}/log/bowtie2_map.time" -v bowtie2 -x "${INDEX_INDICATOR}_bowtie" -U ${FASTQ_PATH} -S "${BENCHMARK_PATH}/output/bowtie2_map.sam" &> "${BENCHMARK_PATH}/log/bowtie2_map.log"

# run bwa to map the reads
echo "Mapping using bwa"
#/usr/bin/time -o "${BENCHMARK_PATH}/log/bwa_map.time" -v bwa mem "${INDEX_INDICATOR}_bwa" ${FASTQ_PATH} > "${BENCHMARK_PATH}/output/bwa_map.sam" 2> "${BENCHMARK_PATH}/log/bwa_map.log"

# run subread to map the reads
echo "Mapping using subread"
#/usr/bin/time -o "${BENCHMARK_PATH}/log/subread_map.time" -v subread-align --SAMoutput -i "${INDEX_INDICATOR}_subread" -r ${FASTQ_PATH} -t 1 -o "${BENCHMARK_PATH}/output/subread_map.sam" &> "${BENCHMARK_PATH}/log/subread_map.log"

# run minimap2
echo "Mapping using minimap2"
#/usr/bin/time -o "${BENCHMARK_PATH}/log/minimap2_map.time" -v minimap2 -a "${INDEX_INDICATOR}_minimap.mmi" ${FASTQ_PATH} > "${BENCHMARK_PATH}/output/minimap2_map.sam" 2> "${BENCHMARK_PATH}/log/minimap2_map.log"

# run bucketmap
echo "Mapping using BucketMap"
/usr/bin/time -o "${BENCHMARK_PATH}/log/bucketmap_map.time" -v ../../../build/bucketmap --version-check 0 -r 300 -s 20 -e 0.6 -l 14 -b 10 -i "${INDEX_INDICATOR}_bucketmap" -q ${FASTQ_PATH} -o "${BENCHMARK_PATH}/output/bucketmap_fracMinHash_map.sam" &> "${BENCHMARK_PATH}/log/bucketmap_fracMinHash_map.log"

echo "Mapping using BucketMap_align"
#/usr/bin/time -o "${BENCHMARK_PATH}/log/bucketmap_align_map.time" -v bucketmap_align --version-check 0 -r 160 -s 30 -e 0.6 -u 30 -i "${INDEX_INDICATOR}_bucketmap" -q ${FASTQ_PATH} -o "${BENCHMARK_PATH}/output/bucketmap_align_map.sam" &> "${BENCHMARK_PATH}/log/bucketmap_align_map.log"
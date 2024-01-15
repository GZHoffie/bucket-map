#!/bin/bash


# the path to the fasta/fastq file is required to be an absolute path
FASTQ_PATH=$1
BENCHMARK_PATH=$2
INDEX_INDICATOR=$3
FASTA_PATH=$4


# go to the index path
cd "${BENCHMARK_PATH}/index"

# run minimap2
echo "Mapping using minimap2"
#/usr/bin/time -o "${BENCHMARK_PATH}/log/minimap2_map.time" -v minimap2 -a "${INDEX_INDICATOR}_minimap.mmi" ${FASTQ_PATH} > "${BENCHMARK_PATH}/output/minimap2_map.sam" 2> "${BENCHMARK_PATH}/log/minimap2_map.log"

# run GraphMap
echo "Mapping using GraphMap"
/usr/bin/time -o "${BENCHMARK_PATH}/log/graphmap_map.time" -v graphmap align -v 0 -r ${FASTA_PATH} -d ${FASTQ_PATH} -o "${BENCHMARK_PATH}/output/graphmap_map.sam" &> "${BENCHMARK_PATH}/log/graphmap_map.log"


# run bucketmap
echo "Mapping using BucketMap"
/usr/bin/time -o "${BENCHMARK_PATH}/log/bucketmap_map.time" -v ../../../build/bucketmap --version-check 0 -s 30 -e 0.9 -n 0.1 -l 12 -p 20 -u 5 -i "${INDEX_INDICATOR}_bucketmap" -q ${FASTQ_PATH} -o "${BENCHMARK_PATH}/output/bucketmap_map.sam" &> "${BENCHMARK_PATH}/log/bucketmap_map.log"

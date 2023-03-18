#!/bin/bash

FASTA_PATH=$1 # the path to the fasta file is required to be an absolute path
BENCHMARK_PATH=$2
INDEX_INDICATOR=$3

# go to the index directory
cd "${BENCHMARK_PATH}/index"

# run the indexing method of bowtie2
echo "Indexing using bowtie2"
/usr/bin/time -o "${BENCHMARK_PATH}/log/bowtie2_index.time" -v bowtie2-build ${FASTA_PATH} "${INDEX_INDICATOR}_bowtie" &> "${BENCHMARK_PATH}/log/bowtie2_index.log"

# run indexing for bwa
echo "Indexing using bwa"
#/usr/bin/time -o "${BENCHMARK_PATH}/log/bwa_index.time" -v bwa index -p "${INDEX_INDICATOR}_bwa" ${FASTA_PATH} &> "${BENCHMARK_PATH}/log/bwa_index.log"

# run indexing for subread
echo "Indexing using subread"
#/usr/bin/time -o "${BENCHMARK_PATH}/log/subread_index.time" -v subread-buildindex -o "${INDEX_INDICATOR}_subread" ${FASTA_PATH} &> "${BENCHMARK_PATH}/log/subread_index.log"

# run indexing for minimap2
echo "Indexing using minimap2"
#/usr/bin/time -o "${BENCHMARK_PATH}/log/minimap2_index.time" -v minimap2 -d "${INDEX_INDICATOR}_minimap.mmi" ${FASTA_PATH} &> "${BENCHMARK_PATH}/log/minimap2_index.log"

# run indexing for BucketMap
echo "Indexing using BucketMap"
#/usr/bin/time -o "${BENCHMARK_PATH}/log/bucketmap_index.time" -v bucketmap -x -i "${INDEX_INDICATOR}_bucketmap" &> "${BENCHMARK_PATH}/log/bucketmap_index.log"

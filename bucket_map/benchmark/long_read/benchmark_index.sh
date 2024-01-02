#!/bin/bash

FASTA_PATH=$1 # the path to the fasta file is required to be an absolute path
BENCHMARK_PATH=$2
INDEX_INDICATOR=$3

# go to the index directory
cd "${BENCHMARK_PATH}/index"


# run indexing for minimap2
echo "Indexing using minimap2"
/usr/bin/time -o "${BENCHMARK_PATH}/log/minimap2_index.time" -v minimap2 -x map-ont -d "${INDEX_INDICATOR}_minimap.mmi" ${FASTA_PATH} &> "${BENCHMARK_PATH}/log/minimap2_index.log"

# run indexing for BucketMap
echo "Indexing using BucketMap"
/usr/bin/time -o "${BENCHMARK_PATH}/log/bucketmap_index.time" -v bucketmap -x -i "${INDEX_INDICATOR}_bucketmap" &> "${BENCHMARK_PATH}/log/bucketmap_index.log"

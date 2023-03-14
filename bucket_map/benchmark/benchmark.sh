#!/bin/bash

GENOME_FILE="/mnt/d/genome/GRCh38_latest_genomic.fasta"
BENCHMARK_PATH="/home/zhenhao/bucket-map/bucket_map/benchmark"
QUERY_FILE="/mnt/d/genome/test/GRCH38_300_1M.fastq"
INDICATOR="GRCh38"

# initialize benchmark directory if it doesnt exist
mkdir -p ${BENCHMARK_PATH}
mkdir -p "${BENCHMARK_PATH}/index"
mkdir -p "${BENCHMARK_PATH}/log"
mkdir -p "${BENCHMARK_PATH}/output"

# Run index benchmarking
#./benchmark_index.sh ${GENOME_FILE} ${BENCHMARK_PATH} ${INDICATOR}

# Run map benchmarking
./benchmark_map.sh ${QUERY_FILE} ${BENCHMARK_PATH} ${INDICATOR}
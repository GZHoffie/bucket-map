#!/bin/bash

GENOME_FILE="/mnt/d/genome/Egu.v3.genome_f.fasta"
BENCHMARK_PATH="/home/zhenhao/bucket-map/bucket_map/benchmark"
QUERY_FILE="/mnt/d/genome/TS1.81.90.001.fq"
INDICATOR="egu"

# initialize benchmark directory if it doesnt exist
mkdir -p ${BENCHMARK_PATH}
mkdir -p "${BENCHMARK_PATH}/index"
mkdir -p "${BENCHMARK_PATH}/log"
mkdir -p "${BENCHMARK_PATH}/output"

# Run index benchmarking
#./benchmark_index.sh ${GENOME_FILE} ${BENCHMARK_PATH} ${INDICATOR}

# Run map benchmarking
./benchmark_map.sh ${QUERY_FILE} ${BENCHMARK_PATH} ${INDICATOR}

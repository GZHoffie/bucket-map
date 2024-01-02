#!/bin/bash

GENOME_FILE="/home/zhenhao/data/mapping/GCF_001890245.1_ASM189024v1_genomic.fna"
BENCHMARK_PATH="/home/zhenhao/bucket-map/bucket_map/benchmark"
QUERY_FILE="/home/zhenhao/data/mapping/ecoli_simulated_read1_renamed.fq"
INDICATOR="EColi"

# initialize benchmark directory if it doesnt exist
mkdir -p ${BENCHMARK_PATH}
mkdir -p "${BENCHMARK_PATH}/index"
mkdir -p "${BENCHMARK_PATH}/log"
mkdir -p "${BENCHMARK_PATH}/output"

# Run index benchmarking
#./benchmark_index.sh ${GENOME_FILE} ${BENCHMARK_PATH} ${INDICATOR}

# Run map benchmarking
./benchmark_map.sh ${QUERY_FILE} ${BENCHMARK_PATH} ${INDICATOR}

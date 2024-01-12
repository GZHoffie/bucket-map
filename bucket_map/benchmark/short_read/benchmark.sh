#!/bin/bash

GENOME_FILE="/home/zhenhao/mapping_data/GCA_004358405.1_ASM435840v1_genomic.fna"
BENCHMARK_PATH="/home/zhenhao/bucket-map/bucket_map/benchmark/short_read"
QUERY_FILE="/home/zhenhao/mapping_data/EColi_sim.bwa.read1.fastq"
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

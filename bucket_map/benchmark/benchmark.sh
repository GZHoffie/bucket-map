#!/bin/bash

GENOME_FILE="/home/guzh/data/mapping/GCA_022991725.1_PDT001286877.1_genomic.fna"
BENCHMARK_PATH="/home/guzh/bucket-map/bucket_map/benchmark"
QUERY_FILE="/home/guzh/neat-genreads/simulated_data_read1.fq"
INDICATOR="EColi"

# initialize benchmark directory if it doesnt exist
mkdir -p ${BENCHMARK_PATH}
mkdir -p "${BENCHMARK_PATH}/index"
mkdir -p "${BENCHMARK_PATH}/log"
mkdir -p "${BENCHMARK_PATH}/output"

# Run index benchmarking
./benchmark_index.sh ${GENOME_FILE} ${BENCHMARK_PATH} ${INDICATOR}

# Run map benchmarking
./benchmark_map.sh ${QUERY_FILE} ${BENCHMARK_PATH} ${INDICATOR}

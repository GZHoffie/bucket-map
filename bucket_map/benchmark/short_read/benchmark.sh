#!/bin/bash

HOME_DIR="/home/zhenhao"
GENOME_FILE="${HOME_DIR}/mapping_data/GRCh38_adjusted.fna"
BENCHMARK_PATH="${HOME_DIR}/bucket-map/bucket_map/benchmark/short_read"
QUERY_FILE="${HOME_DIR}/mapping_data/GRCh38_sim.bwa.read1.fastq"
INDICATOR="GRCh38_fracMinHash"
BUCKET_LEN=65536

BUILD_BUCKETMAP=1

# Build bucketmap
if [ $BUILD_BUCKETMAP -eq 1 ]
then
    mkdir -p ../../../build/
    cd ../../../build/
    cmake ../bucket_map/ -DBM_FASTA_FILE=${GENOME_FILE} -DBM_BUCKET_LEN=${BUCKET_LEN}
    cmake --build . --target bucketmap
    cd ../bucket_map/benchmark/short_read
fi


# initialize benchmark directory if it doesnt exist
mkdir -p ${BENCHMARK_PATH}
mkdir -p "${BENCHMARK_PATH}/../index"
mkdir -p "${BENCHMARK_PATH}/log"
mkdir -p "${BENCHMARK_PATH}/output"

# Run index benchmarking
./benchmark_index.sh ${GENOME_FILE} ${BENCHMARK_PATH} ${INDICATOR}

# Run map benchmarking
./benchmark_map.sh ${QUERY_FILE} ${BENCHMARK_PATH} ${INDICATOR}

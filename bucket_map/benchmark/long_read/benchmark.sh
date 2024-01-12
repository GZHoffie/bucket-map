#!/bin/bash
HOME_DIR="/home/guzh"
GENOME_FILE="${HOME_DIR}/data/mapping/GCA_022991725.1_PDT001286877.1_genomic.fna"
BENCHMARK_PATH="${HOME_DIR}/bucket-map/bucket_map/benchmark/long_read"
QUERY_FILE="${HOME_DIR}/data/mapping/lr_simulated/sd_0001.fastq"
INDICATOR="EColi"
BUCKET_LEN=65536

# Build bucketmap
cd ../../../build/
cmake ../bucket_map/ -DBM_FASTA_FILE=${FASTA_FILE} -DBM_BUCKET_LEN=${BUCKET_LEN}
cmake --build . --target bucketmap


# initialize benchmark directory if it doesnt exist
mkdir -p ${BENCHMARK_PATH}
mkdir -p "${BENCHMARK_PATH}/index"
mkdir -p "${BENCHMARK_PATH}/log"
mkdir -p "${BENCHMARK_PATH}/output"

# Run index benchmarking
#./benchmark_index.sh ${GENOME_FILE} ${BENCHMARK_PATH} ${INDICATOR}

# Run map benchmarking
./benchmark_map.sh ${QUERY_FILE} ${BENCHMARK_PATH} ${INDICATOR} ${GENOME_FILE}

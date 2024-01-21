#!/bin/bash

HOME_DIR="/home/zhenhao"
GENOME_FILE="${HOME_DIR}/mapping_data/GRCh38_adjusted.fna"
BENCHMARK_PATH="${HOME_DIR}/bucket-map/bucket_map/experiments/distinguishability_quality_filter"
QUERY_FILE="${HOME_DIR}/mapping_data/GRCh38_sim_short.bwa.read1.fastq"
INDICATOR="GRCh38_65536"
BUCKET_LEN=65536

BUILD_BUCKETMAP=0

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
mkdir -p "${BENCHMARK_PATH}/../../benchmark/index"
mkdir -p "${BENCHMARK_PATH}/log"
mkdir -p "${BENCHMARK_PATH}/output"

cd "${BENCHMARK_PATH}/../../benchmark/index"

# Build index
#/usr/bin/time -o "${BENCHMARK_PATH}/log/bucketmap_index.time" -v ../../../build/bucketmap -x -i "${INDEX_INDICATOR}_bucketmap" &> "${BENCHMARK_PATH}/log/bucketmap_index.log"


for dist in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
    /usr/bin/time -o "${BENCHMARK_PATH}/log/bucketmap_map.time" -v ../../../build/bucketmap --version-check 0 -r 300 -s 20 -e 0.6 -b 15 -d ${dist} -i "${INDICATOR}_bucketmap" -q ${QUERY_FILE} -o "${BENCHMARK_PATH}/output/bucketmap_${dist}_1_map.sam" &> "${BENCHMARK_PATH}/log/bucketmap_${dist}_1_map.log"
    #/usr/bin/time -o "${BENCHMARK_PATH}/log/bucketmap_map.time" -v ../../../build/bucketmap --version-check 0 -r 300 -s 20 -e 0.6 -b 0 -d ${dist} -i "${INDICATOR}_bucketmap" -q ${QUERY_FILE} -o "${BENCHMARK_PATH}/output/bucketmap_${dist}_0_map.sam" &> "${BENCHMARK_PATH}/log/bucketmap_${dist}_0_map.log"
done

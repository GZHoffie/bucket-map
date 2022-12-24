#!/bin/bash


# the path to the fasta/fastq file is required to be an absolute path
FASTQ_PATH=$1
INDEX_NAME=$2


# go to the index path
cd ./index/
echo $(pwd)

# run bowtie2 to map the reads
/usr/bin/time -o ../log/bowtie2_map.time -v bowtie2 -x ${INDEX_NAME} -U ${FASTQ_PATH} -S ../output/bowtie2_map.sam > ../log/bowtie2_map.log
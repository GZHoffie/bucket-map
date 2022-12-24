#!/bin/bash

FASTA_PATH=$1 # the path to the fasta file is required to be an absolute path
INDEX_NAME=$2

# go to the index directory
cd ./index

# run the indexing method of bowtie2
/usr/bin/time -o ../log/bowtie2_index.time -v bowtie2-build ${FASTA_PATH} ${INDEX_NAME} > ../log/bowtie2_index.log
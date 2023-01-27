#!/bin/bash

# Script to copy the local files to Computing clusters of SoC.

CLUSTER_ADDRESS="192.168.51.113"
USER_NAME="guzh@comp.nus.edu.sg"

BUCKET_MAP_SRC="/home/zhenhao/bucket-map"
GENOME_DATA_PATH="/mnt/d/genome"

# copy bucket map source files to the cluster
#scp -r ${BUCKET_MAP_SRC} "${USER_NAME}@${CLUSTER_ADDRESS}:~/bucket-map"

# copy the data to the cluster
scp -r ${GENOME_DATA_PATH} "${USER_NAME}@${CLUSTER_ADDRESS}:~/genome_data"
#!/bin/bash

## Install mapping tools
# Need zlib, gcc, g++, make, cmake
sudo apt install -y zlib1g-dev gcc g++ make cmake

# Bowtie2
git clone https://github.com/BenLangmead/bowtie2.git
cd bowtie2; make
echo "export PATH=\${PATH}:$(pwd)" >> ~/.bashrc
cd ..

# BWA
git clone https://github.com/lh3/bwa.git
cd bwa; make
echo "export PATH=\${PATH}:$(pwd)" >> ~/.bashrc
cd ..

# Subread
git clone https://github.com/ShiLab-Bioinformatics/subread.git
cd subread/src
make -f Makefile.Linux
cd ..
echo "export PATH=\${PATH}:$(pwd)/bin" >> ~/.bashrc
cd ..

# Minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make
echo "export PATH=\${PATH}:$(pwd)" >> ~/.bashrc
cd ..

# GraphMap
git clone https://github.com/isovic/graphmap.git  
cd graphmap  
make modules  
make  
echo "export PATH=\${PATH}:$(pwd)/bin/Linux-x64" >> ~/.bashrc
cd ..

# enable search in PATH
source ~/.bashrc


# Samtools, for converting bam/sam files
sudo apt install -y samtools

## Import datasets
mkdir -p mapping_data

# Install `datasets`, NCBI command line tool
mkdir tools
cd tools
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
chmod +x ./datasets
echo "export PATH=\${PATH}:$(pwd)" >> ~/.bashrc
cd ..

# enable search in PATH
source ~/.bashrc

# Human genome GRCh38 (3.1 Gb), E Coli genome (5.8 Mb)
cd mapping_data
for NCBI_ID in GCF_000001405.26 GCA_000005845.2
do
    datasets download genome accession ${NCBI_ID}
    unzip ncbi_dataset.zip
    mv ncbi_dataset/data/${NCBI_ID}/*.fna ./
    rm -r README.md ncbi_dataset ncbi_dataset.zip
done
cd ..
# Pseudomonas Genus



# Note: need to use ./delete_invalid_bases.sh to remove `N` in the genome


## Install simulation tools
# DWGSIM
# Need libcurses
sudo apt-get install -y libcurses-ocaml-dev

git clone https://github.com/nh13/DWGSIM.git
cd DWGSIM
git submodule init
git submodule update
make
echo "export PATH=\${PATH}:$(pwd)" >> ~/.bashrc
cd ..

## Generate simulated datasets
cd mapping_data

# Remove the `N` characters
~/bucket-map/bucket_map/benchmark/delete_invalid_bases.sh GCF_000001405.26_GRCh38_genomic.fna GRCh38_adjusted.fna

dwgsim -N 1000000 -1 300 -2 0 GCA_004358405.1_ASM435840v1_genomic.fna EColi_sim
zcat EColi_sim.bwa.read1.fastq.gz >> EColi_sim.bwa.read1.fastq
rm *.fastq.gz

dwgsim -N 1000000 -1 300 -2 0 GRCh38_adjusted.fna GRCh38_sim
zcat GRCh38_sim.bwa.read1.fastq.gz >> GRCh38_sim.bwa.read1.fastq
rm *.fastq.gz

dwgsim -N 10000 -1 300 -2 0 GRCh38_adjusted.fna GRCh38_sim_short
zcat GRCh38_sim_short.bwa.read1.fastq.gz >> GRCh38_sim_short.bwa.read1.fastq
rm *.fastq.gz

cd ..


# pbsim3
git clone https://github.com/yukiteruono/pbsim3.git
cd pbsim3
./configure
make
sudo make install
cd ..


## Download real read datasets
# Install sratools, if haven't done so
sudo apt-get install sra-toolkit

# Download datasets
# Illumina MiSeq, Homo sapiens
prefetch -v SRR25821753	

# NextSeq 2000, Escherichia coli K-12
prefetch -v SRR24524113
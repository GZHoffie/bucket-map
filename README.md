![](https://www.comp.nus.edu.sg/templates/t3_nus2015/images/assets/logos/logo.png)

# **BucketMap**: New Cache-efficient Alignment-free Short Read Mapping Tool

This repository contains my notes and work on the MComp Dissertation at School of Computing, National University of Singapore.

## Installation

This tool is build with [Seqan3](https://docs.seqan.de/seqan/3-master-user/index.html). To properly build the package, you need to have GCC >= 11.3, G++ and CMake installed.

```bash
git clone --recurse-submodules https://github.com/GZHoffie/bucket-map.git
cd ./bucket-map/bucket_map/seqan3
git submodule update --init
```

## Usage
# Benchmarking of DNA Read Mapping Tools

To find time and space usage, we can use the following command.

```shell
/usr/bin/time -v build/bucket_map_benchmark
```

## To-dos

1. Try with bowtie2, BWA.
2. Benchmark
   - running time, 
   - peak memory usage,
   - number of page faults during the run,
   - mapping accuracy,
   - mapping precision (false positive rate).
3. Try to find the dependency on
   - read length, 
   - reference genome length.
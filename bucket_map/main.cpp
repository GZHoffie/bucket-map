#include <seqan3/core/debug_stream.hpp>
#include "bucket_map.h"


int main()
{
    bucket_mapper * mapper = new bucket_mapper("test", 10000, 100);
    mapper->index("/home/zhenhao/mcomp-dissertation/sequence_sample.fasta", "test");



    seqan3::debug_stream << "Hello World!\n";
    return 0;
}
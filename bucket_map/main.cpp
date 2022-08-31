#include "q_gram_map.h"

#include <seqan3/search/kmer_index/shape.hpp>

using seqan3::operator""_shape;

int main()
{
    q_gram_mapper<700> map(10000, 100, 0b11101001010011_shape);
    map.read("/home/zhenhao/mcomp-dissertation/sequence_sample.fasta", "/home/zhenhao/mcomp-dissertation/build/test");
}
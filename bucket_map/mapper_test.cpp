#include "indexer/bucket_hash_indexer.h"
#include "mapper/q_gram_mapper.h"
#include "tools/short_read_simulator.h"
#include "./utils.h"

#include <map>
#include <chrono>

int main() {
    using namespace seqan3::literals;

    std::vector<seqan3::dna4> text{"ACGTAGC"_dna4};
    unsigned int k = 4;

    auto hashes = text | seqan3::views::kmer_hash(seqan3::shape{0b1101_shape});
    auto hashes_rev_comp_seqan = text | std::views::reverse | seqan3::views::complement | seqan3::views::kmer_hash(seqan3::shape{0b1101_shape});
    seqan3::debug_stream << hashes_rev_comp_seqan << "\n";

    auto hashes_rev_comp = hashes | std::views::transform([&](unsigned int hash) {
        return hash_reverse_complement(hash, k);
    });

    seqan3::debug_stream << hashes << "\n" << hashes_rev_comp << "\n";


    auto r = std::ranges::iota_view{0, 10};
    seqan3::debug_stream << r << "\n";
    


    
    

}
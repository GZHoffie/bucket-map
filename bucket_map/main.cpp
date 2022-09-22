#include "mapper/q_gram_map.h"

#include <seqan3/search/kmer_index/shape.hpp>

#include <cstdlib>
#include <ctime>

using seqan3::operator""_shape;

class short_read_simulator {
private:
    std::vector<std::vector<seqan3::dna4>> bucket_sequence;
    int bucket_length, read_length;

    // TODO: Add error in the short read


public:
    short_read_simulator(int bucket_len, int read_len) {
        bucket_length = bucket_len;
        read_length = read_len;
    }

    void read(std::filesystem::path const & fasta_file_name) {
        seqan3::sequence_file_input reference_genome{fasta_file_name};
        for (auto && record : reference_genome) {
            // Divide the record into buckets
            float total_length = (float) record.sequence().size();
            int num_buckets = (int) ceil(total_length / bucket_length);
            
            // read each bucket
            for (int i = 0; i < num_buckets; i++) {
                int start = i * bucket_length;
                int end = start + bucket_length + read_length;
                if (end > record.sequence().size()) {
                    end = record.sequence().size();
                }
                std::vector<seqan3::dna4> seq(&record.sequence()[start], &record.sequence()[end]);
                bucket_sequence.push_back(seq);
            }
        }
    }

    std::pair<std::vector<seqan3::dna4>, int> sample() {
        int bucket = rand() % bucket_sequence.size();
        std::vector<seqan3::dna4> current_bucket = bucket_sequence[bucket];
        int start = rand() % bucket_length;
        int end = start + read_length;
        if (end > current_bucket.size()) {
            end = current_bucket.size();
        }
        
        std::vector<seqan3::dna4> sample_sequence(current_bucket.begin() + start, current_bucket.begin() + end);
        return std::make_pair(sample_sequence, bucket);
    }
};


int main() {
    srand(time(NULL));

    q_gram_mapper<700> map(10000, 100, 0b11101001010011_shape, 10, 2);
    map.read("/home/zhenhao/mcomp-dissertation/sequence_sample.fasta");
    map.store("/home/zhenhao/mcomp-dissertation/build/sequence_sample");

    short_read_simulator sim(10000, 100);
    sim.read("/home/zhenhao/mcomp-dissertation/sequence_sample.fasta");
    int correct = 0;
    int total_size = 0;
    for (int i = 0; i < 100000; i++) {
        auto sample = sim.sample();
        int bucket = std::get<1>(sample);
        std::vector<seqan3::dna4> sequence = std::get<0>(sample);
        std::vector<int> buckets = map.query_sequence(sequence);

        if (std::find(buckets.begin(), buckets.end(), bucket) != buckets.end()) {
            correct++;
        }
        total_size += buckets.size();

    }
    std::cout << correct << std::endl;
    std::cout << ((float) total_size) / 100000 << std::endl;

}
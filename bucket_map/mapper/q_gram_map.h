#include "../index/bucket_index.h"

#include <string>
#include <vector>
#include <bitset>
#include <cmath>

#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

using seqan3::operator""_shape;


template<unsigned int NUM_BUCKETS>
class fault_tolerate_filter {
private:
    unsigned int num_fault_tolerance;
    std::vector<std::bitset<NUM_BUCKETS>> filters;

    std::vector<int> _set_bits(int index) {
        /**
         * @brief Find all set bits in filters[index]
         * @param index indicates which filter we want to check.
         */
        std::vector<int> res;
        if (index >= num_fault_tolerance) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The input index "
                                 << index << " exceeds the fault tolerance level." << '\n';
            return res;
        } else if (filters[index].none()) {
            return res;
        }
        for (int i = filters[index]._Find_first(); i < NUM_BUCKETS; i = filters[index]._Find_next(i)) {
            res.push_back(i);
        }
        return res;
    }

public:
    fault_tolerate_filter(unsigned int fault) {
        num_fault_tolerance = fault;
        for (int i = 0; i < num_fault_tolerance; i++) {
            std::bitset<NUM_BUCKETS> filter;
            filters.push_back(filter.set());
        }
    }

    void reset() {
        for (int i = 0; i < num_fault_tolerance; i++) {
            filters[i].set();
        }
    }

    void read(std::bitset<NUM_BUCKETS> input) {
        /**
         * @brief Read a set of bits and filter out the indices at which they are 0.
         * @param input a bitset with size NUM_BUCKETS indicating whether a certain k-mer appear
         *              in each bucket. (0 at position i indicates that the k-mer doesn't exist 
         *              in bucket i)
         */
        // Record how many time we see a 0 for each bucket
        for(unsigned int i = 0; i < num_fault_tolerance - 1; i++) {
            filters[i] &= (filters[i+1] | input);
        }
        // Modify filters[size-1]
        filters[num_fault_tolerance - 1] &= input;
    }

    std::vector<int> best_results() {
        /**
         * @brief Return the buckets that contains the most number of k-mers.
         */
        std::vector<int> res;
        for (int i = num_fault_tolerance-1; i >= 0; i--) {
            res = _set_bits(i);
            if (!res.empty()) {
                return res;
            }
        }
        return res;
    }



};





template<unsigned int NUM_BUCKETS>
class q_gram_mapper {
private:
    // q-gram index related information
    seqan3::shape q_gram_shape;
    std::vector<std::bitset<NUM_BUCKETS>> q_grams_index;
    unsigned int q;
    unsigned int size;

    // bucket index related information
    unsigned int bucket_length;
    unsigned int read_length;

    // mapper related information
    unsigned int num_samples;
    unsigned int num_fault_tolerance;

    // filter that filter out the most possible bucket
    fault_tolerate_filter<NUM_BUCKETS>* filter;


    void _insert_into_bucket(std::vector<seqan3::dna4> sequence, unsigned int bucket_num) {
        /**
         * @brief Read the sequence, extract all the q-grams and store in `q_grams_index`.
         * @param sequence the corresponding sequence for the bucket.
         * @param bucket_num an integer indicating the position of the bucket. Should be
         *                   between 0 and NUM_BUCKETS - 1.
         */
        auto q_gram = seqan3::views::kmer_hash(q_gram_shape);

        // Extract all q_grams from sequence
        for (auto && value : sequence | q_gram) {
            q_grams_index[value].set(bucket_num);
        }
    }

    std::vector<unsigned char> _bitset_to_bytes(const std::bitset<NUM_BUCKETS>& bs){
        /**
         * @brief Convert bitset into 8-byte chars.
         * Adopted from https://stackoverflow.com/a/7463972.
         */
        std::vector<unsigned char> result((NUM_BUCKETS + 7) >> 3);
        for (unsigned int j = 0; j < NUM_BUCKETS; j++)
            result[j>>3] |= (bs[j] << (j & 7));
        return result;
    }

    std::bitset<NUM_BUCKETS> _bitset_from_bytes(const std::vector<unsigned char>& buf) {
        /**
         * @brief Convert 8-byte chars to bitset.
         * Adopted from https://stackoverflow.com/a/7463972.
         */
        assert(buf.size() == ((NUM_BUCKETS + 7) >> 3));
        std::bitset<NUM_BUCKETS> result;
        for (unsigned int j = 0; j < NUM_BUCKETS; j++)
            result[j] = ((buf[j>>3] >> (j & 7)) & 1);
        return result;
    }



public:
    q_gram_mapper(unsigned int bucket_len, unsigned int read_len, seqan3::shape shape, 
                  unsigned int samples, unsigned int fault) {
        // initialize private variables
        bucket_length = bucket_len;
        read_length = read_len;
        
        q_gram_shape = shape;
        size = std::ranges::size(shape);
        q = shape.count();
        seqan3::debug_stream << "[INFO]\t\t" << "Set q-gram shape to be: " 
                             << shape << " with number of effective characters: " << q << '\n';
        
        num_samples = samples;
        num_fault_tolerance = fault;

        // initialize filter
        filter = new fault_tolerate_filter<NUM_BUCKETS>(num_fault_tolerance);
    }

    ~q_gram_mapper() {
        delete filter;
    }

    void read(std::filesystem::path const & fasta_file_name) {
        /**
         * @brief Read the fasta file, store the q_gram of each bucket in `q_grams_index`.
         *        Store the information in the `index_directory`.
         * @param fasta_file_name the name of the file containing reference genome.
         */
        // q_grams_index should be empty
        if (!q_grams_index.empty()) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The q-gram index is not empty. Terminating read.\n";
            return;
        }
        // initialize q_gram index
        int total_q_grams = (int) pow(4, q);
        for (int i = 0; i < total_q_grams; i++) {
            std::bitset<NUM_BUCKETS> q_gram_bucket;
            q_grams_index.push_back(q_gram_bucket);
        }
        // Read the genome
        seqan3::sequence_file_input reference_genome{fasta_file_name};
        unsigned int bucket_num = 0;
        for (auto && record : reference_genome) {
            // Divide the record into buckets
            float total_length = (float) record.sequence().size();
            int num_buckets = (int) ceil(total_length / bucket_length);
            seqan3::debug_stream << "[INFO]\t\t" << record.id() << " with length " << (int) total_length
                                 << " divided into " << num_buckets << " buckets.\n";
            
            // read each bucket
            for (int i = 0; i < num_buckets; i++) {
                int start = i * bucket_length;
                int end = start + bucket_length + read_length;
                if (end > record.sequence().size()) {
                    end = record.sequence().size();
                }
                std::vector<seqan3::dna4> bucket_sequence(&record.sequence()[start], &record.sequence()[end]);
                _insert_into_bucket(bucket_sequence, bucket_num);
                bucket_num++;
            }
        }
        seqan3::debug_stream << "[INFO]\t\t" << "Total number of buckets: " 
                             << bucket_num << "." << '\n';
    }

    void store(std::filesystem::path const & index_directory) {
        /**
         * @brief Store the created index inside the index_directory.
         * @param index_directory the directory to store the q_gram_count_file.
         */
        // q_grams_index should not be empty
        if (q_grams_index.empty()) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The q-gram index is empty. No file is created.\n";
            return;
        }
        // Create directory if directory is not created yet.
        // Return if the index files already exist.
        if (!std::filesystem::create_directories(index_directory)) {
            seqan3::debug_stream << "[WARNING]\t" << "The specified index directory "
                                 << index_directory << " is already created." << '\n';
            for (const auto& entry : std::filesystem::directory_iterator(index_directory)) {
                if (entry.path().extension() == ".qgram" || entry.path().extension() == ".pattern") {
                    seqan3::debug_stream << "[ERROR]\t\t" << "The q-gram file " << entry.path() << " already exists" 
                                         << " in the specified directory. Terminating store." << '\n';
                    return;
                }
            }
        }
        // Store the q-gram index in the directory
        std::ofstream index_file(index_directory / "index.qgram");
        for (const auto &i : q_grams_index) {
            for (const auto &c : _bitset_to_bytes(i)) {
                index_file << c;
            }
            index_file << "\n";
        }
        seqan3::debug_stream << "[INFO]\t\t" << "The bucket q-gram index is stored in: " 
                             << index_directory / "index.qgram" << ".\n";
        
        // Store the q-gram pattern in the directory
        std::ofstream pattern_file(index_directory / "index.pattern");
        pattern_file << q_gram_shape;
        seqan3::debug_stream << "[INFO]\t\t" << "The q-gram shape is stored in: " 
                             << index_directory / "index.pattern" << "." << '\n';
    }


    void load(std::filesystem::path const & index_directory) {
        /**
         * @brief Look for index and pattern file inside the index_directory,
         *        read the files and store the values in class attribute.
         * @param index_directory the directory that stores the q_gram_count_file.
         */
        // q_grams_index should be empty
        if (!q_grams_index.empty()) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The q-gram index is not empty. Terminating load.\n";
            return;
        }
        // initialize q_gram index
        int total_q_grams = (int) pow(4, q);
        for (int i = 0; i < total_q_grams; i++) {
            std::bitset<NUM_BUCKETS> q_gram_bucket;
            q_grams_index.push_back(q_gram_bucket);
        }
        // Read the index file
        std::ifstream file(index_directory / "index.qgram");
        std::string line;
        while (std::getline(file, line)) {
            std::vector<unsigned char> index(line.begin(), line.end());
            q_grams_index.push_back(_bitset_from_bytes(index));
        }
        seqan3::debug_stream << "[INFO]\t\t" << "Successfully loaded " 
                             << index_directory / "index.qgram" << "." << '\n';
    }


    std::vector<int> query(std::vector<int> q_gram_hash) {
        /**
         * @brief From `q_grams_index`, determine where the sequence may be coming from.
         * @param q_gram_hash the vector containing all hash values of q-grams in the
         *                    query sequence.
         * @returns a vector of integers indicating the possible regions that the sequence
         *          may belong to.
         */
        // q_grams_index should not be empty
        if (q_grams_index.empty()) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The q-gram index is empty. Cannot accept query.\n";
            std::vector<int> res;
            return res;
        }
        // Reset the filter
        filter->reset();

        // insert those samples into the filter
        for (int h : q_gram_hash) {
            filter->read(q_grams_index[h]);
        }
        return filter->best_results();
    }

    std::vector<int> query_sequence(std::vector<seqan3::dna4> sequence) {
        /**
         * @brief From `q_grams_index`, determine where the sequence may be coming from.
         * @param sequence A dna4 vector containing the query sequence.
         * @returns a vector of integers indicating the possible regions that the sequence
         *          may belong to.
         * TODO: implement query with argument being a file containing short reads, together with quality.
         */
        auto q_gram = seqan3::views::kmer_hash(q_gram_shape);
        auto hash_values = sequence | q_gram;
        std::vector<int> samples;
        // Randomly sample `num_samples` q-grams for query.
        std::sample(hash_values.begin(), hash_values.end(), 
                    std::back_inserter(samples), num_samples,
                    std::mt19937{std::random_device{}()});
        return query(samples);
    }


};
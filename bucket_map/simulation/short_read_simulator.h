#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/alphabet/all.hpp>

#include <cstdlib>
#include <ctime>
#include <random>

using seqan3::operator""_dna4;

class short_read_simulator {
    /**
     * @brief A tool to randomly generate short reads based on a reference genome.
     *        It allows setting of substitution/indel error rates and random generation
     *        of errors.
     */
private:
    // variables storing the sequence information
    std::vector<std::vector<seqan3::dna4>> bucket_sequence;
    int bucket_length, read_length;

    // Error generation
    std::vector<seqan3::dna4> neucleotides{ 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4 };
    std::poisson_distribution<int>* substitution_dist; 
    std::poisson_distribution<int>* insertion_dist;
    std::poisson_distribution<int>* deletion_dist;

    // Random number generator
    std::mt19937* gen;

    // Error generating functions
    void _add_substitution(std::vector<seqan3::dna4>& sequence) {
        sequence[rand() % sequence.size()] = neucleotides[rand() % 4];
    }

    void _add_insertion(std::vector<seqan3::dna4>& sequence) {
        sequence.insert(sequence.begin() + rand() % sequence.size(), neucleotides[rand() % 4]);
    }

    void _add_deletion(std::vector<seqan3::dna4>& sequence) {
        sequence.erase(sequence.begin() + rand() % sequence.size());
    }

public:
    short_read_simulator(int bucket_len, int read_len, float substitution_rate = 0, 
                         float insertion_rate = 0, float deletion_rate = 0) {
        // set a random seed
        srand(time(NULL));
        
        // Set indexing-related numbers
        bucket_length = bucket_len;
        read_length = read_len;
        
        // Set the error rates
        substitution_dist = new std::poisson_distribution<int>(substitution_rate * read_length);
        insertion_dist = new std::poisson_distribution<int>(insertion_rate * read_length);
        deletion_dist = new std::poisson_distribution<int>(deletion_rate * read_length);
    }

    ~short_read_simulator() {
        delete substitution_dist, insertion_dist, deletion_dist;
    }

    void add_errors(std::vector<seqan3::dna4>& sequence, int substitutions, int deletions, int insertions) {
        /**
         * @brief Add errors to a short read sequence.
         * @param sequence the short read sequence to be modified.
         * @param substitutions number of substitutions to make.
         * @param deletions number of deletions to make.
         * @param insertions number of insertions to make.
         */
        for (int i = 0; i < deletions; i++) _add_deletion(sequence);
        for (int i = 0; i < insertions; i++) _add_insertion(sequence);
        for (int i = 0; i < substitutions; i++) _add_substitution(sequence);
    }

    void simulate_errors(std::vector<seqan3::dna4>& sequence) {
        /**
         * @brief Do simulation and generate the number of substitution, deletion and insertion randomly
         *        from poisson distributions with lamda = substitution_rate, deletion_rate and 
         *        insertion_rate, respectively.
         * @note the number of errors generated might always be 0 if error rates are too small. Might consider
         *       setting them higher than actual error rate.
         */
        std::random_device rd;
        std::mt19937 gen(rd());
        add_errors(sequence,
                   (*substitution_dist)(gen),
                   (*deletion_dist)(gen),
                   (*insertion_dist)(gen));
    }

    void read(std::filesystem::path const & fasta_file_name) {
        /**
         * @brief Read the fasta file and store the sequences for each bucket. The reference
         *        genome is then used for short read generation.
         * @param fasta_file_name the name of the sequence file to be read.
         */
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

    std::pair<std::vector<seqan3::dna4>, int> sample(bool simulate_error = true) {
        /**
         * @brief Take a sample short read from the reference genome and add errors to it.
         * @param simulate_error whether we add random error to the sequence or not.
         * @return a pair storing <the generated short read, the bucket id it actually belongs to>.
         */
        int bucket = rand() % bucket_sequence.size();
        std::vector<seqan3::dna4> current_bucket = bucket_sequence[bucket];
        int size = current_bucket.size();
        int start = 0;
        if (current_bucket.size() > read_length + 1) {
            start = rand() % (current_bucket.size() - read_length - 1);
        }
        int end = start + read_length;
        if (end > current_bucket.size()) {
            end = current_bucket.size();
        }
        std::vector<seqan3::dna4> sample_sequence(current_bucket.begin() + start, current_bucket.begin() + end);
        if (simulate_error) {
            simulate_errors(sample_sequence);
        }
        return std::make_pair(sample_sequence, bucket);
    }
};
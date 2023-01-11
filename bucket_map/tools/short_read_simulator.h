#ifndef BUCKET_MAP_SHORT_READ_SIMULATOR_H
#define BUCKET_MAP_SHORT_READ_SIMULATOR_H

#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/alphabet/all.hpp>

#include <cstdlib>
#include <ctime>
#include <random>
#include <tuple>

#include "../utils.h"

using seqan3::operator""_dna4;
using seqan3::operator""_cigar_operation;

class short_read_simulator {
    /**
     * @brief A tool to randomly generate short reads based on a reference genome.
     *        It allows setting of substitution/indel error rates and random generation
     *        of errors.
     */
private:
    // variables storing the sequence information
    std::vector<std::vector<seqan3::dna4>> bucket_sequence;
    std::vector<std::pair<int, int>> bucket_ids;
    int bucket_length, read_length;

    // Error generation
    std::vector<seqan3::dna4> neucleotides{ 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4 };
    std::poisson_distribution<int>* substitution_dist; 
    std::poisson_distribution<int>* insertion_dist;
    std::poisson_distribution<int>* deletion_dist;

    // Random number generator
    std::mt19937* gen;

    // Error generating functions
    void _add_substitution(std::vector<seqan3::dna4>& sequence, CIGAR& cigar) {
        unsigned int index = rand() % sequence.size();
        seqan3::dna4 orig_nt = sequence[index];
        seqan3::dna4 new_nt = neucleotides[rand() % 4];
        while (new_nt == orig_nt) {
            new_nt = neucleotides[rand() % 4];
        }
        sequence[index] = new_nt;
        cigar.replace(index, 'X'_cigar_operation);
    }

    void _add_insertion(std::vector<seqan3::dna4>& sequence, CIGAR& cigar) {
        unsigned int index = rand() % sequence.size();
        sequence.insert(sequence.begin() + index, neucleotides[rand() % 4]);
        cigar.insert(index, 'I'_cigar_operation);
    }

    void _add_deletion(std::vector<seqan3::dna4>& sequence, CIGAR& cigar) {
        unsigned int index = rand() % sequence.size();
        sequence.erase(sequence.begin() + index);
        cigar.insert(index, 'D'_cigar_operation);
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

    void add_errors(std::vector<seqan3::dna4>& sequence, CIGAR& cigar, 
                    int substitutions, int deletions, int insertions) {
        /**
         * @brief Add errors to a short read sequence.
         * @param sequence the short read sequence to be modified.
         * @param cigar the ground truth cigar string.
         * @param substitutions number of substitutions to make.
         * @param deletions number of deletions to make.
         * @param insertions number of insertions to make.
         */
        for (int i = 0; i < deletions; i++) _add_deletion(sequence, cigar);
        for (int i = 0; i < insertions; i++) _add_insertion(sequence, cigar);
        for (int i = 0; i < substitutions; i++) _add_substitution(sequence, cigar);
    }

    void simulate_errors(std::vector<seqan3::dna4>& sequence, CIGAR& cigar) {
        /**
         * @brief Do simulation and generate the number of substitution, deletion and insertion randomly
         *        from poisson distributions with lamda = substitution_rate, deletion_rate and 
         *        insertion_rate, respectively.
         * @note the number of errors generated might always be 0 if error rates are too small. Might consider
         *       setting them higher than actual error rate.
         */
        std::random_device rd;
        std::mt19937 gen(rd());
        add_errors(sequence, cigar,
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
        int index = 0;
        int reference_id = -1;
        std::string last_id = "";
        auto operation = [&](std::vector<seqan3::dna4> seq, std::string id) {
            bucket_sequence.push_back(seq);
            if (id != last_id) {
                index = 0;
                last_id = id;
                reference_id++;
            } 
            bucket_ids.push_back(std::make_pair(reference_id, index));
            index++;
        };
        iterate_through_buckets(fasta_file_name, bucket_length, read_length, operation);
    }

    std::tuple<std::vector<seqan3::dna4>, int, int, CIGAR> sample(bool simulate_error = true) {
        /**
         * @brief Take a sample short read from the reference genome and add errors to it.
         * @param simulate_error whether we add random error to the sequence or not.
         * @return a tuple storing 3 values: <the generated short read, the bucket id it actually belongs to,
         *         the starting point of the short read>.
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
        // assign the original cigar sequence
        CIGAR cigar(sample_sequence.size(), '='_cigar_operation);

        // insert errors
        if (simulate_error) {
            simulate_errors(sample_sequence, cigar);
        }
        return std::make_tuple(sample_sequence, bucket, start, cigar);
    }

    void generate_fastq_file(std::filesystem::path output_path, std::string indicator, unsigned int size,  
                             bool simulate_error = true) {
        /**
         * @brief Generate a fastq file containing short reads.
         * @param output_path the directory where we output the fastq file and the answer file.
         * @param indicator a string that indicate the name of output fastq/answer file.
         * @param size the number of short reads to be simulated.
         */
        // Create directory if directory is not created yet.
        // Return if the index files already exist.
        if (!std::filesystem::create_directories(output_path)) {
            seqan3::debug_stream << "[WARNING]\t" << "The specified output directory "
                                 << output_path << " is already created." << '\n';
            for (const auto& entry : std::filesystem::directory_iterator(output_path)) {
                if (entry.path() == indicator + ".fastq" || entry.path() == indicator + ".ground_truth") {
                    seqan3::debug_stream << "[ERROR]\t\t" << "The fastq file or ground truth file " << entry.path() << " already exists" 
                                         << " in the specified directory. Terminating generation." << '\n';
                    return;
                }
            }
        }

        std::ofstream fastq_file(output_path / (indicator + ".fastq"));
        std::ofstream bucket_gt_file(output_path / (indicator + ".bucket_ground_truth"));
        std::ofstream pos_gt_file(output_path / (indicator + ".position_ground_truth"));
        for (unsigned int i = 0; i < size; i++) {
            fastq_file << "@" << i << "\n";
            // generate sequence
            auto res = sample(simulate_error);
            auto & [sequence, bucket, offset, cigar] = res;
            for (auto nt : sequence) {
                fastq_file << nt.to_char();
            }
            // insert a quality string.
            fastq_file << "\n+\n" << std::string(sequence.size(), 'E') << "\n";
            // record the ground truth.
            bucket_gt_file << bucket << " " << offset << " " << cigar.to_string() << "\n";
            // record the true locations
            auto true_position = bucket_ids[bucket];
            pos_gt_file << std::get<0>(true_position) << " " << std::get<1>(true_position) * bucket_length + offset + 1 << " " << cigar.to_string() << "\n";
        }

        seqan3::debug_stream << "[INFO]\t\t" << "The generated fastq file is stored in: " 
                             << output_path / (indicator + ".fastq") << ".\n";
        seqan3::debug_stream << "[INFO]\t\t" << "The ground truth for buckets and offsets is stored in: " 
                             << output_path / (indicator + ".bucket_ground_truth") << ".\n";
        seqan3::debug_stream << "[INFO]\t\t" << "The ground truth for exact locations is stored in: " 
                             << output_path / (indicator + ".position_ground_truth") << ".\n";
    }

};

#endif
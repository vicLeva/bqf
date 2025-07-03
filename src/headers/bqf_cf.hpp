#ifndef BQF_CF
#define BQF_CF

#include "bqf_ec.hpp"
#include <gtest/gtest.h>

typedef enum {
    text, binary, stream
} output_mode_t;

class Bqf_cf : public Bqf_ec {
    uint64_t counter = 0;
    output_mode_t mode;
    std::string str_buffer = "";
    std::vector<uint64_t> bin_buffer;
    uint64_t buffer_iterator = 0;

public:
    Bqf_cf(){};
    /** 
     * \brief Constructor that instantiates a BackpackCQF from quotient and remainder sizes
     * \param q_size The desired size of quotient, will induce filter's size
     * \param k The length of a k_mer
     * \param z The length difference between a k-mer and the inserted s-mer
     * \param mode with which the output shall be written (text : human-readable, binary : binary encoding of dna, stream : not in a file not yet added)
     * \param verb to print on-going operations in stdout (default: false)
     */
    Bqf_cf(uint64_t q_size, uint64_t k, output_mode_t mode = text, bool verb=false);

    Bqf_cf(uint64_t max_memory, output_mode_t mode = text, bool verb = false);


    /**
     * @brief filters a fasta/q file : writes all k-mers present more than once in an output file
     * \param files a vector of all the fastx files from which the kmers should be read
     * \param output is the file in which to write redundant kmers
     */
    void filter_fastx_file(std::vector<std::string> files, std::string output);


    /**
     * @brief if the stream mode is the one chosen, this function can be called to get the next kmer
     * \returns the following solid kmer
     */
    uint64_t yield_kmer(void);

private :
    //unitary tests
    friend class BqfCfTest;
    FRIEND_TEST(BqfCfTest, SimpleInsert);

    /**
     * @brief inserts all kmers from a DNA sequence and writes all k-mers present more than once in a stream
     * \param sequence is the sequence from which the kmers are read
     */
    void insert_from_sequence(std::string sequence);

    /**
     * \brief adds an occurence of an element in a certain position
     * \param position where the element should be counted once more
     */
    bool add_one_to_counter(uint64_t position);
    /**
     * \brief adds an occurence of a number in the BQF
     * \param number that should be inserted
     */
    bool is_second_insert(uint64_t number);
    /**
     * \brief adds an occurence of an encoded in the BQF, and if it has 
     * already been inserted exactly once, writes the kmer in a buffer
     * \param coded_kmer that should be inserted
     */
    void insert_kmer(uint64_t coded_kmer);
    /**
     * @brief enumerates all hashes of solid (appearing at least twice) kmers in the BQF
     * \returns a vector of all hashes of the solid kmers
     */
    std::vector<uint64_t> enumerate_solid();

    
    
};

#endif
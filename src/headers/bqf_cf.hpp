#ifndef BQF_CF
#define BQF_CF

#include "bqf_ec.hpp"

class Bqf_cf : public Bqf_ec {
    //int counter = 0;
public:
    Bqf_cf(){};
    /** 
     * \brief Constructor that instantiates a BackpackCQF from quotient and remainder sizes
     * \param q_size The desired size of quotient, will induce filter's size
     * \param k The length of a k_mer
     * \param z The length difference between a k-mer and the inserted s-mer
     * \param verb to print on-going operations in stdout (default: false)
     */
    Bqf_cf(uint64_t q_size, uint64_t k, bool verb=false);


    /**
     * @brief filters a fasta/q file : writes all k-mers present more than once in an output file
     * \param files a vector of all the fastx files from which the kmers should be read
     * \param output is the file in which to write redundant kmers
     */
    void filter_fastx_file(std::vector<std::string> files, std::string output);

    /**
     * @brief inserts all kmers from a DNA sequence and writes all k-mers present more than once in a stream
     * \param sequence is the sequence from which the kmers are read
     * \param output is the stream in which kmers present more than once can be written
     */
    void insert_from_sequence(std::string sequence, std::string& to_write);

    

    /**
     * \brief adds an occurence of an element in a certain position
     * \param position where the element should be counted once more
     * \returns if the element has already been inserted exactly once before
     */
    bool is_second_add_to_counter(uint64_t position);
    /**
     * \brief adds an occurence of a number in the BQF
     * \param number that should be inserted
     * \returns if the number has already been inserted exactly once before
     */
    bool is_second_insert(uint64_t number);
    /**
     * \brief adds an occurence of a kmer (in string format or already encoded) in the BQF, and if it has 
     * already been inserted exactly once, writes the kmer in output
     * \param kmer that should be inserted
     * \param output where the kmer will be written if it is already present once in the BQF
     */
    void is_second_insert(std::string kmer, std::ofstream& output);
    void is_second_insert(uint64_t coded_kmer, std::string& to_write);


    
};

#endif
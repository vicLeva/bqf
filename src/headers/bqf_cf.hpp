#ifndef BQF_CF
#define BQF_CF

#include "bqf_ec.hpp"

class Bqf_cf : public Bqf_ec {



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
     * @brief inserts a list of smers in the BQF and copies in a file all smers that appear more than twice
     * \param input is the file from which to read the smers
     * \param output is the file in which to write redundant smers
     */
    void insert_from_file_and_filter(std::string input, std::string output);

    void filter_fastx_file(std::vector<std::string> files, std::string output);

    void insert_from_sequence(std::string sequence, std::string output);

    

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
     * \brief adds an occurence of an smer in the BQF, and if it has already been inserted exactly once,
     * writes the smer in output
     * \param smer that should be inserted
     * \param output where the smer will be written if it is already present once in the BQF
     */
    void is_second_insert(std::string smer, std::ofstream& output);


    
};

#endif
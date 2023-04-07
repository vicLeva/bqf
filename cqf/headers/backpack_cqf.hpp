#ifndef BACKPACK_CQF_HPP
#define BACKPACK_CQF_HPP

#include <stdint.h> 
#include <map>

#include "filter.hpp"


/**
 * \brief Represents a variant of the Counting quotient filter
 *
 * This class implements a coutning rank & select quotient filter, it supports insert, remove, query and
 * enumerate operations. It can be instantiated by giving the constructor a memory limit (in MB)
 * or a precise size for quotient and remainder aswell as a bit number for the counter. The counter for each inserted kmer
 * is a few bits at the end of its remainder. These bits can be an exact count or an order of magnitude.
 * Inherited from RSQF class.
 */
class Backpack_cqf : public Rsqf{
    
    public:
    /*  
        ================================================================
        CONSTRUCTORS
        ================================================================
    */ 

    Backpack_cqf();

    /** 
     * \brief Constructor that instantiates a BackpackCQF from quotient and remainder sizes
     * \param q_size The desired size of quotient, will induce filter's size
     * \param r_size The desired size of quotient (usually 64-q_size)
     * \param c_size The desired size of remainders counters
     * \param verbose to print on-going operations in stdout
     */
    Backpack_cqf(uint64_t q_size, uint64_t r_size, uint64_t c_size, bool verbose);

    /** 
     * \brief Constructor that deduces quotient and remainder sizes from the desired struct size
     * Max memory is an upper bound, in reality the filter will 
     * adjust its size and will probably be smaller than the given value
     * \param max_memory The desired (maximum) size of the Rsqf (in MBytes)
     * \param c_size The desired size of remainders counters
     * \param verbose to print on-going operations in stdout
     */
    Backpack_cqf(uint64_t max_memory, uint64_t c_size, bool verbose);


    void insert(uint64_t number, uint64_t count = 1);

    uint64_t query(uint64_t number);

    bool remove(uint64_t number, uint64_t count = 1);

    std::map<uint64_t, uint64_t> enumerate();


    protected:

    uint64_t count_size;

    uint64_t find_quotient_given_memory(uint64_t max_memory, uint64_t count_size);

    void add_to_counter(uint64_t position, uint64_t rem_count);

    void sub_to_counter(uint64_t position, uint64_t count);

    uint64_t get_remainder(uint64_t position, bool w_counter = false);

    uint64_t remainder(uint64_t num) const;
};

#endif
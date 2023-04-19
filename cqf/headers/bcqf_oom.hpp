#ifndef BACKPACK_CQF_OOM_HPP
#define BACKPACK_CQF_OOM_HPP

#include <stdint.h> 
#include <map>

#include "filter.hpp"


/**
 * \brief Represents a variant of the CQF: The Backpack counting filter
 *
 * This class implements a counting rank & select quotient filter, it supports insert, remove, query and
 * enumerate operations. It can be instantiated by giving the constructor a memory limit (in MB)
 * or precise sizes for quotient and remainder aswell as a bit number for the counter. The counter for each inserted kmer
 * is a few bits at the end of its remainder. These bits code for an interval (see Bcqf_ec for exact counting instead).
 * Inherited from RSQF class.
 */
class Bcqf_oom : public Rsqf{
    
    public:
    /*  
        ================================================================
        CONSTRUCTORS
        ================================================================
    */ 

    Bcqf_oom();

    /** 
     * \brief Constructor that instantiates a BackpackCQF from quotient and remainder sizes
     * \param q_size The desired size of quotient, will induce filter's size
     * \param r_size The desired size of quotient (usually 64-q_size)
     * \param c_size The desired size of remainders counters
     * \param verbose to print on-going operations in stdout (default: false)
     */
    Bcqf_oom(uint64_t q_size, uint64_t r_size, uint64_t c_size, bool verbose=false);

    /** 
     * \brief Constructor that deduces quotient and remainder sizes from the desired struct size
     * Max memory is an upper bound, in reality the filter will 
     * adjust its size and will probably be smaller than the given value
     * \param max_memory The desired (maximum) size of the Rsqf (in MBytes)
     * \param c_size The desired size of remainders counters
     * \param verbose to print on-going operations in stdout (default: false)
     */
    Bcqf_oom(uint64_t max_memory, uint64_t c_size, bool verbose=false);

    /** 
     * \brief Insert a NEW number in the filter alongside with his count (reduced to the power of 2 just below). 
     * 
     * This function inserts a NEW number in the BackpackCQF. 
     * If a number has the same quotient but different remainders, it stores remainders in a monothonic way
     * (i.e. each remainder in a run is greater or equal than the predecessor).
     * If the filter is full, new insertions are authomatically discarded.
     * When adding a new distinct element, all remainders and runend bits are shifted right of 1 position.
     * 
     * \param number to insert
     * \param count number of occurences of the element to insert (default: 1)
     */
    void insert(uint64_t number, uint64_t count = 1);

    /** 
     * \brief query a number from the filter.
     * 
     * It queries a number. It first checks if the occupied bit of the quotient is set to 1. If so it scans
     * in a linear way the remainders of the run associated to this element. If it find the remainder it returns its count else 0.
     * Stops immediately if the filter is empty
     * Because the count is a power of 2, the real count is something between 2^count and 2^(count+1) 
     * 
     * \param number Number to query
     * \return the abundance of the given number in the filter (the lower bound )
     */
    uint64_t query(uint64_t number);

    /** 
     * \brief Removes (if present) a number from the filter
     * 
     * This method removes an element from the filter. At the beginning it works like a query. If it finds the 
     * searched element, it find the rightmost remainder it has to shift to mantain the QF organized and then shift
     * everything to remove the element. The rightmost remainder can be the one before the FUS (First Unused Slot) or another remainder 
     * saved in a slot before it. 
     * 
     * \param number element to remove
     * \return 1 if the value has been found in the process, 0 if the element was absent
     **/
    bool remove(uint64_t number);

    /** 
     * \brief Enumerate every element that has been inserted in the filter (possibly hashes)
     * 
     * This method iterates over every slot in the filter, for occupied ones it gets the positions of the corresponding run.
     * Then it computes for every remainder in the run, the original number inserted (by concatenating the remainder value
     * and the quotient value (of the run)) and pushes it into the unordered_set alongside with its abundance
     * 
     * \return a uint_64t map, linking every originally inserted hash to its abundance in the filter
     **/ 
    std::map<uint64_t, uint64_t> enumerate();


    private:

    /** 
     * \brief size in bits of the counter, will determine filter's size in addition to remainder's size
     */
    uint64_t count_size;

    /** 
     * \brief Deduce a quotient size from the memory occupation limit and counter's size
     * 
     * The filter is divided into blocks: each block contains 64 quotients, and occupies 3 words of metadata and 'r' words of remainders 
     * (since ther are 64 remainders inside) --> (r + c + 3) words == ((64 - q + c) + 3) words;
     * 
     * \param max_memory Max size to occupy with the RSQF (in MBytes). 
     * \param count_size size of the counter in bits. "c" in the filter size calculus 
     * 
     * \return Quotient size in bits
     **/
    uint64_t find_quotient_given_memory(uint64_t max_memory, uint64_t count_size);

    /** 
     * \brief returns the remainder slot associated to the requested quotient
     * 
     * The difference with RSQF.get_remainder() is that sometimes we want only the remainder (for finding the element),
     * and sometimes the remainder alongside with its counter (to update the value)
     * 
     * \param position quotient 
     * \param w_counter bool value to know if we retrieve only remainder (false) or remainder + count value (true) (default: false)
     * 
     * \return an uint64 with the value stored in the slot
     */
    uint64_t get_remainder(uint64_t position, bool w_counter = false);

    /** 
     * \brief Computes the lower bound of the original count as a power of 2
     *
     * If the computed value overflows the number of bits allocated for the counter,
     * then (2^counter_size)-1 is inserted as the abundance. When queried, this will return 2^((2^counter_size)-1)
     * 
     * \param count the abundance we want to insert/query
     * \return the highest power of 2 that is still less than count (or (2^count_size)-1 if overflow) 
     */
    uint64_t inf_bound(uint64_t count);
};

#endif
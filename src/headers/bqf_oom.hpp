#ifndef BQF_OOM_HPP
#define BQF_OOM_HPP

#include <stdint.h> 
#include <string> 
#include <map>

#include "abstract_bqf.hpp"


/**
 * \brief Represents a variant of the CQF: The Backpack counting filter
 *
 * This class implements a counting rank & select quotient filter, it supports insert, remove, query and
 * enumerate operations. It can be instantiated by giving the constructor a memory limit (in MB)
 * or precise sizes for quotient and remainder aswell as a bit number for the counter. The counter for each inserted kmer
 * is a few bits at the end of its remainder. These bits code for an interval (see Bcqf_ec for exact counting instead).
 * Inherited from RSQF class.
 */
class Bqf_oom : public Bqf{
    
    public:
    /*  
        ================================================================
        CONSTRUCTORS
        ================================================================
    */ 

    Bqf_oom();

    /** 
     * \brief Constructor that instantiates a BackpackCQF from quotient and remainder sizes
     * \param q_size The desired size of quotient, will induce filter's size
     * \param r_size The desired size of quotient (usually 64-q_size)
     * \param c_size The desired size of remainders counters
     * \param verbose to print on-going operations in stdout (default: false)
     */
    Bqf_oom(uint64_t q_size, uint64_t c_size, uint64_t k, uint64_t z, bool verb=false);

    /** 
     * \brief Constructor that deduces quotient and remainder sizes from the desired struct size
     * Max memory is an upper bound, in reality the filter will 
     * adjust its size and will probably be smaller than the given value
     * \param max_memory The desired (maximum) size of the Rsqf (in MBytes)
     * \param c_size The desired size of remainders counters
     * \param verbose to print on-going operations in stdout (default: false)
     */
    Bqf_oom(uint64_t max_memory, uint64_t c_size, bool verbose=false);


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
     * \brief Removes (if present) a kmer from the filter
     * 
     * \param kmer element to remove
     * \return 1 if the value has been found in the process, 0 if the element was absent
     **/
    bool remove(std::string kmer);


    private:

    /** 
     * \brief Adds the insertion abundance to the filter abundance
     * 
     * In the oom count case, we skip because we cant add to an order of magnitude
     * 
     * \param position the position of the remainder to update
     * \param rem_count the new abundance to add
     **/
    void add_to_counter(uint64_t position, uint64_t rem_count) override;

    /** 
     * \brief Computes the lower bound of the original count as a power of 2
     *
     * If the computed value overflows the number of bits allocated for the counter,
     * then (2^counter_size)-1 is inserted as the abundance. When queried, this will return 2^((2^counter_size)-1)
     * 
     * \param count the abundance we want to insert/query
     * \return the highest power of 2 that is still less than count (or (2^count_size)-1 if overflow) 
     */
    uint64_t process_count(uint64_t count) override;

    static Bqf_oom load_from_disk(const std::string& filename);
};

#endif
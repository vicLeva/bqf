#ifndef BACKPACK_CQF_EC_HPP
#define BACKPACK_CQF_EC_HPP

#include <stdint.h> 
#include <string.h> 
#include <map>

#include "filter.hpp"


/**
 * \brief Represents a variant of the CQF: The Backpack counting filter
 *
 * This class implements a counting rank & select quotient filter, it supports insert, remove, query and
 * enumerate operations. It can be instantiated by giving the constructor a memory limit (in MB)
 * or precise sizes for quotient and remainder aswell as a bit number for the counter. The counter for each inserted kmer
 * is a few bits at the end of its remainder. These bits code for an exact count (see Bcqf_oom for order of magnitude instead).
 * Inherited from RSQF class.
 */
class Bcqf_ec : public Rsqf{
    
    public:
    /*  
        ================================================================
        CONSTRUCTORS
        ================================================================
    */ 

    Bcqf_ec();

    /** 
     * \brief Constructor that instantiates a BackpackCQF from quotient and remainder sizes
     * \param q_size The desired size of quotient, will induce filter's size
     * \param r_size The desired size of quotient (usually 64-q_size)
     * \param c_size The desired size of remainders counters
     * \param verbose to print on-going operations in stdout (default: false)
     */
    Bcqf_ec(uint64_t q_size, uint64_t r_size, uint64_t c_size, bool verbose=false);

    /** 
     * \brief Constructor that deduces quotient and remainder sizes from the desired struct size
     * Max memory is an upper bound, in reality the filter will 
     * adjust its size and will probably be smaller than the given value
     * \param max_memory The desired (maximum) size of the Rsqf (in MBytes)
     * \param c_size The desired size of remainders counters
     * \param verbose to print on-going operations in stdout (default: false)
     */
    Bcqf_ec(uint64_t max_memory, uint64_t c_size, bool verbose=false);

    /** 
     * \brief Insert a number in the filter alongside with his count. 
     * 
     * This function inserts a number in the BackpackCQF. IF the number is already present, it just adds the count to the already present one.
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
     * TODO: allow user to choose hash function
     * \brief Insert a kmer in the filter alongside with his count. 
     * 
     * This function inserts a kmer in the BackpackCQF. IF the number is already present, it just adds the count to the already present one.
     * It is advised that the kmer is in a canonical form, it will be hashed then inserted.
     * 
     * \param kmer to insert
     * \param count kmer abundance
     */
    void insert(std::string kmer, uint64_t count);

    /** 
     * \brief Insert every kmer + abundance of a kmer count software output file (eg KMC)
     * 
     * This function will read every line of the file to insert every pair <kmer, count> into the BCQF
     * 
     * \param file path to file
     */
    void insert(std::string file);

    /** 
     * \brief query a number from the filter.
     * 
     * It queries a number. It first checks if the occupied bit of the quotient is set to 1. If so it scans
     * in a linear way the remainders of the run associated to this element. If it find the remainder it returns its exact count else 0.
     * Stops immediately if the filter is empty
     * 
     * \param number Number to query
     * \return the abundance of the given number in the filter
     */
    uint64_t query(uint64_t number);

    /** 
     * \brief query a kmer from the filter.
     * 
     * Every smer of the kmer will be queried, and the smallest count amongst them will 
     * be returned (see fimpera)
     * 
     * \param kmer the kmer to query
     * \param k the kmer size, k-s+1 smers will be effectively queried
     * \return the abundance of the given kmer in the filter
     */
    uint64_t query(std::string kmer, uint64_t k);

    /** 
     * \brief Removes (if present) a number from the filter
     * 
     * This method removes an element from the filter. At the beginning it works like a query. If it finds the 
     * searched element, it find the rightmost remainder it has to shift to mantain the QF organized and then shift
     * everything to remove the element. The rightmost remainder can be the one before the FUS (First Unused Slot) or another remainder 
     * saved in a slot before it. 
     * A simple substraction is performed if count < filter_abundance, else the element is removed from the Backpack_CQF.
     * 
     * \param number element to remove
     * \param count abundance value to remove (default: 1)
     * \return 1 if the value has been found in the process, 0 if the element was absent
     **/
    bool remove(uint64_t number, uint64_t count = 1);
    
    /** 
     * \brief Removes (if present) a kmer from the filter
     * 
     * This method removes a kmer from the filter. 
     * A simple substraction is performed if count < filter_abundance, else the element is removed from the Backpack_CQF.
     * 
     * \param kmer element to remove
     * \param count abundance value to remove (default: 1)
     * \return 1 if the value has been found in the process, 0 if the element was absent
     **/
    bool remove(std::string kmer, uint64_t count = 1);

    /** 
     * \brief Enumerate every element that has been inserted in the filter (possibly hashes)
     * 
     * This method iterates over every slot in the filter, for occupied ones it gets the positions of the corresponding run.
     * Then it computes for every remainder in the run, the original number inserted (by concatenating the remainder value
     * and the quotient value (of the run)) and pushes it into the unordered_set alongside with its abundance
     * 
     * \return a string to uint_64t map, linking every originally inserted kmer to its abundance in the filter
     **/ 
    std::map<std::string, uint64_t> enumerate();


    private:

    /** 
     * \brief size in bits of the counter, will determine filter's size in addition to remainder's size
     */
    uint64_t count_size;

    /** 
     * \brief size in bits of the hashes that will be inserted
     */
    uint64_t hash_size;
    
    /** 
     * \brief k, supposed to be hash_size/2
     */
    uint64_t kmer_size;


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
     * \brief Adds the insertion abundance to the filter abundance
     * 
     * In the exact count case, if an element is already present, we simply add the count newly inserted
     * to the one already present in the filter. If it overflows the counter size, then the abundance is
     * set to (2^counter_size)-1
     * 
     * \param position the position of the remainder to update
     * \param rem_count the new abundance to add
     **/
    void add_to_counter(uint64_t position, uint64_t rem_count);

    /** 
     * \brief substractes an abundance to originally present abundance
     * 
     * In the exact count case, we can remove only a certain amount of already present abundance
     * without removing the element itself, by substracting the abundance already present by "count"
     * 
     * \param position the position of the remainder to update
     * \param count the value to substract
     **/
    void sub_to_counter(uint64_t position, uint64_t count);

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
};

#endif
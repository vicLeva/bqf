#ifndef BQF_HPP
#define BQF_HPP

#include <algorithm>
#include <limits>
#include <map>

#include "rsqf.hpp"

/**
 * \brief Abstract class of BQF, user can use else Bqf_ec or Bqf_oom
 *
 */
class Bqf : public Rsqf {
private:
    /** 
     * \brief size in bits of the counter, will determine filter's size in addition to remainder's size
     */
    uint64_t count_size;
    
    /** 
     * \brief s, supposed to be hash_size/2
     */
    uint64_t smer_size;

    uint64_t kmer_size;


    /** 
     * \brief Insert every kmer + abundance of a kmer count software output file (eg KMC)
     * 
     * This function will read every line of the file to insert every pair <kmer, count> into the BQF
     * 
     * \param file path to file
     */
    void insert(std::string file);

    /** 
     * TODO: allow user to choose hash function
     * \brief Insert a kmer in the filter alongside with his count. 
     * 
     * This function inserts a kmer in the BackpackCQF.
     * It is advised that the kmer is in a canonical form, it will be hashed then inserted.
     * 
     * \param kmer to insert
     * \param count kmer abundance
     */
    void insert(std::string kmer, uint64_t count);

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

    void query(std::ifstream& infile, std::ofstream& outfile);

    /** 
     * \brief query a sequence from the filter.
     * 
     * Every abundance of each kmer of the query sequence will be queried. This is done following Fimpera scheme
     * and the smallest count amongst them will be returned. 
     * 
     * \param seq the sequence to query
     * \param k the kmer size, k-s+1 smers will be effectively queried
     * \return the abundance of the given kmer in the filter
     */
    result_query query(std::string seq);

    /** 
     * \brief query a number from the filter.
     * 
     * It queries a number. It first checks if the occupied bit of the quotient is set to 1. If so it scans
     * in a linear way the remainders of the run associated to this element. If it 
     * finds the remainder it returns its exact count (or order of magnitude, depending of the filter used ) else 0.
     * Stops immediately if the filter is empty
     * 
     * \param number Number to query
     * \return the abundance of the given number in the filter
     */
    uint64_t query(uint64_t number);
    

    /** 
     * \brief Enumerate every element that has been inserted in the filter (possibly hashes)
     * 
     * This method iterates over every slot in the filter, for occupied ones it gets the positions of the corresponding run.
     * Then it computes for every remainder in the run, the original number inserted (by concatenating the remainder value
     * and the quotient value (of the run)) and pushes it into the unordered_set alongside with its abundance
     * 
     * \return a string to uint_64t map, linking every originally inserted kmer to its abundance in the filter
     **/ 
    std::map<uint64_t, uint64_t> enumerate();

    //tmp public
    void resize(int n); 


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
    virtual void add_to_counter(uint64_t position, uint64_t rem_count) = 0;

    virtual uint64_t insert_process_count(uint64_t count) = 0;
    virtual uint64_t query_process_count(uint64_t count) = 0;

    /** 
     * \brief Serializes BQF on disk without compression, can be loaded later for queries.
     * 
     **/
    void save_on_disk(const std::string& filename);
};

#endif
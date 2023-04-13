#ifndef FILTER_HPP
#define FILTER_HPP

#include <vector>
#include <stdint.h> 
#include <unordered_set>

#include "additional_methods.hpp"


class Backpack_cqf;


/**
 * \brief Represents a RSQF
 *
 * This class implements a rank & select quotient filter, it supports insert, remove, query and
 * enumerate operations. It can be instantiated by giving the constructor a memory limit (in MB)
 * or a precise size for quotient and remainder.
 */
class Rsqf {
    
    public:
    /*  
        ================================================================
        CONSTRUCTORS
        ================================================================
    */ 

    Rsqf();

    /** 
     * \brief Constructor that instantiates a QF from quotient and remainder sizes
     * \param q_size The desired size of quotient, will induce filter's size
     * \param r_size The desired size of remainders (usually 64-q_size)
     * \param verbose to print on-going operations in stdout
     */
    Rsqf(uint64_t q_size, uint64_t r_size, bool verbose = false);

    /** 
     * \brief Constructor that deduces quotient and remainder sizes from the desired struct size
     * Max memory is an upper bound, in reality the filter will 
     * adjust its size and will probably be smaller than the given value
     * \param max_memory The desired (maximum) size of the Rsqf (in MBytes)
     * \param verbose to print on-going operations in stdout
     */
    Rsqf(uint64_t max_memory, bool verbose = false);
    

    /*  
        ================================================================
        HIGH LEVEL PUBLIC OPERATIONS 
        ================================================================
    */ 

    /** 
     * \brief Insert a new number in the filter. 
     * 
     * This function inserts a number in the RSQF. IF the number is already present, it just adds another copy of the number.
     * If a number has the same quotient but different remainders, it stores remainders in a monothonic way
     * (i.e. each reminder in a run is greater or equal than the predecessor).
     * If the filter is full, new insertions are authomatically discarded.
     * When adding a new element, all remainders and runend bits are shifted right of 1 position.
     * 
     * \param number to insert
     */
    void insert(uint64_t number);
    
    /** 
     * \brief Query a number from the filter.
     * 
     * it queries a number. It first checks if the occupied bit of the quotient is set to 1. If so it scans
     * in a linear way the remainders of the run associated to this element. If it find the remainder it returns 1 else 0.
     * Stops immediately if the filter is empty
     * 
     * \param number Number to query
     * \return the presence (1) or the absence (0) of the given number in the filter
     */
    bool query(uint64_t number);

    /** 
     * \brief Remove (if present) a number from the filter
     * 
     * This method removes an element from the filter. At the beginning it works like a query. If it finds the 
     * searched element, it find the rightmost remainder it has to shift to mantain the QF organized and then shift
     * everything to remove the element. The rightmost remainder can be the one before the FUS (First Unused Slot) or another remainder 
     * saved in a slot before it. There are some edge cases where the remainder slots are all occupied up to a certain
     * slot but I cannot shift left part of the remainders beacuse they would be moved in a position that is smaller
     * than the one of their quotient. This vould violate one of the rules of the QF.
     * 
     * \param number value to remove
     * \return 1 if the value has been found (and deleted), 0 if the element was absent
     */
    bool remove(uint64_t number);


    /** 
     * \brief Enumerate every number that has been inserted in the filter (possibly hashes)
     * 
     * This method iterates over every slot in the filter, for occupied ones it gets the positions of the corresponding run.
     * Then it computes for every remainder in the run, the original number inserted (by concatenating the remainder value
     * and the quotient value (of the run)) and pushes it into the unordered_set
     * 
     * \return a uint_64t unordered set of number, present in the filter
     */
    std::unordered_set<uint64_t> enumerate(); 


    /*  
        ================================================================
        DISPLAYING 
        ================================================================
    */ 
    
    /** 
     * \brief Generate a string to represent a block. 
     * The 3 first lines are offset, occ and runend metadata. Then the other lines are remainders.
     * \param block_id The block to print.
     * \param bitformat a flag to print remainders as bitvector instead of numbers.
     * \return The block represented as string
     */
    std::string block2string(size_t block_id, bool bit_format = false);



    protected:
    // ATTRIBUTES

    /** 
     * \brief the filter itself, a list of 64bits words
     * 
     * The filter is a vector of 64bits unsigned integers, the offset, occupied and runends vectors are 64bits
     * integers and remainders are written into 1 64bits integer or over multiple ones.
     */
    std::vector<uint64_t> filter; 
    
    /** 
     * \brief max number of bits occupied by the qf to check no memory leaks
     */
    uint64_t m_num_bits; 

    /** 
     * \brief size in bits of quotients, will determine filter's size
     */
    uint64_t quotient_size; 

    /** 
     * \brief size in bits of remainders
     */
    uint64_t remainder_size; 

    /** 
     * \brief Block size (in machine words (64 bits multiples))
     */
    uint64_t block_size;

    /** 
     * \brief number of blocks the qf is divided into
     */
    uint64_t number_blocks;

    /** 
     * \brief number of elements (non distinct) inserted in the filter
     */
    uint64_t elements_inside;

    /** 
     * \brief used mainly for debug, displays in stdout what's going on in the QF
     */
    bool verbose; 

    /** 
     * \brief unittests
     */
    friend class RsqfTest;
    
    // CLASS METHODS

    /*  
        ================================================================
        GETTERS 
        ================================================================
    */  


    /** 
     * \brief returns the size in bits of the quotient
     */    
    uint64_t get_quot_size();

    /** 
     * \brief returns the number of numbers inserted 
     */ 
    uint64_t get_num_el_inserted();

    /** 
     * \brief check if the occupied bit of the given quotient is 1 or 0.
     * It extracts the related bit of the occupied word in the block of the quotient
     * \param position quotient to query
     * \return true if 1, false if 0.
     */
    bool is_occupied(uint64_t position);

    /** 
     * \brief returns the remainder slot associated to the requested quotient
     * \param position quotient 
     * \return an uint64 with the value stored in the slot
     */
    uint64_t get_remainder(uint64_t position);

    /** 
     * \brief gets the word where the slot starts
     * returns the id of the word in the filter that contains (or starts containing) the
     * remainder slot associated to a certain quotient (by this I mean quotient 155 is associated to
     * the 155th slot (since there is a 0th)).
     * \param quotient the *th quotient
     * \return uint64 returning the id of *th word in the vector.
     */
    uint64_t get_remainder_word_position(uint64_t quotient);

    /** 
     * \brief It returns the bit of word where the slot starts.
     * complementary to "get_remainder_word_position", it gives back the bit in the word where the requested slot
     * starts.
     * \param quotient the *th slot requested
     * \return uint64 with the shift in the word where the slot starts
     */
    uint64_t get_remainder_shift_position(uint64_t quotient);

    /** 
     * \brief start and end positions of quotient's run
     * returns the start and the end of the run that contains all the remainders of the numbers inserted that have 
     * the given quotient. I.e. I have quotient 1 and want to see where the run that contains the remainders that are linked to it in the qf
     * start and ends.
     * I have the start and end of the run of the remainders of the numbers that have that quotient.
     * \param quotient quotient of the numbers where we want to have the associated
     * \return a std::pair with start and end position of the run
     */
    std::pair<uint64_t,uint64_t> get_run_boundaries(uint64_t quotient); //const

    /** 
     * \brief gets the end of the previous run of the selected quotient (or its own if it exists). 
     * It is the equivalent of doing  select(runend_vector, rank(occupieds_vector, quotient))
     * \param quotient quotient to give to the function
     * \return - the position of the end of the previous run
     *         - true if the pos is out of quotient block, else false (useful for first_unused_slot)
     */
    std::pair<uint64_t, bool> get_runend(uint64_t quotient);

    /** 
     * \brief It returns the start of the run associated with quotient in parameter
     * If quotient is not occupied it returns the first slot where  would be inserted a remainder linked to this quotient
     * Counts how many runs begin before the quotient using rank(), then find the position of the end of the last
     * run. The next slot is the beginning of our interest run, if it is before our quotient, then the slot "quotient"
     * is the beginning of the run
     * \param quotient quotient of which we want the beginning of the run 
     * \param occ_bit the information of occupation of the quotient (1=the quotient is occupied) 
     * \return a uint64_t, the position of the beginning of the run
     */
    uint64_t get_runstart(uint64_t quotient, bool occ_bit=1);

    /** 
     * \brief Used in get runstart()
     * for the special case where the quotient if the first slot of a block
     * Because we start counting in the block how many runs start before our quotient (the occupieds
     * we can't do it if it is in first position, so we start by counting the runends.
     * \param quotient quotient of which we want the beginning of the run 
     * \param paj the position after the jump done using the offset of the block 
     * \param offset the information of the offset of the block (might be different of paj)
     * \param occ_bit the information of occupation of the quotient (1=the quotient is occupied) 
     * \return a uint64_t, the position of the beginning of the ru
     */
    uint64_t get_runstart_shift0(uint64_t quotient, uint64_t paj, uint64_t offset, bool occ_bit);


    /** 
     * \brief gives the next remainder 
     * If the one given is at the end, it gives back the one at the beginning of the filter. 
     * It also skips the metadata words. It is used forremainder slot shifting
     * \param current_word current quotient
     * \return uint64 with position of the one after
     */
    uint64_t get_next_remainder_word(uint64_t current_word) const;


    /** 
     * \brief It gives back the quotient after the one given. 
     * If the one given is at the end, it gives back the one at the begining of the filter (0). 
     * \param current_quot current quotient
     * \return uint64 with position of the one after
     */    
    uint64_t get_next_quot(uint64_t current_quot) const;

    /** 
     * \brief it gives back the quotient before the one given. 
     * If the one given is 0, it gives back the one at the end of the filter. 
     * \param current_quot current quotient
     * \return uint64 with position of the one before
     */
    uint64_t get_prev_quot(uint64_t current_quot) const;

    /** 
     * \brief it gives back the block_id before the one given
     * A block is made by the metadata and remainder slots words
     * If the one given is 0, it gives back the id of the one at the end of the filter. 
     * \param current_quot current quotient
     * \return uint64 with position of the one before
     */  
    uint64_t get_prev_block_id(uint64_t current_block) const;

    /** 
     * \brief it gives back the block_id after the one given
     * A block is made by the metadata and remainder slots words
     * If the one given is 0, it gives back the id of the one at the end of the filter. 
     * \param current_quot current quotient
     * \return uint64 with position of the one before
     */  
    uint64_t get_next_block_id(uint64_t current_block) const;

    /** 
     * \brief It returns the word containing the slice of the runends bitvector corresponding to the requested block. 
     * \param current_block chosen block 
     * \return uint64 with the runend word of the particular block.
     */
    uint64_t get_runend_word(uint64_t current_block) const;

    /** 
     * \brief It returns the word containing the slice of the occupieds bitvector corresponding to the requested block. 
     * \param current_block chosen block 
     * \return uint64 with the occupieds word of the particular block.
     */
    uint64_t get_occupied_word(uint64_t current_block) const;

    /** 
     * \brief It returns the word containing the slice of the occupieds bitvector corresponding to the requested block.  
     * \param current_block chosen block 
     * \return uint64 with the offset word of the particular block.
     */
    uint64_t get_offset_word(uint64_t current_block) const;



    /*  
        ================================================================
        SETTERS 
        ================================================================
    */  

    /** 
     * \brief sets the remainder slot associated to the requested quotient to a given value
     * If I want to set the 155th remainder slot to 1, I will use position=155 and value=1.
     * \param position quotient associated to the slot to set 
     * \param value value to store in the slot
     */
    void set_remainder(uint64_t position, uint64_t value);

    /** 
     * \brief It sets the word containing the slice of the runends bitvector corresponding to the requested block to a chosen value
     * \param current_block chosen block
     * \param value new runend word value
     */
    void set_runend_word(uint64_t current_block, uint64_t value);

    /** 
     * \brief It sets the word containing the offset corresponding o the requested block to a chosen value
     * \param current_block chosen block;
     * \param value new offset value;
     * \return uint64 with the runend word of the particular block.
     */
    void set_offset_word(uint64_t current_block, uint64_t value);

    /** 
     * \brief Decrease the chosen block offset value by 1
     * \param current_block chosen block;
     */
    void decrement_offset(uint64_t current_block);

    /** 
     * \brief It sets the bit of the occupieds bitvector corresponding to the chosen block and the bit position in the word.
     * e.g. set 5th bit of 3rd block to 1 or 0
     * \param current_block chosen block 
     * \param value has to be 0 or 1 as it sets a bit
     * \param bit_pos the position of the bit in the occupied word 
     * \return uint64 with the runend word of the particular block.
     */
    void set_occupied_bit(uint64_t current_block, uint64_t value, uint64_t bit_pos);

    /** 
     * \brief It sets the bit of the runend bitvector corresponding to the chosen block and (bit position in the word).
     * e.g. set 5th bit of 3rd block to 1 or 0
     * \param current_block chosen block 
     * \param value has to be 0 or 1 as it sets a bit
     * \param bit_pos the position of the bit in the occupied word 
     * \return uint64 with the runend word of the particular block.
     */
    void set_runend_bit(uint64_t current_block, uint64_t value ,uint64_t bit_pos);

    


    /*  
        ================================================================
        FILTER OPERATIONS 
        ================================================================
    */  

    /** 
     * \brief Deduce a quotient size from the memory occupation limit
     * The filter is divided into blocks: each block contains 64 quotients, and occupies 3 words of metadata and 'r' words of remainders 
     * (since ther are 64 remainders inside) --> (r + 3) words == ((64 - q) + 3) words;
     * There are 2^(q) quotients in the filter. Each block contains 64 quotients --> #blocks = 2^(q)/64. 
     * Since I cannot analitycally get q from the equation MAX_MEM = (2^(q) / 64)((64 - q) + 3),I try values in a for loop. 
     * Starting with q = 63, I lower the size of q till I get the mem of the filter <= of max memory.
     * P.S.: To solve some problems with the comparison of numbers, for large values of q, I divide the 2^(q) by the 
     * number of bytes in a MB. For the smaller values I multiply the MAX_MEM value by the same amount.
     * \param max_memory Max size to occupy with the RSQF (in MBytes). The memory of the RSQF is given by 1 parameter: the quotient size(q), the remainder size(r) is (64 - q).
     * \return Quotient size in bits
     **/
    uint64_t find_quotient_given_memory(uint64_t max_memory);

    /** 
     * \brief Computes the remainder from a number to insert/find in the filter
     * \param num the number (hash) we want to insert/query
     * \return the remainder part that will be inserted or queried in the filter
     */
    uint64_t remainder(uint64_t num) const;

    /** 
     * \brief Computes the quotient from a number to insert/find in the filter
     * \param num the number (hash) we want to insert/query
     * \return the quotient (slot index) where the remainder of the number will be inserted (or found if it was inserted)
     */
    uint64_t quotient(uint64_t num) const;
    

    /** 
     * \brief shift the remainders on 1 slot from start to end and insert a new remainder in the start position
     * It shifts (left) bits in words to leave space for a new insertion. Then inserts the value.
     * \param start_quotient quotient where to start the shifting
     * \param end_quotient quotient where to end the shifting
     * \param next_remainder the value to store in the remainder slot at 'start_quotient' position
     */
    void shift_left_and_set_circ(uint64_t start_quotient,uint64_t end_quotient, uint64_t next_remainder);

    /** 
     * \brief shift the remainders of 1 slot from end to start and clears the 'end' slot. 
     * It is used to remove remainders.
     * It shifts bits right in the words to overwrite the 'start_quotient' slot that has to be removed.
     * \param start_quotient position that will be overwritten when shifting
     * \param end_quotient position that will be zeroed when shifting
     */
    void shift_right_and_rem_circ(uint64_t start_quotient,uint64_t end_quotient);


    /** 
     * \brief used in shift_bits_right_metadata() during a deletion
     */
    void shift_runend_right(uint64_t start_shift, uint64_t end_shift, uint64_t block);

    /** 
     * \brief Metadata function that shifts the bits right in the runend word when there is a deletion.
     * It also handles cases when occupieds and offsets have to be adjourned.
     * \param quotient quotient of the number to delete.
     * \param flag_bit flags (with 0) if the occupied bit has to be zeroed, else it is 1.
     * \param start_position starting position of the shifting.
     * \param end_position end position of the shifting.
     */
    void shift_bits_right_metadata(uint64_t quotient, uint64_t start_position, uint64_t end_position);

    /** 
     * \brief Metadata function that shifts the bits left in the runend word when there is an insertion.
     * It also handles cases when occupieds and offsets have to be adjourned.
     * \param quotient quotient of the number to add.
     * \param flag_bit flags (with 1) if the occupied bit has to be set to 1, else it is 0.
     * \param start_position starting position of the shifting.
     * \param end_position end position of the shifting.
     */
    void shift_bits_left_metadata(uint64_t quotient, uint64_t flag_bit, uint64_t start_position, uint64_t end_position);


    /** 
     * \brief Returns the first unused slot as described in the paper. 
     * Used on insertion and to find the shifting window of deletions.
     * It extracts the related bit of the occupied word in the block of the quotient
     * \param curr_quotient the quotient from where to find the first unused slot
     * \return uint64 of the position of the first unused slot.
     */
    uint64_t first_unused_slot(uint64_t curr_quotient); //const

    /** 
     * \brief Returns the first slot that can't be shifted during a deletion. 
     * It iterates over all the slots after curr_quotient, and stops when it finds else an unused slot
     * or an occupied slot whose run starts at this precise slot. (the run of a quotient can't start before the quotient itself)
     * \param curr_quotient the quotient from where to find the first unshiftable slot
     * \return uint64 of the position of the last shiftable slot.
     */
    uint64_t first_unshiftable_slot(uint64_t curr_quotient); //const

    


    


};  


#endif
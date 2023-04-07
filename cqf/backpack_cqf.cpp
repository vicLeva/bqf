#include <iostream>
#include <cmath>
#include <cassert>

#include "backpack_cqf.hpp" 

// STATIC VARIABLES 
#define MEM_UNIT 64ULL
#define MET_UNIT 3ULL
#define OFF_POS 0ULL
#define OCC_POS 1ULL
#define RUN_POS 2ULL
#define SCALE_INPUT 8388608ULL

using namespace std;

Backpack_cqf::Backpack_cqf(){}

Backpack_cqf::Backpack_cqf(uint64_t q_size, uint64_t r_size, uint64_t c_size, bool verbose){
    assert (c_size < q_size);

    elements_inside = 0;
    quotient_size = q_size;
    remainder_size = r_size + c_size;
    count_size = c_size;
    cout << "quotient_size " << quotient_size << " remainder_size " << remainder_size << " count_size " << count_size << endl; 

    uint64_t num_quots = 1ULL << quotient_size; 
    uint64_t num_of_words = num_quots * (MET_UNIT + remainder_size) / MEM_UNIT; 

    // In machine words
    this->block_size = (3 + this->remainder_size);
    number_blocks = std::ceil(num_quots / MEM_UNIT);

    filter = std::vector<uint64_t>(num_of_words);
    m_num_bits = num_of_words*MEM_UNIT;
}


Backpack_cqf::Backpack_cqf(uint64_t max_memory, uint64_t c_size, bool verbose){
    elements_inside = 0;
    
    // Size of the quotient/remainder to fit into max_memory MB
    quotient_size = find_quotient_given_memory(max_memory, c_size);
    assert (c_size < quotient_size);
    remainder_size = MEM_UNIT - quotient_size + c_size;
    count_size = c_size;

    // Number of quotients must be >= MEM_UNIT
    uint64_t num_quots = 1ULL << quotient_size; //524.288
    uint64_t num_of_words = num_quots * (MET_UNIT + remainder_size) / MEM_UNIT; //393.216

    // In machine words
    this->block_size = (3 + this->remainder_size);
    number_blocks = std::ceil(num_quots / MEM_UNIT);
    
    filter = std::vector<uint64_t>(num_of_words);
    m_num_bits = num_of_words*MEM_UNIT;
}




void Backpack_cqf::insert(uint64_t number, uint64_t count){
    if (elements_inside == number_blocks*MEM_UNIT - 1) return; //100%-1 is max number

    //get quotient q and remainder r
    uint64_t quot = quotient(number);
    uint64_t rem = remainder(number);
    uint64_t rem_count = (rem << count_size) | count;

    if (verbose){
        std::cout << "[INSERT] quot " << quot << std::endl;
        std::cout << "[INSERT] rem " << rem << std::endl;
    }
    

    // GET FIRST UNUSED SLOT
    uint64_t fu_slot = first_unused_slot(quot);
    //assert(fu_slot.remainder == 0)
    assert(get_remainder(fu_slot) == 0);
    
    if (verbose) {
        std::cout << "[INSERT] FUS " << fu_slot << std::endl;
    }

    if (!is_occupied(quot)){
        uint64_t starting_position = get_runstart(quot, 0);

        if (verbose) {
            cout << "starting_position " << starting_position << endl; 
        }
        
        set_occupied_bit(get_block_id(quot), 1, get_shift_in_block(quot));
        shift_bits_left_metadata(quot, 1, starting_position, fu_slot);
        elements_inside++;
        return shift_left_and_set_circ(starting_position, fu_slot, rem_count);
    }
    // IF THE QUOTIENT HAS BEEN USED BEFORE
    // GET POSITION WHERE TO INSERT TO (BASED ON VALUE) IN THE RUN (INCREASING ORDER)
    else{
        if (verbose){
            cout << "occupied" << endl;
        }
        //getting boundaries of the run
        std::pair<uint64_t,uint64_t> boundary = get_run_boundaries(quot);

        if (verbose){
            cout << "boundaries " << boundary.first << " || " << boundary.second << endl;
        }
        
        //find the place where the remainder should be inserted / all similar to a query
        //getting position where to start shifting right
        uint64_t starting_position = boundary.first;
        uint64_t remainder_in_filter = get_remainder(starting_position); 
        
        if (starting_position == boundary.second){ //1 element run
            if (remainder_in_filter < rem) {
                starting_position = get_next_quot(starting_position);
            }
            else if (remainder_in_filter == rem){ 
                    add_to_counter(starting_position, rem_count);
                    return;
            }
        }
        else{
            while(starting_position != boundary.second){
                remainder_in_filter = get_remainder(starting_position); 
                if (remainder_in_filter > rem) {
                    break;
                }
                else if (remainder_in_filter == rem){ 
                    add_to_counter(starting_position, rem_count);
                    return;
                }
                starting_position = get_next_quot(starting_position);
            }


            //last iter, before or after last element
            remainder_in_filter = get_remainder(starting_position); 
            if (remainder_in_filter < rem) {
                starting_position = get_next_quot(starting_position);
            }
            else if (remainder_in_filter == rem){ 
                add_to_counter(starting_position, rem_count);
                return;
            }
        }        
        
        shift_bits_left_metadata(quot, 0, boundary.first, fu_slot);
        // SHIFT EVERYTHING RIGHT AND INSERTING THE NEW REMAINDER
        elements_inside++;
        return shift_left_and_set_circ(starting_position, fu_slot, rem_count);
    }
}


uint64_t Backpack_cqf::query(uint64_t number){
    if (elements_inside == 0) return 0;
    uint64_t quot = quotient(number);
    uint64_t rem = remainder(number);

    if (!is_occupied(quot)) return 0;

    std::pair<uint64_t,uint64_t> boundary = get_run_boundaries(quot);

    // TODO:
    // OPTIMIZE TO LOG LOG SEARCH ?

    uint64_t position = boundary.first;

    while(position != boundary.second){
        uint64_t remainder_in_filter = get_remainder(position);

        if (remainder_in_filter == rem) return get_remainder(position, true) & mask_right(count_size);
        else if (remainder_in_filter > rem) return 0;
        position = get_next_quot(position);
    }

    uint64_t remainder_in_filter = get_remainder(boundary.second); 
    if (remainder_in_filter == rem) return get_remainder(position, true) & mask_right(count_size);

    return 0;
}


bool Backpack_cqf::remove(uint64_t number, uint64_t count){
    if (elements_inside == 0) return 0;
    //get quotient q and remainder r
    uint64_t quot = quotient(number);
    uint64_t rem = remainder(number);

    if (verbose){
        std::cout << "[REMOVE] quot " << quot << std::endl;
        std::cout << "[REMOVE] rem " << rem << std::endl;
    }

    if (is_occupied(quot) == false) return 0;

    std::pair<uint64_t,uint64_t> boundary = get_run_boundaries(quot);

    // TODO:
    // OPTIMIZE TO LOG LOG SEARCH ?
    uint64_t pos_element;
    uint64_t position = boundary.first;
    bool found = false;
    uint64_t remainder_in_filter;

    // GET POSITION
    while(position != boundary.second){
        remainder_in_filter = get_remainder(position); 
        if (remainder_in_filter == rem) {
            pos_element = position;
            found = true;
            break;
        }
        else if (remainder_in_filter > rem) return 0;
        position = get_next_quot(position);
    }
    //last iter
    remainder_in_filter = get_remainder(boundary.second); 
    if (remainder_in_filter == rem) {
        pos_element = position;
        found = true;
        }
    if (!found) return 0; //not found

    //remove some of the ones present
    if (count < (get_remainder(pos_element, true) & mask_right(count_size))){
        sub_to_counter(pos_element, count);
        return 1;
    }

    //remove everything (delete remainder)
    // GET FIRST UNSHIFTABLE SLOT (first quot = runend(quot) or unused)
    uint64_t end_slot = first_unshiftable_slot(quot);

    if (verbose){
        cout << "[REMOVE] FUS " << end_slot << endl;
    }

    if (pos_element == end_slot) {}

    //METADATA
    if(boundary.first == boundary.second) {
        // LAST ELEMENT, ZERO OCC
        set_occupied_bit(get_block_id(quot), 0, get_shift_in_block(quot));

        if (boundary.second == end_slot){
            //MANUALLY HANDLING METADATA (ISOLATED RUN)
            if (get_shift_in_block(quot) == 0){
                decrement_offset(get_block_id(quot));
            }
            set_runend_bit(get_block_id(end_slot), 0, get_shift_in_block(end_slot));
        }
        else{
            //THROUGH SHIFTING EVERYTHING (LAST ELEMENT DELETION STILL SHIFTS NEXT RUNS)
            shift_bits_right_metadata(quot, pos_element, end_slot);
        }
    }
    else {
        //REMOVE ONE ELEMENT IN A RUN, CLASSIC
        shift_bits_right_metadata(quot, pos_element, end_slot);
    }


    // REMAINDERS
    shift_right_and_rem_circ(pos_element, end_slot);

    elements_inside--;
    return 1;
}


std::map<uint64_t, uint64_t> Backpack_cqf::enumerate(){
    std::map<uint64_t, uint64_t> finalSet;
    uint64_t curr_occ;
    
    std::pair<uint64_t, uint64_t> bounds;
    uint64_t cursor;

    uint64_t quotient;
    uint64_t number;

   
    for(uint block = 0; block < number_blocks; ++block){
        curr_occ = get_occupied_word(block);
        if (curr_occ == 0) continue;

        for (uint64_t i=0; i<MEM_UNIT; i++){
            if (curr_occ & 1ULL){ //occupied
                quotient = block*MEM_UNIT + i;
                bounds = get_run_boundaries(quotient);
                cursor = bounds.first;
                while (cursor != (bounds.second)){ //every remainder of the run
                    number = rebuild_number(quotient, get_remainder(cursor), quotient_size);
                    finalSet[number] = get_remainder(cursor, true) & mask_right(count_size);
                    cursor = get_next_quot(cursor);
                }

                number = rebuild_number(quotient, get_remainder(cursor), quotient_size);
                finalSet[number] = get_remainder(cursor, true) & mask_right(count_size);
            }

            curr_occ >>= 1ULL; //next bit of occupied vector
        }
    }

    return finalSet;
}



uint64_t Backpack_cqf::find_quotient_given_memory(uint64_t max_memory, uint64_t count_size){
    uint64_t quotient_size;
    uint64_t curr_m;
    
    for (int i = MEM_UNIT - 1; i > 0; --i){
        quotient_size = i;

        if (quotient_size >= (MEM_UNIT - MET_UNIT)){
            curr_m = ((1ULL << quotient_size)/SCALE_INPUT) * (MEM_UNIT + MET_UNIT - quotient_size + count_size);
            if (max_memory >= curr_m) return quotient_size;
        }
        else{
            curr_m = ((1ULL << quotient_size)) * (MEM_UNIT + MET_UNIT - quotient_size + count_size);
            if (max_memory*SCALE_INPUT >= curr_m) return quotient_size;
        }

    }
    return 0;

}

void Backpack_cqf::add_to_counter(uint64_t position, uint64_t remainder_w_count){
    uint64_t old_rem = get_remainder(position, true);

    uint64_t sum = (old_rem & mask_right(count_size)) + (remainder_w_count & mask_right(count_size));
    if (! (sum < 1ULL << (count_size+1))){
        sum = (1ULL << (count_size+1)) - 1;
    }  
    
    sum |= old_rem & mask_left(MEM_UNIT - count_size);

    set_bits(filter, 
            get_remainder_word_position(position) * MEM_UNIT + get_remainder_shift_position(position), 
            sum, 
            remainder_size);
}

void Backpack_cqf::sub_to_counter(uint64_t position, uint64_t count){
    uint64_t old_rem = get_remainder(position, true);

    cout << "sub " << (old_rem & mask_right(count_size)) << endl;

    uint64_t sub = (old_rem & mask_right(count_size)) - count;

    
    
    sub |= old_rem & mask_left(MEM_UNIT - count_size);

    set_bits(filter, 
            get_remainder_word_position(position) * MEM_UNIT + get_remainder_shift_position(position), 
            sub, 
            remainder_size);
}


uint64_t Backpack_cqf::get_remainder(uint64_t position, bool w_counter ){ //default=false
    uint64_t block = get_block_id(position);
    uint64_t pos_in_block = get_shift_in_block(position);
    uint64_t pos = (block*((MET_UNIT+remainder_size)*MEM_UNIT)+MET_UNIT*MEM_UNIT+pos_in_block*remainder_size); 

    if (w_counter) return get_bits(filter, pos, remainder_size);
    else return get_bits(filter, pos, remainder_size) >> count_size;
}

uint64_t Backpack_cqf::remainder(uint64_t num) const{
    return num >> (MEM_UNIT - remainder_size + count_size);
}


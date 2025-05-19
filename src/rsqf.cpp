#include "rsqf.hpp" 

//FILTER CLASS

using namespace std;

Rsqf::Rsqf(){
}

Rsqf::Rsqf(uint64_t q_size, uint64_t r_size, bool verbose) : verbose(verbose) {
    assert(q_size >= 7);

    quotient_size = q_size;
    remainder_size = r_size;

    uint64_t num_quots = 1ULL << quotient_size; 
    uint64_t num_of_words = num_quots * (MET_UNIT + remainder_size) / MEM_UNIT;
    
    size_limit = num_quots * 0.95;

    // In machine words
    number_blocks = ceil(num_quots / BLOCK_SIZE);

    filter = vector<uint64_t>(num_of_words);
}


Rsqf::Rsqf(uint64_t max_memory, bool verbose) : verbose(verbose) { 
    
    // Size of the quotient/remainder to fit into max_memory MB
    quotient_size = find_quotient_given_memory(max_memory);
    assert(quotient_size >= 7);
    remainder_size = MEM_UNIT - quotient_size;

    // Number of quotients must be >= MEM_UNIT
    uint64_t num_quots = 1ULL << quotient_size; //524.288
    uint64_t num_of_words = num_quots * (MET_UNIT + remainder_size) / MEM_UNIT; 

    // In machine words
    number_blocks = ceil(num_quots / BLOCK_SIZE);
    
    filter = vector<uint64_t>(num_of_words);
}

string Rsqf::block2string(size_t block_id, bool bit_format) {
    stringstream stream;

    // Init the position in filter to the first machine word of the block
    uint64_t position = block_id * (MET_UNIT + this->remainder_size);

    stream << "BLOCK " << block_id;
    stream << "   first quotient: " << block_id*64 << endl;

    // --- Offset ---
    stream << "offset : " << this->filter[position++] << endl;
    
    // --- Occurences ---
    stream << "occ    : ";
    uint64_t occ_register = this->filter[position++];
    for (size_t bit=0 ; bit<MEM_UNIT ; bit++) {
        stream << ((occ_register & 0b1) ? '1' : '0');
        occ_register >>= 1;
        if (bit % 4 == 3)
            stream << "  ";
    } stream << endl;

    // --- Runends ---
    stream << "runends: ";
    uint64_t runs_register = this->filter[position++];
    for (size_t bit=0 ; bit<MEM_UNIT ; bit++) {
        stream << ((runs_register & 0b1) ? '1' : '0');
        runs_register >>= 1;
        if (bit % 4 == 3)
            stream << "  ";
    } stream << endl;


    stream << "BLOCK " << block_id << endl;
    
    stream << "BLOCK " << block_id << endl;
    
    // --- Remainders ---
    // Read all bits of the block one by one
    uint64_t bit_position = 0;
    bool end_remainder = false;
    for (size_t remainder_id=0 ; remainder_id<MEM_UNIT ; remainder_id++) {
        uint64_t current_remainder = 0;

        // Get the remainder bit by bit
        for (size_t remainder_bit=0 ; remainder_bit<this->remainder_size ; remainder_bit++) {
            // Add one bit to the remainder
            size_t byte_position = bit_position / MEM_UNIT;
            uint64_t bit_value = (this->filter[position + byte_position] >> (bit_position % MEM_UNIT)) & 0b1;
            current_remainder += bit_value << remainder_bit;

            // Pretty print the bit format
            if (bit_format) {
                // Beginning of the line
                if (bit_position % MEM_UNIT == 0) {
                    stream << "         ";
                }
                // Remainder separation symbol
                if (end_remainder and remainder_bit == 0)
                    stream << '|';
                // Current bit value
                stream << (bit_value ? '1' : '0');
                // End of the line
                if (bit_position % MEM_UNIT == MEM_UNIT - 1) {
                    stream << endl;
                    end_remainder = false;
                }
                // Separation between bits
                else if (bit_position % 4 == 3) {
                    stream << (end_remainder ? " " : "  ");
                    end_remainder = false;
                }
            }

            // Update the loop
            bit_position += 1;
        }
        // Take into account the end of the remainder for bit pretty printing
        end_remainder = true;

        // Format the string using decimal numbers
        if (! bit_format) {
            if (remainder_id % 4 == 0)
                stream << "         ";
            stream << current_remainder << '\t';
            if (remainder_id % 4 == 3)
                stream << endl;
        }
    }

    return stream.str();
}


/*
HIGH LEVEL OPERATIONS
*/
uint64_t Rsqf::get_quot_size(){
    return quotient_size;
}

uint64_t Rsqf::get_num_el_inserted(){
    return elements_inside;
}

uint64_t Rsqf::find_quotient_given_memory(uint64_t max_memory){ //TO CHANGE, MISS HASH INFORMATION
    uint64_t quotient_size;
    uint64_t curr_m;
    
    for (int i = MEM_UNIT - 1; i > 0; --i){
        quotient_size = i;

        if (quotient_size >= (MEM_UNIT - MET_UNIT)){
            curr_m = ((1ULL << quotient_size)/SCALE_INPUT) * (MEM_UNIT + MET_UNIT - quotient_size);
            if (max_memory >= curr_m) return quotient_size;
        }
        else{
            curr_m = ((1ULL << quotient_size)) * (MEM_UNIT + MET_UNIT - quotient_size);
            if (max_memory*SCALE_INPUT >= curr_m) return quotient_size;
        }

    }
    return 0;

}

using namespace std;

void Rsqf::insert(uint64_t number){
    if (elements_inside+1 == size_limit){
        if (verbose){
            cout << "RESIZING" << endl;
        }
        this->resize(1);    
    }

    //get quotient q and remainder r
    uint64_t quot = quotient(number);
    uint64_t rem = remainder(number);

    if (verbose){
        cout << "[INSERT] quot " << quot << endl;
        cout << "[INSERT] rem " << rem << endl;
    }
    

    // GET FIRST UNUSED SLOT
    uint64_t fu_slot = first_unused_slot(quot);
    assert(get_remainder(fu_slot) == 0);
    
    if (verbose) {
        cout << "[INSERT] FUS " << fu_slot << endl;
    }

    if (!is_occupied(quot)){
        uint64_t starting_position = get_runstart(quot, 0);

        if (verbose) {
            cout << "starting_position " << starting_position << endl; 
        }
        
        set_occupied_bit(get_block_id(quot), 1, get_shift_in_block(quot));
        shift_bits_left_metadata(quot, 1, starting_position,fu_slot);
        elements_inside++;
        return shift_left_and_set_circ(starting_position, fu_slot, rem);
    }
    // IF THE QUOTIENT HAS BEEN USED BEFORE
    // GET POSITION WHERE TO INSERT TO (BASED ON VALUE) IN THE RUN (INCREASING ORDER)
    else{
        if (verbose){
            cout << "occupied" << endl;
        }
        //getting boundaries of the run
        pair<uint64_t,uint64_t> boundary = get_run_boundaries(quot);

        if (verbose){
            cout << "boundaries " << boundary.first << " || " << boundary.second << endl;
        }
        
        //find the place where the remainder should be inserted / all similar to a query
        //getting position where to start shifting right
        uint64_t starting_position = boundary.first;
        uint64_t remainder_in_filter = get_remainder(starting_position); 
        
        if (starting_position == boundary.second){
            if (remainder_in_filter < rem) {
                starting_position = get_next_quot(starting_position);
            }
        }
        else{
            while(starting_position != boundary.second){
                remainder_in_filter = get_remainder(starting_position); 
                if (remainder_in_filter > rem) {
                    break;
                }
                starting_position = get_next_quot(starting_position);
            }
            //last iter, before or after last present element
            remainder_in_filter = get_remainder(starting_position); 
            if (remainder_in_filter <= rem) {
                starting_position = get_next_quot(starting_position);
            }
        }
        
        uint64_t metadata_starting_position = boundary.first; 
        
        
        shift_bits_left_metadata(quot, 0, metadata_starting_position, fu_slot);
        
        // SHIFT EVERYTHING RIGHT AND INSERTING THE NEW REMINDER
        elements_inside++;
        shift_left_and_set_circ(starting_position, fu_slot, rem);
    }
}


bool Rsqf::query(uint64_t number){
    if (elements_inside == 0) return 0;
    //get quotient q and remainder r
    uint64_t quot = quotient(number);
    uint64_t rem = remainder(number);

    if (!is_occupied(quot)) return 0;

    pair<uint64_t,uint64_t> boundary = get_run_boundaries(quot);

    // TODO:
    // OPTIMIZE TO LOG LOG SEARCH ?

    uint64_t position = boundary.first;

    while(position != boundary.second){
        uint64_t remainder_in_filter = get_remainder(position);

        if (remainder_in_filter == rem) return 1;
        else if (remainder_in_filter > rem) return 0;
        position = get_next_quot(position);
    }

    uint64_t remainder_in_filter = get_remainder(boundary.second); 
    if (remainder_in_filter == rem) return 1;

    return 0;
}


bool Rsqf::remove(uint64_t number){

    if (elements_inside == 0) return 0;
    //get quotient q and remainder r
    uint64_t quot = quotient(number);
    uint64_t rem = remainder(number);

    if (verbose){
        cout << "[REMOVE] quot " << quot << endl;
        cout << "[REMOVE] rem " << rem << endl;
    }

    if (is_occupied(quot) == false) return 0;

    pair<uint64_t,uint64_t> boundary = get_run_boundaries(quot);

    // TODO:
    // OPTIMIZE TO LOG LOG SEARCH ?
    uint64_t pos_element;
    uint64_t position = boundary.first;
    bool found = false;
    uint64_t remainder_in_filter;

    // GET POSITION
    while(position != boundary.second){
        remainder_in_filter = get_remainder(position); 
        cout << "rem_infilt " << remainder_in_filter << endl;
        if (remainder_in_filter == rem) {
            pos_element = position;
            found = true;
            break;
        }
        else if (remainder_in_filter > rem) return 0;
        position = get_next_quot(position);
    }
    
    remainder_in_filter = get_remainder(boundary.second); 
    cout << "rem_infilt " << remainder_in_filter << endl;
    if (remainder_in_filter == rem) {
        pos_element = position;
        found = true;
        }
    if (!found) return 0; //not found

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


unordered_set<uint64_t> Rsqf::enumerate(){
    unordered_set<uint64_t> finalSet;
    uint64_t curr_occ;
    
    pair<uint64_t, uint64_t> bounds;
    uint64_t cursor;

    uint64_t quotient;
    uint64_t number;

   
    for(uint block = 0; block < number_blocks; ++block){
        curr_occ = get_occupied_word(block);
        if (curr_occ == 0) continue;

        for (uint64_t i=0; i<BLOCK_SIZE; i++){
            if (curr_occ & 1ULL){ //occupied
                quotient = block*BLOCK_SIZE + i;
                bounds = get_run_boundaries(quotient);
                cursor = bounds.first;
                while (cursor != (bounds.second)){ //every remainder of the run
                    number = rebuild_number(quotient, get_remainder(cursor), quotient_size);
                    finalSet.insert(number);
                    cursor = get_next_quot(cursor);
                }

                number = rebuild_number(quotient, get_remainder(cursor), quotient_size);
                finalSet.insert(number); 
            }

            curr_occ >>= 1ULL; //next bit of occupied vector
        }
    }

    return finalSet;
}


void Rsqf::resize(int n){
    unordered_set<uint64_t> inserted_elements = this->enumerate();

    this->quotient_size += n;
    this->remainder_size -= n;
    
    uint64_t num_quots = 1ULL << this->quotient_size; 
    uint64_t num_of_words = num_quots * (MET_UNIT + this->remainder_size) / MEM_UNIT; 

    this->size_limit = num_quots * 0.95;

    // In machine words
    this->number_blocks = ceil(num_quots / BLOCK_SIZE);

    this->filter = vector<uint64_t>(num_of_words);

    for (const uint64_t& hash: inserted_elements) {
        this->insert(hash);
    }
}



//
uint64_t Rsqf::quotient(uint64_t num) const{
    return num & mask_right(quotient_size);
}

uint64_t Rsqf::remainder(uint64_t num) const{
    return num >> (quotient_size);
}

// REMAINDER OPERATIONS

uint64_t Rsqf::get_remainder(uint64_t position){
    uint64_t block = get_block_id(position);
    uint64_t pos_in_block = get_shift_in_block(position);

    uint64_t pos = block * ((MET_UNIT+remainder_size)*BLOCK_SIZE) + MET_UNIT*BLOCK_SIZE + pos_in_block*remainder_size; 
    // ----------| all remainders slots before current block -----| all metadata slots -| all remainder slots before current slot 

    return get_bits(filter, pos, remainder_size);
}


void Rsqf::set_remainder(uint64_t position, uint64_t value){
    assert(position < number_blocks*BLOCK_SIZE);

    uint64_t block = get_block_id(position);
    uint64_t pos_in_block = get_shift_in_block(position);

    uint64_t pos = block * ((MET_UNIT+remainder_size)*BLOCK_SIZE) + MET_UNIT*BLOCK_SIZE + pos_in_block*remainder_size;

    set_bits(filter, pos, value, remainder_size);
}

uint64_t Rsqf::get_remainder_word_position(uint64_t quotient){
    return (get_block_id(quotient) * (MET_UNIT + remainder_size) + MET_UNIT + ((remainder_size * get_shift_in_block(quotient)) / BLOCK_SIZE)); 
}

uint64_t Rsqf::get_remainder_shift_position(uint64_t quotient){
    return (get_shift_in_block(quotient) * remainder_size ) % BLOCK_SIZE;
}

void Rsqf::shift_left_and_set_circ(uint64_t start_quotient,uint64_t end_quotient, uint64_t next_remainder){
    if (verbose){
        cout << "shift_left_and_set_circ start_quotient " << start_quotient << " end_quotient " << end_quotient << " next_remainder " << next_remainder << endl;
        cout << "next_rem: "; 
        print_bits(next_remainder);
    } 
    assert(start_quotient < ( 1ULL << quotient_size));
    assert(end_quotient < ( 1ULL << quotient_size)); 

    uint64_t mask;
    uint64_t keep;
    uint64_t to_shift;

    // to verify if we are inserting a remainder on 2 words
    uint64_t left_word_shift;
    uint64_t complementary_shift;

    uint64_t curr_word_pos = get_remainder_word_position(start_quotient);
    uint64_t curr_word_shift = get_remainder_shift_position(start_quotient);

    const uint64_t end_word_pos = get_remainder_word_position(get_next_quot(end_quotient));
    const uint64_t end_word_shift = get_remainder_shift_position(get_next_quot(end_quotient));

    while (curr_word_pos != end_word_pos){
        // left bits we will lose after the shift
        to_shift = (filter[curr_word_pos] & mask_left(remainder_size)) >> (64 - remainder_size);

        // mask of what we keep
        mask = mask_right(curr_word_shift);
        keep = filter[curr_word_pos] & mask;
        // shifting operation and keeping what we dont want to be shifted
        filter[curr_word_pos] &= ~mask;
        filter[curr_word_pos] <<= remainder_size;
        filter[curr_word_pos] |= keep;
        filter[curr_word_pos] |= ((next_remainder & mask_right(remainder_size)) << curr_word_shift);

        left_word_shift = 64 - curr_word_shift;

        // if the next_remainder is on 2 words
        // this can only happen if the words are consecutives
        if(left_word_shift < remainder_size){
            complementary_shift = remainder_size - left_word_shift;

            // mask to keep the rest of the 2nd word
            mask = mask_right(complementary_shift);
            keep = filter[curr_word_pos + 1] & mask;
            
            // inserting the missing part of the next_remainder 
            filter[curr_word_pos + 1] &= ~mask;
            filter[curr_word_pos + 1] |= ((next_remainder & mask_right(remainder_size)) >> left_word_shift);

            // next_remainder is the previous remainder that was on the 2 words
            to_shift = (to_shift >> complementary_shift) | (keep << left_word_shift);

            curr_word_shift = complementary_shift;
        } else {
            curr_word_shift = 0;
        }

        // updating the next remainder
        next_remainder = to_shift;

        // we advance by 1 position
        curr_word_pos += 1;
        if (curr_word_pos % (MET_UNIT + remainder_size) == 0)
            curr_word_pos += 3;
        curr_word_pos %= filter.size();
    }

    // mask of what we keep
    mask = mask_left(64ULL - end_word_shift) | mask_right(curr_word_shift);
    keep = filter[curr_word_pos] & mask;
    // shifting operation and keeping what we dont want to be shifted
    filter[curr_word_pos] &= ~mask;
    filter[curr_word_pos] <<= remainder_size;
    filter[curr_word_pos] |= keep;
    filter[curr_word_pos] |= ((next_remainder & mask_right(remainder_size)) << curr_word_shift);
}

void Rsqf::shift_right_and_rem_circ(uint64_t start_quotient,uint64_t end_quotient){
    assert(start_quotient < ( 1ULL << quotient_size));
    assert(end_quotient < ( 1ULL << quotient_size)); 

    uint64_t curr_word_pos = get_remainder_word_position(start_quotient);
    uint64_t curr_word_shift = get_remainder_shift_position(start_quotient);

    uint64_t to_copy_pos = curr_word_pos;
    uint64_t to_copy_shift = curr_word_shift;

    uint64_t to_shift;


    while (start_quotient != end_quotient){
        start_quotient = get_next_quot(start_quotient);

        if (curr_word_shift + remainder_size >= BLOCK_SIZE){
            to_copy_shift = remainder_size - (BLOCK_SIZE - curr_word_shift); //(curr_word_shift + rem_size) % BLOCK_SIZE 
            to_copy_pos = get_next_remainder_word(curr_word_pos);
        }
        else{
            to_copy_shift += remainder_size;
        }

        to_shift = get_bits(filter, to_copy_pos*BLOCK_SIZE + to_copy_shift, remainder_size);
        set_bits(filter, curr_word_pos*BLOCK_SIZE + curr_word_shift, to_shift, remainder_size);

        curr_word_pos = to_copy_pos;
        curr_word_shift = to_copy_shift;
    }

    set_bits(filter, curr_word_pos*BLOCK_SIZE + curr_word_shift, 0, remainder_size);
}

// CIRCULAR FILTER OPERATIONS

uint64_t Rsqf::get_next_remainder_word(uint64_t current_word) const{
    uint64_t current_block = current_word / (MET_UNIT + remainder_size);
    uint64_t pos_in_block = current_word % (MET_UNIT + remainder_size);

    if (pos_in_block != (MET_UNIT + remainder_size - 1)) return ++current_word;
    else{
        uint64_t next_block = get_next_block_id(current_block);
        return next_block * (MET_UNIT + remainder_size) + (MET_UNIT); 
    }
}


uint64_t Rsqf::get_next_quot(uint64_t current_quot) const{
    if (current_quot < number_blocks*BLOCK_SIZE - 1) return ++current_quot;
    else return 0;
}

uint64_t Rsqf::get_prev_quot(uint64_t current_quot) const{
    if (current_quot > 0 ) return --current_quot;
    else return number_blocks*BLOCK_SIZE - 1;
}


uint64_t Rsqf::get_prev_block_id(uint64_t current_block) const{
    if (current_block > 0 ) return --current_block;
    else return number_blocks - 1;
}

uint64_t Rsqf::get_next_block_id(uint64_t current_block) const{
    if (current_block < number_blocks - 1 ) return ++current_block;
    else return 0;
}

uint64_t Rsqf::get_runend_word(uint64_t current_block) const{
    uint64_t runend_id = (current_block *(MET_UNIT+remainder_size)) + RUN_POS;
    return filter[runend_id];
}

uint64_t Rsqf::get_occupied_word(uint64_t current_block) const{
    uint64_t occupied_id = (current_block *(MET_UNIT+remainder_size)) + OCC_POS;
    return filter[occupied_id];
}

uint64_t Rsqf::get_offset_word(uint64_t current_block) const{
    uint64_t offset_id = (current_block *(MET_UNIT+remainder_size)) + OFF_POS;
    assert(filter[offset_id] >= 0);
    return filter[offset_id];
}

void Rsqf::set_runend_word(uint64_t current_block, uint64_t value ){
    uint64_t runend_id = (current_block *(MET_UNIT+remainder_size)) + RUN_POS;
    filter[runend_id] = value;
}

void Rsqf::set_offset_word(uint64_t current_block, uint64_t value ){
    uint64_t offset_id = (current_block *(MET_UNIT+remainder_size)) + OFF_POS;
    filter[offset_id] = value;
}

void Rsqf::decrement_offset(uint64_t current_block){
    uint64_t old_offset = get_offset_word(current_block);
    set_offset_word(current_block, (0 == old_offset)? 0 : old_offset-1);
}

void Rsqf::set_occupied_bit(uint64_t current_block, uint64_t value ,uint64_t bit_pos){
    uint64_t occupied_id = (current_block *(MET_UNIT+remainder_size)) + OCC_POS;
    uint64_t occ_word = get_occupied_word(current_block);

    value &= mask_right(1ULL);
    value <<= bit_pos;
    uint64_t out_value = ((occ_word & mask_right(bit_pos)) | value);

    out_value |= (occ_word & mask_left(MEM_UNIT-bit_pos-1));
    filter[occupied_id] = out_value;

}

void Rsqf::set_runend_bit(uint64_t current_block, uint64_t value ,uint64_t bit_pos){
    uint64_t runend_id = (current_block *(MET_UNIT+remainder_size)) + RUN_POS;
    uint64_t rend_word = get_runend_word(current_block);

    value &= mask_right(1ULL);
    value <<= bit_pos;
    uint64_t out_value = ((rend_word & mask_right(bit_pos)) | value);

    out_value |= (rend_word & mask_left(MEM_UNIT-bit_pos-1));
    filter[runend_id] = out_value;

}

//BITVECTOR AND METADATA OPERATIONS

bool Rsqf::is_occupied(uint64_t position){
    uint64_t block = get_block_id(position);
    uint64_t pos_in_block = get_shift_in_block(position);
    return get_bit_from_word(get_occupied_word(block) ,pos_in_block);
}


uint64_t Rsqf::first_unshiftable_slot(uint64_t curr_quotient){ //const
    pair<uint64_t, bool> rend_pos = get_runend(curr_quotient);
    
    if (verbose){
        cout << "[FUnshftbleS] first runend " << rend_pos.first << " | " <<  rend_pos.second << " (quot " << curr_quotient << ")" << endl;
    }

    if (get_shift_in_block(rend_pos.first) == 0 && !rend_pos.second && get_offset_word(get_block_id(rend_pos.first)) == 0){
        //case we jumped at quot (shift) 0 and runend(shift0) = 0 but 3 cases possible:
        //run overflows by 1 || run of 1 starts at this quotient || no run at the beginning of this block (offset==0)
        //in 3rd case, we return the curr_quot  
        return curr_quotient;
    }

    // rend_pos.second == true means we left the quotient block while looking for runend (handles toricity)
    while( rend_pos.second || curr_quotient < rend_pos.first ){ 
        curr_quotient = get_next_quot(rend_pos.first);

        if (get_runstart(curr_quotient, is_occupied(curr_quotient)) == curr_quotient) return rend_pos.first; //previous one

        rend_pos = get_runend(curr_quotient);

        if (get_shift_in_block(rend_pos.first) == 0 && !rend_pos.second && get_offset_word(get_block_id(rend_pos.first)) == 0){
            break;
        }

        if (verbose){
            cout << "[FUnshftbleS] loop runend " << rend_pos.first << " | " <<  rend_pos.second << " (quot " << curr_quotient << ")" << endl;
        }
    }

    return curr_quotient;
}

uint64_t Rsqf::first_unused_slot(uint64_t curr_quotient){ //const
    pair<uint64_t, bool> rend_pos = get_runend(curr_quotient);
    
    if (verbose){
        cout << "[FUS] first runend " << rend_pos.first << " | " <<  rend_pos.second << " (quot " << curr_quotient << ")" << endl;
    }

    if (get_shift_in_block(rend_pos.first) == 0 && !rend_pos.second && get_offset_word(get_block_id(rend_pos.first)) == 0){
        //case we jumped at quot (shift) 0 and runend(shift0) = 0 but 3 cases possible:
        //run overflows by 1 || run of 1 starts at this quotient || no run at the beginning of this block (offset==0)
        //in 3rd case, we return the curr_quot  
        return curr_quotient;
    }

    // rend_pos.second == true means we left the quotient block while looking for runend (handles toricity)
    while( rend_pos.second || curr_quotient <= rend_pos.first ){ 
        curr_quotient = get_next_quot(rend_pos.first);
        rend_pos = get_runend(curr_quotient);

        if (get_shift_in_block(rend_pos.first) == 0 && !rend_pos.second && get_offset_word(get_block_id(rend_pos.first)) == 0){
            break;
        }

        if (verbose){
            cout << "[FUS] loop runend " << rend_pos.first << " | " <<  rend_pos.second << " (quot " << curr_quotient << ")" << endl;
        }
    }

    return curr_quotient;
}


pair<uint64_t, bool> Rsqf::get_runend(uint64_t quotient){
    uint64_t current_block = get_block_id(quotient);
    uint64_t current_shift = get_shift_in_block(quotient);
    uint64_t offset = get_offset_word(current_block);

    if (current_shift == 0) {
        if (offset <= 1){
            return make_pair(quotient, false);
        }
        else {
            return make_pair((quotient + offset - 1) % (1ULL << quotient_size), 
                                  (offset-1 >= BLOCK_SIZE));
        }
    }

    //paper : d = rank(occupied[i+1, ..., j], j-i-1)
    uint64_t nb_runs = bitrankasm(get_occupied_word(current_block) & mask_left(MEM_UNIT-1),
                                  current_shift);

    

    if (nb_runs == 0) { 
        uint64_t offset_orig = (offset==0) ? 0 : offset-1;
        return make_pair(quotient - current_shift + offset_orig, (offset_orig >= BLOCK_SIZE));
    }

    
    uint64_t pos_after_jump = ((current_block * BLOCK_SIZE) + offset) % (1ULL << quotient_size);
    bool left_block = (current_block != get_block_id(pos_after_jump));

    current_block = get_block_id(pos_after_jump);
    current_shift = get_shift_in_block(pos_after_jump);

    uint64_t runend_mask = mask_left(MEM_UNIT-current_shift); 
    //mask runend[pos_after_jump : end]

    //paper : t = select(runends[i+O(i)+1, ..., end], d)
    uint64_t select_val = bitselectasm(get_runend_word(current_block) & runend_mask, 
                                       nb_runs);
                                       
    
    nb_runs -= bitrankasm(get_runend_word(current_block) & runend_mask, 
                          MEM_UNIT-1);


    while (select_val == BLOCK_SIZE){   
        left_block = true;
        current_block = get_next_block_id(current_block);
        select_val = bitselectasm(get_runend_word(current_block), nb_runs);
        nb_runs -= bitrankasm(get_runend_word(current_block), MEM_UNIT-1);
    }

    return make_pair(current_block*BLOCK_SIZE + select_val, left_block);
}


uint64_t Rsqf::get_runstart(uint64_t quotient, bool occ_bit){
    uint64_t current_block = get_block_id(quotient);
    uint64_t current_shift = get_shift_in_block(quotient);
    uint64_t offset = get_offset_word(current_block);
    offset = (offset==0) ? 0 : offset-1;

    // === JUMP W/ OFFSET ===
    uint64_t select_val;
    uint64_t pos_after_jump = ((current_block * BLOCK_SIZE) + offset) % (1ULL << quotient_size);
    // (pos of block[0] + offset) % nbQuotsMax


    if (current_shift == 0){
        //particular case of shift == 0 
        //(offset jumps to the end of this run)
        return get_runstart_shift0(quotient, pos_after_jump, offset, occ_bit);
    } 
    
    else {
        // === COMPUTE d ===
        uint64_t nb_runs; 
        //nb_runs = d in the paper, number of runs that have to end before ours can start
        nb_runs = bitrankasm(get_occupied_word(current_block),
                             current_shift-1);   

        if (verbose){
            //cout << "[RUNSTART] d " << nb_runs << endl; 
        }
            

        current_shift = get_shift_in_block(pos_after_jump);
        
        uint64_t runend_mask = mask_left(MEM_UNIT - current_shift - (!is_occupied(current_block*BLOCK_SIZE))); 
        // if slot 0 of block is not occupied, look runend[paj+1:end]
        // if slot 0 of block is occupied, look runend[paj:end]

        // === FIRST BLOCK WE JUMP INTO ===
        if (offset < BLOCK_SIZE){
            // jumped into the same block
            if (nb_runs == 0){
                return pos_after_jump < quotient ? quotient : get_next_quot(pos_after_jump); // changed <= to < (might bug, maybe check if block[0] is occupied)
            }

            select_val = bitselectasm(get_runend_word(current_block) & runend_mask, 
                                      nb_runs);            

            //cout << "select_val " << select_val << endl;
                                            
            if (select_val < get_shift_in_block(quotient)){ return quotient; }
        }
        else {
            // jumped w/ offset into further block
            if (nb_runs == 0){
                return get_next_quot(pos_after_jump);
            }

            current_block = get_block_id(pos_after_jump);
            select_val = bitselectasm(get_runend_word(current_block) & runend_mask, 
                                    nb_runs);
        }

        nb_runs -= bitrankasm(get_runend_word(current_block) & runend_mask, 
                            MEM_UNIT-1);

        // === OTHER BLOCKS TO FIND THE END OF dth RUN === 
        while (select_val == BLOCK_SIZE){     
                current_block = get_next_block_id(current_block);
                select_val = bitselectasm(get_runend_word(current_block), nb_runs);
                nb_runs -= bitrankasm(get_runend_word(current_block), MEM_UNIT-1);
        }

        return (current_block*BLOCK_SIZE + select_val + 1) % (1ULL << quotient_size); 
        //+1 because select_val is the end of the dth run before our quot's one
    }
} 

uint64_t Rsqf::get_runstart_shift0(uint64_t quotient, uint64_t paj, uint64_t offset, bool occ_bit){
    if (verbose){
        cout << "get_runstart_shift0 " << quotient << " " << paj << " " << offset << " " << occ_bit << endl; 
    } 
    //Used in get_runstart, NOT OPTI AT ALL
    //example of case where this func is needed :
    // run 1 : quotient50; start50; end70
    // run 2 : quotient64; start71; end80;  get_runstart(64)?
    // run 2' : quotient64; start71; end150;  get_runstart(64)?
    uint64_t nth_runend = 0;
    uint64_t current_block = get_block_id(quotient); 
    if (offset < BLOCK_SIZE){ //jump same block
        nth_runend = bitrankasm(get_runend_word(current_block), get_shift_in_block(paj));

        if (nth_runend-occ_bit <= 0){ return quotient; } //only one run in beginning of block, ours
        
        else { //more runs before the interesting one, we find the end of the last then next quot
            uint64_t last_run_runend_shift = bitselectasm(get_runend_word(current_block), nth_runend-occ_bit);
            return get_next_quot(current_block*64 + last_run_runend_shift); 
        }
        
    }

    //jumped into further block
    uint64_t final_block = get_block_id(paj);
    while (current_block != final_block){
        nth_runend += bitrankasm(get_runend_word(current_block), MEM_UNIT-1);
        current_block = get_next_block_id(current_block);
    }
    nth_runend += bitrankasm(get_runend_word(current_block), get_shift_in_block(paj));
    
    if (nth_runend-occ_bit <= 0){ return quotient; } //only one run until pos_after_jump, ours
    
    else { //other runs before the interesting one, we have to find in which block the last one ends up
        current_block = get_block_id(quotient);
        while ( nth_runend > BLOCK_SIZE || bitselectasm(get_runend_word(current_block), nth_runend-occ_bit) == BLOCK_SIZE ){
            nth_runend -= bitrankasm(get_runend_word(current_block), MEM_UNIT-1);
            current_block = get_next_block_id(current_block);
        }
        uint64_t last_run_runend_shift = bitselectasm(get_runend_word(current_block), nth_runend-occ_bit);
        return get_next_quot(current_block*64 + last_run_runend_shift); 
    }
}


pair<uint64_t,uint64_t> Rsqf::get_run_boundaries(uint64_t quotient){ //const
    //SUB OPTI
    assert(is_occupied(quotient));
    
    pair<uint64_t, uint64_t> boundaries;

    boundaries.first = get_runstart(quotient, 1);
    boundaries.second = get_runend(quotient).first;

    return boundaries;
}


void Rsqf::shift_bits_left_metadata(uint64_t quotient, uint64_t overflow_bit, uint64_t start_position, uint64_t end_position){
    if (verbose){
        cout << "shift_bits_left_metadata quotient " << quotient << " SP " << start_position << " EP " << end_position << endl;
    }
    
    assert(start_position < (1ULL << quotient_size));
    assert(end_position < (1ULL << quotient_size));
    
    // METHOD FOR INSERTION
    uint64_t current_block = get_block_id(quotient);
    uint64_t current_shift_in_block = get_shift_in_block(start_position);

    uint64_t start_block = get_block_id(start_position);
    
    uint64_t end_block = get_block_id(end_position);
    uint64_t end_shift_in_block = get_shift_in_block(end_position);
    //(end_position == first unused slot)

    uint64_t word_to_shift;
    uint64_t save_right;
    uint64_t to_shift;
    uint64_t next_block;
    uint64_t save_left;

    //OFFSET case quotient is in first slot (TEST_F(RsqfTest, offset2))
    if (get_shift_in_block(quotient) == 0) {
        set_offset_word(current_block, get_offset_word(current_block)+1);
    }

    //OFFSET case everyblock we cross until runstart 
    while (current_block != start_block){
        current_block = get_next_block_id(current_block);
        set_offset_word(current_block, get_offset_word(current_block)+1);
    } 
    
    // current_block == start_block
    //RUNEND case runstart at the end of filter (shift almost all filter)
    if ((current_block == end_block) && (start_position > end_position)) {

        word_to_shift = get_runend_word(current_block);
        save_right = word_to_shift & mask_right(current_shift_in_block);
        to_shift = shift_left(shift_right(word_to_shift,current_shift_in_block), (current_shift_in_block + 1));

        overflow_bit = shift_left(overflow_bit,current_shift_in_block);
        to_shift |= (save_right | overflow_bit);
        set_runend_word(current_block, to_shift); 
        overflow_bit = shift_right(word_to_shift,(MEM_UNIT - 1));

        current_block = get_next_block_id(current_block);
        current_shift_in_block = 0;
        set_offset_word(current_block, get_offset_word(current_block)+1);
    }


    while ( current_block != end_block ) {

        word_to_shift = get_runend_word(current_block);
        save_right = word_to_shift & mask_right(current_shift_in_block);
        to_shift = shift_left(shift_right(word_to_shift,current_shift_in_block), (current_shift_in_block + 1));

        overflow_bit = shift_left(overflow_bit,current_shift_in_block);
        to_shift |= (save_right | overflow_bit);
        set_runend_word(current_block, to_shift); 
        overflow_bit = shift_right(word_to_shift,(MEM_UNIT - 1));

        next_block = get_next_block_id(current_block);
        current_block = next_block;

        //OFFSET case everyblock we cross until runend 
        set_offset_word(current_block, get_offset_word(current_block)+1);
        current_shift_in_block = 0;
    }

    word_to_shift = get_runend_word(current_block);
    save_left = (word_to_shift & mask_left(MEM_UNIT-end_shift_in_block));
    save_right = word_to_shift & mask_right(current_shift_in_block);

    to_shift = (((word_to_shift & mask_right(end_shift_in_block)) & mask_left(MEM_UNIT - current_shift_in_block)) << 1);
    to_shift |= ((save_left | shift_left(overflow_bit,current_shift_in_block)) | save_right);

    set_runend_word(current_block, to_shift); 
}


void Rsqf::shift_runend_right(uint64_t start_shift, uint64_t end_shift, uint64_t block){
    uint64_t word_to_shift = get_runend_word(block);

    uint64_t save_right = word_to_shift & mask_right(start_shift); //start_shift in pos @5, we want to save the 5 first elem (0,1,2,3,4)
    uint64_t save_left = word_to_shift & mask_left(MEM_UNIT - end_shift - 1);

    word_to_shift &= mask_left(MEM_UNIT - start_shift);
    word_to_shift &= mask_right(end_shift+1);
    word_to_shift >>= 1;
    
    uint64_t overflow_bit = 0;
    if (end_shift == MEM_UNIT-1){
        overflow_bit = shift_left(get_runend_word(get_next_block_id(block)) & 1ULL, MEM_UNIT-1);
    }
    
    word_to_shift |= save_right | save_left | overflow_bit;
    set_runend_word(block, word_to_shift);
}

void Rsqf::shift_bits_right_metadata(uint64_t quotient, uint64_t start_position, uint64_t end_position){
    if (verbose){
        cout << "shift_bits_RIGHT_metadata quotient " << quotient << " SP " << start_position << " EP " << end_position << endl;
    }
    // METHOD FOR DELETION

    uint64_t current_block = get_block_id(quotient);
    uint64_t current_shift_in_block = get_shift_in_block(quotient);

    uint64_t start_block = get_block_id(start_position); 
    uint64_t end_block = get_block_id(end_position);
    uint64_t end_shift_in_block = get_shift_in_block(end_position);

    //OFFSET FROM QUOTIENT TO POS_ELEM
    if (current_shift_in_block == 0) {
        decrement_offset(current_block);
    }

    //OFFSET FROM QUOTIENT TO POS_ELEM
    while (current_block != start_block){
        current_block = get_next_block_id(current_block);
        decrement_offset(current_block);
    }
  
    current_shift_in_block = get_shift_in_block(start_position);

    if (current_shift_in_block == 0) {
        //unitary test TEST_F(RsqfTest, shift_bits_right_metadata)
        uint_fast64_t prev_block = get_prev_block_id(current_block);
        uint64_t overflow_bit = shift_left(get_runend_word(current_block) & 1ULL, MEM_UNIT-1);
        set_runend_word(prev_block, get_runend_word(prev_block) | overflow_bit);
    }


    //OFFSET+RUNEND CASE SHIFT ALMOST ALL FILTER
    if ((current_block == end_block) && (start_position > end_position)){
        shift_runend_right(current_shift_in_block, MEM_UNIT-1, current_block);

        current_block = get_next_block_id(current_block);
        current_shift_in_block = 0;
        decrement_offset(current_block);        
    }

    //OFFSET+RUNEND EVERYBLOCK UNTIL END
    while ( current_block != end_block ){
        shift_runend_right(current_shift_in_block, MEM_UNIT-1, current_block);

        current_block = get_next_block_id(current_block);
        decrement_offset(current_block);
        current_shift_in_block = 0;
    }

    //OFFSET+RUNEND CURRENT BLOCK == END BLOCK
    shift_runend_right(current_shift_in_block, end_shift_in_block, current_block);
}


void Rsqf::save_on_disk(const string& filename) { //remove 5 
    ofstream file(filename, ios::out | ios::binary);
    if (file.is_open()) {
        file.write(reinterpret_cast<const char*>(&this->quotient_size), sizeof(uint64_t));
        file.write(reinterpret_cast<const char*>(&this->remainder_size), sizeof(uint64_t));
        uint64_t num_words = (1ULL<<this->quotient_size) * (MET_UNIT + remainder_size) / MEM_UNIT;
        file.write(reinterpret_cast<const char*>(this->filter.data()), sizeof(uint64_t) * num_words);
        file.close();
    } else {
        cerr << "Unable to open file for writing: " << filename << endl;
        exit( EXIT_FAILURE );
    }
}

Rsqf Rsqf::load_from_disk(const string& filename){
    Rsqf qf;
    ifstream file(filename, ios::in | ios::binary);
    if (file.is_open()) {
        file.read(reinterpret_cast<char*>(&qf.quotient_size), sizeof(int64_t));
        file.read(reinterpret_cast<char*>(&qf.remainder_size), sizeof(int64_t));
        uint64_t num_words = (1ULL<<qf.quotient_size) * (MET_UNIT + qf.remainder_size) / MEM_UNIT;
        qf.filter.resize(num_words);
        file.read(reinterpret_cast<char*>(qf.filter.data()), sizeof(int64_t) * num_words);
        file.close();
    } else {
        cerr << "Unable to open file for reading: " << filename << endl;
        return 1;
    }
    return qf;
}
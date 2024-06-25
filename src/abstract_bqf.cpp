#include "abstract_bqf.hpp" 
#include <stdio.h>


using namespace std;

void Bqf::insert(string kmc_input){
    try {
        ifstream infile(kmc_input);

        if (!infile) {
            throw std::runtime_error("File not found: " + kmc_input);
        }

        string smer; 
        uint64_t count;

        //1st elem, check s == smer_size
        infile >> smer >> count;
        if (smer.size() == this->smer_size){
            this->insert(smer, count);
        } else {
            std::cerr << "BQF has been configured to welcome " << this->smer_size << "mers but trying to insert " << smer.size() << "mers, end of insertions" << std::endl;
            return;
        }
        

        while (infile >> smer >> count) {
            this->insert(smer, count);
        }

        infile.close();
    } catch (const std::exception &e) {
        // Gérez l'exception ici, par exemple, affichez un message d'erreur
        std::cerr << "Error: " << e.what() << std::endl;
    }
}   


void Bqf::insert(string kmer, uint64_t count){
    this->insert(kmer_to_hash(kmer, smer_size), count);
}

void Bqf::insert(uint64_t number, uint64_t count){
    if (elements_inside+1 == size_limit){
        if (verbose){
            cout << "RESIZING, nbElem: " << elements_inside << endl;
        }
        this->resize(1);    
    }

    //get quotient q and remainder r
    uint64_t quot = quotient(number);
    uint64_t rem = remainder(number);
    uint64_t rem_count = (rem << count_size) | insert_process_count(count); //PROCESS count
    //handles count > 2^c 

    if (verbose){
        cout << "[INSERT] quot " << quot << " from hash " << number << endl;
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
        pair<uint64_t,uint64_t> boundary = get_run_boundaries(quot);

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
        shift_left_and_set_circ(starting_position, fu_slot, rem_count);
    }
}

void Bqf::query(std::ifstream& infile, std::ofstream& outfile){
    std::string seq; 
    uint64_t i = 1;

    std::getline(infile, seq); 
    if (seq[0] == '>'){ //conventionnal fasta
        std::getline(infile, seq);
        outfile << "Sequence" << i << " : " << this->query(seq) << "\n";
        while (std::getline(infile, seq)) { // skip first header line
            std::getline(infile, seq);

            outfile << "Sequence" << i++ << " : " << this->query(seq) << "\n";
        }
    } else { //1 seq / line
        outfile << "Sequence" << i << " : " << this->query(seq) << "\n";
        while (std::getline(infile, seq)) {
            outfile << "Sequence" << i++ << " : " << this->query(seq) << "\n";
        }
    }

    infile.close();
    outfile.close();
}



result_query Bqf::query(string seq){
    int s = this->smer_size;
    int k = this->kmer_size;
    int n = seq.length();
    
    if (k == s && s == n) { 
        uint64_t res = this->query(bfc_hash_64(flip(canonical(flip(encode(seq), 2*s), 2*s), 2*s), mask_right(s*2)));
        return result_query {(int)res, (int)res, (float)res, (float)(res!=0)};
    }
    int z = k-s;
    int last_smers_abundances[z+1];
    int* kmer_abundance;
    int nb_presence = 0;
    int avg = 0;
    int minimum = numeric_limits<int>::max();
    int maximum = 0;

    uint64_t current_smer = 0;
    
    //build current_smer (s first chars)
    for (auto i = 0; i < s-1; i++){
        current_smer <<= 2;
        current_smer |= nucl_encode(seq[i]);
    } 

    //1st kmer (s+z first chars), skipped if k==s
    for (auto i = s-1; i < k-1; i++){
        current_smer <<= 2;
        current_smer = (current_smer | nucl_encode(seq[i])) & mask_right(2*s);

        last_smers_abundances[i-(s-1)] = this->query(bfc_hash_64(flip(canonical(current_smer, 2*s), 2*s), mask_right(s*2)));
    }


    //all kmers
    for (auto i = k-1; i < n; i++){
        current_smer <<= 2;
        current_smer = (current_smer | nucl_encode(seq[i])) & mask_right(2*s);

        last_smers_abundances[(i-s+1)%(z+1)] = this->query(bfc_hash_64(flip(canonical(current_smer, 2*s), 2*s), mask_right(s*2)));

        kmer_abundance = min_element(last_smers_abundances, last_smers_abundances+z+1);
        if (*kmer_abundance == 0){
            minimum = 0;
        } else {
            minimum = std::min(minimum, *kmer_abundance);
            maximum = std::max(maximum, *kmer_abundance);
            avg = avg + *kmer_abundance;
            nb_presence ++;
        }
    }

    return result_query {minimum, maximum, (float)(avg / (n-k+1)), (float)nb_presence/(n-k+1)};
}

uint64_t Bqf::query(uint64_t number){
    if (elements_inside == 0) return 0;
    uint64_t quot = quotient(number);
    uint64_t rem = remainder(number);
    if (!is_occupied(quot)) return 0;

    pair<uint64_t,uint64_t> boundary = get_run_boundaries(quot);

    // TODO:
    // OPTIMIZE TO LOG LOG SEARCH ?

    uint64_t position = boundary.first;

    while(position != boundary.second){
        uint64_t remainder_in_filter = get_remainder(position);

        if (remainder_in_filter == rem) return query_process_count(get_remainder(position, true) & mask_right(count_size));
        else if (remainder_in_filter > rem) return 0;
        position = get_next_quot(position);
    }

    uint64_t remainder_in_filter = get_remainder(boundary.second); 
    if (remainder_in_filter == rem) return query_process_count(get_remainder(position, true) & mask_right(count_size));

    return 0;
}

std::map<uint64_t, uint64_t> Bqf::enumerate(){
    std::map<uint64_t, uint64_t> finalSet;
    uint64_t curr_occ;
    
    std::pair<uint64_t, uint64_t> bounds;
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
                    finalSet[number] = query_process_count(get_remainder(cursor, true) & mask_right(count_size));
                    cursor = get_next_quot(cursor);
                }

                number = rebuild_number(quotient, get_remainder(cursor), quotient_size);
                finalSet[number] = query_process_count(get_remainder(cursor, true) & mask_right(count_size));
            }

            curr_occ >>= 1ULL; //next bit of occupied vector
        }
    }

    assert(finalSet.size() == this->elements_inside);

    return finalSet;
}

void Bqf::resize(uint n){
    if(n == 0) return;
    
    assert(n <= this->remainder_size);

    const uint64_t new_quotient_size = this->quotient_size + n;
    const uint64_t new_remainder_size = this->remainder_size - n;

    const uint64_t num_quots = 1ULL << new_quotient_size; 
    const uint64_t num_of_words = num_quots * (MET_UNIT + new_remainder_size) / MEM_UNIT;

    // In machine words
    const uint64_t new_number_blocks = std::ceil(num_quots / BLOCK_SIZE);
    std::vector<uint64_t> new_filter = std::vector<uint64_t>(num_of_words);

    uint64_t current_block;
    uint64_t curr_occ;

    std::pair<uint64_t, uint64_t> bounds;
    uint64_t cursor;
    uint64_t endCursor;

    uint64_t quotient;
    uint64_t new_quotient;
    uint64_t added_quotient_bits;
    uint64_t remainder;
    uint64_t count;

    uint64_t new_block;
    uint64_t pos_in_block;
    uint64_t pos;

    // init the offsets
    const uint offsets_size = (1 << n);
    std::pair<uint32_t, bool> offsets[offsets_size];
    std::fill(offsets, offsets + offsets_size, std::make_pair(0, false));

    // starting position
    const uint64_t start = first_unused_slot(0);
    const uint64_t block_start = get_block_id(start);

    for(uint block = 0; block <= this->number_blocks; ++block){
        // getting the true block
        current_block = block_start + block;
        if (current_block >= this->number_blocks){
            current_block -= this->number_blocks;
        }

        // shifting the offsets when going back to 0
        if(block != 0 && current_block == 0){
            pos = offsets[offsets_size - 1].first;
            for(uint i = offsets_size - 1; i > 0; --i){
                offsets[i].first = offsets[i - 1].first;
            }
            offsets[0].first = pos;
        }

        // the occupied word of the block
        curr_occ = this->get_occupied_word(current_block);

        // skip if whole block is empty
        if (curr_occ == 0) {
            for(uint i = 0; i < offsets_size; ++i){
                // setting the offsets
                new_block = ((i << this->quotient_size) / BLOCK_SIZE) | current_block;
                pos = (new_block *(MET_UNIT+new_remainder_size)) + OFF_POS;
                new_filter[pos] = offsets[i].first;

                // substracting / resetting values
                if(offsets[i].first >= BLOCK_SIZE){
                    offsets[i].first -= BLOCK_SIZE;
                } else {
                    offsets[i].first = 0;
                }
                offsets[i].second = false;
            }
            continue;
        }

        // inserting the elements in the new filter
        for (uint i=0; i<BLOCK_SIZE; i++){
            // start
            if(block == 0 && i == 0){
                pos = get_shift_in_block(start);
                i = pos;
                curr_occ >>= pos;
            }

            quotient = current_block*BLOCK_SIZE + i;

            // end
            if(block == this->number_blocks && quotient == start){
                break;
            }
            
            if (curr_occ & 1ULL){ //occupied
                bounds = this->get_run_boundaries(quotient);
                cursor = bounds.first;
                endCursor = this->get_next_quot(bounds.second);

                while (cursor != endCursor){ // every remainder of the run
                    remainder = this->get_remainder(cursor, true);
                    count = remainder & mask_right(this->count_size);
                    remainder >>= this->count_size;
                    added_quotient_bits = remainder & mask_right(n);

                    assert(added_quotient_bits >= 0 && added_quotient_bits < offsets_size);

                    // new remainder and new quotient
                    remainder = ((remainder >> n) << this->count_size) | count;
                    new_quotient = ((added_quotient_bits << this->quotient_size) | quotient) + offsets[added_quotient_bits].first;
                    if (new_quotient >= num_quots){
                        new_quotient -= num_quots;
                    }

                    // set the remainder
                    new_block = get_block_id(new_quotient);
                    pos_in_block = get_shift_in_block(new_quotient);
                    pos = new_block * ((MET_UNIT+new_remainder_size)*BLOCK_SIZE) + MET_UNIT*BLOCK_SIZE + pos_in_block*new_remainder_size;
                    set_bits(new_filter, pos, remainder, new_remainder_size);

                    offsets[added_quotient_bits].first++;
                    offsets[added_quotient_bits].second = true;

                    cursor = this->get_next_quot(cursor);
                }


                for(uint i = 0; i < offsets_size; ++i){
                    if (offsets[i].second) {
                        // setting the occupied
                        new_quotient = (i << this->quotient_size) | quotient;
                        new_block = get_block_id(new_quotient);
                        pos = (new_block *(MET_UNIT+new_remainder_size)) + OCC_POS;
                        new_filter[pos] |= (1ULL << new_quotient);

                        // setting the runend
                        new_quotient += offsets[i].first - 1;
                        if (new_quotient >= num_quots){
                            new_quotient -= num_quots;
                        }
                        new_block = get_block_id(new_quotient);
                        pos = (new_block *(MET_UNIT+new_remainder_size)) + RUN_POS;
                        new_filter[pos] |= (1ULL << new_quotient);
                    }
                }
            }

            // setting the offset 
            if (i == 0){
                for(uint i = 0; i < offsets_size; ++i){
                    new_block = ((i << this->quotient_size) / BLOCK_SIZE) | current_block;
                    pos = (new_block *(MET_UNIT+new_remainder_size)) + OFF_POS;
                    new_filter[pos] = offsets[i].first;
                }
            }

            // substracting and resetting the offsets
            for(uint i = 0; i < offsets_size; ++i){
                if(offsets[i].first > 0){
                    offsets[i].first--;
                }
                offsets[i].second = false;
            }

            curr_occ >>= 1ULL; //next bit of occupied vector
        }
    }

    this->quotient_size = new_quotient_size;
    this->remainder_size = new_remainder_size;

    this->size_limit = num_quots * 0.95;

    this->number_blocks = new_number_blocks;
    
    this->filter.swap(new_filter);
}

uint64_t Bqf::get_remainder(uint64_t position, bool w_counter ){ //default=false
    uint64_t block = get_block_id(position);
    uint64_t pos_in_block = get_shift_in_block(position);
    uint64_t pos = block * ((MET_UNIT+remainder_size)*BLOCK_SIZE) + MET_UNIT*BLOCK_SIZE + pos_in_block*remainder_size; 

    if (w_counter) return get_bits(filter, pos, remainder_size);
    else return get_bits(filter, pos, remainder_size) >> count_size;
}


uint64_t Bqf::find_quotient_given_memory(uint64_t max_memory, uint64_t count_size){
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


void Bqf::save_on_disk(const std::string& filename) { 
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (file.is_open()) {
        file.write(reinterpret_cast<const char*>(&this->quotient_size), sizeof(uint64_t));
        file.write(reinterpret_cast<const char*>(&this->remainder_size), sizeof(uint64_t));
        file.write(reinterpret_cast<const char*>(&this->count_size), sizeof(uint64_t));
        file.write(reinterpret_cast<const char*>(&this->kmer_size), sizeof(uint64_t));
        file.write(reinterpret_cast<const char*>(&this->smer_size), sizeof(uint64_t));
        file.write(reinterpret_cast<const char*>(&this->size_limit), sizeof(uint64_t));
        file.write(reinterpret_cast<const char*>(&this->number_blocks), sizeof(uint64_t));
        file.write(reinterpret_cast<const char*>(&this->elements_inside), sizeof(uint64_t));
        uint64_t num_words = (1ULL<<this->quotient_size) * (MET_UNIT + remainder_size) / MEM_UNIT;
        file.write(reinterpret_cast<const char*>(this->filter.data()), sizeof(uint64_t) * num_words);
        file.close();
    } else {
        std::cerr << "Unable to open file for writing: " << filename << std::endl;
    }
}

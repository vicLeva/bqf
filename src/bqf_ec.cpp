#include "bqf_ec.hpp" 

using namespace std;

Bqf_ec::Bqf_ec(){}

Bqf_ec::Bqf_ec(uint64_t q_size, uint64_t c_size, uint64_t k, uint64_t z, bool verb) :
    Bqf(c_size, k, k-z, q_size, verb)
{
    const uint64_t hash_size = 2*smer_size;   

    remainder_size = hash_size - q_size + c_size;
    
    if (verbose){ 
        cout << "q " << quotient_size << " r " << (hash_size - q_size) << " remainder_size " << remainder_size << " count_size " << count_size << endl; 
    }
}


Bqf_ec::Bqf_ec(uint64_t max_memory, uint64_t c_size, bool verb):
    Bqf(max_memory, c_size, verbose){}



bool Bqf_ec::remove(string kmer, uint64_t count){
    return this->remove(kmer_to_hash(kmer, smer_size), count);
}

bool Bqf_ec::remove(uint64_t number, uint64_t count){
    if (elements_inside == 0) return 0;
    //get quotient q and remainder r
    const uint64_t quot = quotient(number);
    const uint64_t rem = remainder(number);

    if (verbose){
        cout << "[REMOVE] quot " << quot << endl;
        cout << "[REMOVE] rem " << rem << endl;
    }

    if (is_occupied(quot) == false) return 0;

    const pair<uint64_t,uint64_t> boundary = get_run_boundaries(quot);

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
    const uint64_t end_slot = first_unshiftable_slot(quot);

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




void Bqf_ec::add_to_counter(uint64_t position, uint64_t remainder_w_count){
    const uint64_t old_rem = get_remainder(position, true);
    uint64_t sum = (old_rem & mask_right(count_size)) + (remainder_w_count & mask_right(count_size));
    if (! (sum < 1ULL << (count_size))){
        sum = (1ULL << (count_size)) - 1;
    }  
    
    sum |= old_rem & mask_left(MEM_UNIT - count_size);

    set_bits(filter, 
            get_remainder_word_position(position) * BLOCK_SIZE + get_remainder_shift_position(position), 
            sum, 
            remainder_size);
}

void Bqf_ec::sub_to_counter(uint64_t position, uint64_t count){
    const uint64_t old_rem = get_remainder(position, true);

    uint64_t sub = (old_rem & mask_right(count_size)) - count;

    sub |= old_rem & mask_left(MEM_UNIT - count_size);

    set_bits(filter, 
            get_remainder_word_position(position) * BLOCK_SIZE + get_remainder_shift_position(position), 
            sub, 
            remainder_size);
}

uint64_t Bqf_ec::insert_process_count(uint64_t count){
    return (count < (1ULL << count_size) ? count : (1ULL << count_size)-1);
}

uint64_t Bqf_ec::query_process_count(uint64_t count){ //nothing to do, used for bqf_oom
    return count;
}

Bqf_ec Bqf_ec::load_from_disk(const string& filename){
    Bqf_ec qf;
    ifstream file(filename, ios::in | ios::binary);
    if (file.is_open()) {
        file.read(reinterpret_cast<char*>(&qf.quotient_size), sizeof(uint64_t));
        file.read(reinterpret_cast<char*>(&qf.remainder_size), sizeof(uint64_t));
        file.read(reinterpret_cast<char*>(&qf.count_size), sizeof(uint64_t));
        file.read(reinterpret_cast<char*>(&qf.kmer_size), sizeof(uint64_t));
        file.read(reinterpret_cast<char*>(&qf.smer_size), sizeof(uint64_t));
        file.read(reinterpret_cast<char*>(&qf.size_limit), sizeof(uint64_t));
        file.read(reinterpret_cast<char*>(&qf.number_blocks), sizeof(uint64_t));
        file.read(reinterpret_cast<char*>(&qf.elements_inside), sizeof(uint64_t));
        const uint64_t num_words = (1ULL<<qf.quotient_size) * (MET_UNIT + qf.remainder_size) / MEM_UNIT;
        qf.filter.resize(num_words);
        file.read(reinterpret_cast<char*>(qf.filter.data()), sizeof(uint64_t) * num_words);
        file.close();

        qf.verbose = false;
    } else {
        cerr << "Unable to open file for reading: " << filename << endl;
    }
    return qf;
}



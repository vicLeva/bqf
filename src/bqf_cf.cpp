#include "bqf_cf.hpp"

using namespace std;

/*  
    ================================================================
    CONSTRUCTORS
    ================================================================
*/ 
Bqf_cf::Bqf_cf(uint64_t q_size, uint64_t k, uint64_t z, bool verb=false) :
    Bqf_ec(q_size, 1, k, z, verb) {};

bool Bqf_cf::add_to_counter(uint64_t position){
    const uint64_t old_rem = get_remainder(position, true);
    if (!(old_rem & 1ULL)){ // 2nd occurence
        set_bits(filter, 
            get_remainder_word_position(position) * BLOCK_SIZE + get_remainder_shift_position(position), 
            old_rem | 1ULL, 
            remainder_size);
        return true;
    }
    return false;
}


/*  
    ================================================================
    HGH-LEVEL METHODS
    ================================================================
*/ 
bool Bqf_cf::is_second_insert(uint64_t number){
    if (elements_inside+1 == size_limit){
        if (verbose){
            cout << "RESIZING, nbElem: " << elements_inside << endl;
        }
        this->resize(1);    
    }

    //get quotient q and remainder r
    const uint64_t quot = quotient(number);
    const uint64_t rem = remainder(number);
    const uint64_t rem_count = rem<<1;
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
        shift_left_and_set_circ(starting_position, fu_slot, rem_count);
        return false;
    }
    // IF THE QUOTIENT HAS BEEN USED BEFORE
    // GET POSITION WHERE TO INSERT TO (BASED ON VALUE) IN THE RUN (INCREASING ORDER)
    else{
        if (verbose){
            cout << "occupied" << endl;
        }

        //getting boundaries of the run
        const pair<uint64_t,uint64_t> boundary = get_run_boundaries(quot);

        if (verbose){
            cout << "boundaries " << boundary.first << " || " << boundary.second << endl;
        }

        // nb of quotients
        const uint64_t quots = (1ULL << this->quotient_size);

        // dichotomous search and insertion
        uint64_t left = boundary.first;
        if (left < quot)   
            left += quots;

        uint64_t right = boundary.second;
        if (right < quot) 
            right += quots;

        uint64_t middle = ceil((left + right) / 2);
        uint64_t position = middle;
        if (position >= quots)
            position -= quots;
        
        uint64_t remainder_in_filter;

        assert(left <= right);

        while (left <= right) {
            middle = ceil((left + right) / 2);
            position = middle;
            if (position >= quots)
                position -= quots;
            remainder_in_filter = get_remainder(position);

            if (remainder_in_filter == rem)
                return add_to_counter(position);
            else if (left == right){
                if (remainder_in_filter < rem)
                    position = get_next_quot(position);
                break;
            }
            else if (remainder_in_filter > rem)
                right = middle;
            else
                left = middle + 1;
            
        }

        shift_bits_left_metadata(quot, 0, boundary.first, fu_slot);
        // SHIFT EVERYTHING RIGHT AND INSERTING THE NEW REMAINDER
        elements_inside++;
        shift_left_and_set_circ(position, fu_slot, rem_count);
        return false;
    }
}

void Bqf_cf::is_second_insert(string kmer, ofstream& output){
    if (this->is_second_insert(kmer_to_hash(kmer, smer_size))) {
        output << kmer << endl;
    }
}

void Bqf_cf::insert_and_filter(string kmc_input, string output) {
    try {
        ifstream infile(kmc_input);
        ofstream outfile(output);

        if (!infile) {
            throw std::runtime_error("File not found: " + kmc_input);
        }

        string smer; 

        //1st elem, check s == smer_size
        infile >> smer;
        if (smer.size() == this->smer_size){
            this->is_second_insert(smer, outfile);
        } else {
            std::cerr << "BQF has been configured to welcome " << this->smer_size << "mers but trying to insert " << smer.size() << "mers, end of insertions" << std::endl;
            return;
        }
        

        while (infile >> smer) {
            this->is_second_insert(smer, outfile);
        }

        infile.close();
    } catch (const std::exception &e) {
        // Gérez l'exception ici, par exemple, affichez un message d'erreur
        std::cerr << "Error: " << e.what() << std::endl;
    }
}
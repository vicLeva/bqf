#include "bqf_cf.hpp"
#pragma push_macro("BLOCK_SIZE")
#undef BLOCK_SIZE
#include "FastxParser.hpp"
#pragma pop_macro ("BLOCK_SIZE")



using namespace std;

/*  
    ================================================================
    CONSTRUCTOR
    ================================================================
*/ 
Bqf_cf::Bqf_cf(uint64_t q_size, uint64_t k, bool verb) :
    Bqf_ec(q_size, 1, k, 0, verb) {};



bool Bqf_cf::is_second_add_to_counter(uint64_t position){
    const uint64_t old_rem = get_remainder(position, true);
    if (verbose) {
        cout << "[ADD] to old_rem " << old_rem << endl;
    }
    bool is_second = !(old_rem & 1ULL);
    if (is_second){
        //flipping count bit to one
        uint64_t pos = position * remainder_size + ((1 + (position/MEM_UNIT))*MEM_UNIT)*MET_UNIT;
        filter[pos/MEM_UNIT] |= (1ULL << (pos&(MEM_UNIT - 1)));
        /* set_bits(filter, 
            get_remainder_word_position(position) * BLOCK_SIZE + get_remainder_shift_position(position), 
            old_rem | 1ULL, 
            remainder_size); */
    }
    return is_second;
}


/*  
    ================================================================
    HIGH-LEVEL METHODS
    ================================================================
*/ 
void Bqf_cf::filter_fastx_file(std::vector<std::string> files, std::string output) {
    string to_write;
    
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(files, 1, 1);
    parser.start();
    auto rg = parser.getReadGroup();
    while (parser.refill(rg)) {
        for (auto& rp : rg) {
            insert_from_sequence(rp.seq, to_write);
        }
    }
    parser.stop();
    ofstream outfile(output);
    if (!outfile.is_open()) {
        throw std::runtime_error("Could not open file " + output);
    }
    outfile << counter << endl;
    outfile << to_write;
    outfile.close();
        
}

void Bqf_cf::insert_from_sequence (std::string sequence, string& to_write) {
    if (verbose){
        cout << "[INSERT] sequence " << sequence << endl;
    }
    uint64_t lgth = sequence.length();
    /*the kmer and its revcomp are created character per character and left_to_compute indicates 
    the nb of characters left to take into account to create a kmer*/ 
    uint64_t kmer = 0;
    uint64_t revcomp = 0;
    uint64_t mask = mask_right(2*kmer_size);
    uint64_t left_to_compute = kmer_size;
    for (uint64_t i = 0; i < lgth; i++) {
        if (is_valid(sequence[i])) {
            const uint64_t encoded = ((sequence[i] >> 1) & 0b11); //quickly encodes the character
            kmer <<= 2;
            kmer |= encoded;
            kmer &= mask;

            revcomp >>= 2;
            revcomp |= ((0b10 ^ encoded) << (2 * (kmer_size - 1)));

            const uint64_t canon  = (kmer < revcomp) ? kmer : revcomp;


            left_to_compute -= left_to_compute ? 1 : 0;
            if (!left_to_compute) {
                is_second_insert(canon, to_write);
            }
        }
        else {
            kmer = 0;
            revcomp = 0;
            left_to_compute = kmer_size;
        }
        
    }
}

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

        pair<uint64_t, bool> pos_and_found = find_insert_position(boundary, rem);
        uint64_t position = pos_and_found.first;

        if (pos_and_found.second) {
            return is_second_add_to_counter(position);
        }
        shift_bits_left_metadata(quot, 0, boundary.first, fu_slot);
        // SHIFT EVERYTHING RIGHT AND INSERT NEW REMAINDER
        elements_inside++;
        shift_left_and_set_circ(position, fu_slot, rem_count);

        return false;
    }
}

void Bqf_cf::is_second_insert(string kmer, ofstream& output){
    if (this->is_second_insert(kmer_to_hash(kmer, kmer_size))) {
        output << kmer << endl;
    }
}

void Bqf_cf::is_second_insert(uint64_t coded_kmer, string &to_write){
    if (this->is_second_insert(kmer_to_hash(coded_kmer, kmer_size))) {
        counter++;
        char rev[4] = {'A', 'C', 'T', 'G'};
        uint64_t mask = mask_right(2);

        for (size_t i=0; i<kmer_size; i++){
            to_write.push_back(rev[coded_kmer>>(2*(kmer_size - 1)) & mask]);
            coded_kmer <<= 2;
        }
        to_write.push_back('\n');
    }
}


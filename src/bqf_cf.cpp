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
Bqf_cf::Bqf_cf(uint64_t q_size, uint64_t k, output_mode_t mode, bool verb) :
    Bqf_ec(q_size, 1, k, 0, verb), mode(mode) {};

Bqf_cf::Bqf_cf(uint64_t max_memory, output_mode_t mode, bool verb) :
    Bqf_ec(max_memory, 1, verb), mode(mode) {};


bool Bqf_cf::is_second_add_to_counter(uint64_t position){
    const uint64_t old_rem = get_remainder(position, true);
    if (verbose) {
        cout << "[ADD] to old_rem " << old_rem << endl;
    }
    bool is_second = !(old_rem & 1ULL);
    if (is_second){
        //flipping count bit to one
        set_bits(filter, 
            get_remainder_word_position(position) * BLOCK_SIZE + get_remainder_shift_position(position), 
            old_rem | 1ULL, 
            remainder_size);
    }
    return is_second;
}


/*  
    ================================================================
    HIGH-LEVEL METHODS
    ================================================================
*/ 
void Bqf_cf::filter_fastx_file(std::vector<std::string> files, std::string output) {
    
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(files, 1, 1);
    parser.start();
    auto rg = parser.getReadGroup();
    while (parser.refill(rg)) {
        for (auto& rp : rg) {
            insert_from_sequence(rp.seq);
        }
    }
    parser.stop();
    ofstream outfile(output, ios::binary);
    if (!outfile.is_open()) {
        throw std::runtime_error("Could not open file " + output);
    }
    switch (mode) {
        case text :
            outfile << counter << endl;
            outfile << str_buffer;
            break;
        case binary:
            outfile.write(reinterpret_cast<const char*>(&counter), sizeof(uint64_t));
            outfile.write(reinterpret_cast<const char*>(bin_buffer.data()), counter*sizeof(uint64_t));
            break;
        case stream:
            std::cout << counter;
            for (auto& kmer : bin_buffer) {
                std::cout << kmer;
            }
            break;
    }
    outfile.close();
        
}

void Bqf_cf::insert_from_sequence (std::string sequence) {
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
                insert_kmer(canon);
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


void Bqf_cf::insert_kmer(uint64_t coded_kmer){
    if (this->is_second_insert(kmer_to_hash(coded_kmer, kmer_size))) {
        counter++;
        char rev[4] = {'A', 'C', 'T', 'G'};
        uint64_t mask = mask_right(2);

        switch (mode) {
        case text :
            for (size_t i=0; i<kmer_size; i++){
                str_buffer.push_back(rev[coded_kmer>>(2*(kmer_size - 1)) & mask]);
                coded_kmer <<= 2;
            }
            str_buffer.push_back('\n');
            break;
        default :
            bin_buffer.push_back(coded_kmer);
            break;  
        }
    }
}

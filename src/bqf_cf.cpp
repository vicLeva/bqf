#include "bqf_cf.hpp"
#pragma push_macro("BLOCK_SIZE")
#undef BLOCK_SIZE
#include "FastxParser.hpp"
#include "FastxParserThreadUtils.hpp"
#include "blockingconcurrentqueue.h"
#include "concurrentqueue.h"
#include "kseq++.hpp"
#include "lightweightsemaphore.h"
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
    HIGH-LEVEL METHODS
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

        pair<uint64_t, bool> pos_and_found = find_insert_position(boundary, quot, rem);
        uint64_t position = pos_and_found.first;

        if (pos_and_found.second) {
            return is_second_add_to_counter(position);
        }
        shift_bits_left_metadata(quot, 0, boundary.first, fu_slot);
        // SHIFT EVERYTHING RIGHT AND INSERTING THE NEW REMAINDER
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

void Bqf_cf::insert_from_file_and_filter(string input, string output) {
    try {
        ifstream infile(input);
        ofstream outfile(output, ios::binary);

        if (!infile) {
            throw std::runtime_error("File not found: " + input);
        }

        if (!outfile.is_open()) {
            throw std::runtime_error("Could not open file " + output);
        }

        string kmer; 

        //1st elem, check k == kmer_size
        infile >> kmer;
        if (kmer.size() == this->kmer_size){
            this->is_second_insert(kmer, outfile);
        } else {
            std::cerr << "BQF has been configured to welcome " << this->kmer_size << "mers but trying to insert " << kmer.size() << "mers, end of insertions" << std::endl;
            return;
        }
        

        while (infile >> kmer) {
            this->is_second_insert(kmer, outfile);
        }

        infile.close();
        outfile.close();
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void Bqf_cf::filter_fastx_file(std::vector<std::string> files, std::string output) {
    try {
        fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(files, 1, 1);
        parser.start();
        auto rg = parser.getReadGroup();
        while (parser.refill(rg)) {
            for (auto& rp : rg) {
                insert_from_sequence(rp.seq, output);
            }
        }
        parser.stop();
        }
    catch (const std::exception &e) {
        std::cerr << "Error :" << e.what() << std::endl;
    }
}

void Bqf_cf::insert_from_sequence (std::string sequence, std::string output) {
    /* uint64_t lgth = sequence.length();
    uint64_t kmer = 0;
    uint64_t revcomp = 0;
    uint64_t mask = mask_right(2*kmer_size);
    for (uint64_t i = 0; i < lgth - kmer_size; i++) {
        if (is_valid(sequence[i])) {

        }
        const uint64_t encoded = ((sequence[i] >> 1) & 0b11); //permet d'encoder rapidement un char selon son code ascii
        kmer <<= 2;
        kmer |= encoded;
        kmer &= mask;

        revcomp >>= 2;
        revcomp |= ( (0x2 ^ encoded) << (2 * (k - 1)));

        const uint64_t canon  = (kmer < revcomp) ? kmer : revcomp;
    } */
}
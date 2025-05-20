#ifndef BQF_CF
#define BQF_CF

#include "bqf_ec.hpp"

class Bqf_cf : public Bqf_ec {

public:
    Bqf_cf();
    /** 
     * \brief Constructor that instantiates a BackpackCQF from quotient and remainder sizes
     * \param q_size The desired size of quotient, will induce filter's size
     * \param k The length of a k_mer
     * \param z The length difference between a k-mer and the inserted s-mer
     * \param verb to print on-going operations in stdout (default: false)
     */
    Bqf_cf(uint64_t q_size, uint64_t c_size, uint64_t k, uint64_t z, bool verb=false);
};

#endif
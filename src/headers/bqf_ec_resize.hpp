#ifndef BQF_EC_RESIZE_HPP
#define BQF_EC_RESIZE_HPP

#include "bqf_ec.hpp"


class Bqf_ec_resize : public Bqf_ec{
    public:
    Bqf_ec_resize(uint64_t q_size, uint64_t c_size, uint64_t k, uint64_t z, bool verb=false);
    void resize(int n) override;
};

#endif
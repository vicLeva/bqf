#include "rsqf.hpp"
#include "abstract_bqf.hpp" 
#include "bqf_ec.hpp"

#include <chrono>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>

int main() {
    std::cout << "acquiring memory..." << std::endl;

    uint64_t* ptr = (uint64_t*) calloc(124000000, sizeof(uint64_t));

    if (ptr == NULL){
        std::cout << "not enough memory" << std::endl;
        exit(1);
    }
    
    std::cout << "memory aquired (<1G)" << std::endl;

    free(ptr);
    
    std::cout << "memory freed" << std::endl;
    std::cout << "acquiring more memory..." << std::endl;

    uint64_t* ptr2 = (uint64_t*) calloc(124000000, sizeof(uint64_t));

    if (ptr2 == NULL){
        std::cout << "not enough memory" << std::endl;
        exit(1);
    }
    
    std::cout << "memory aquired (>1G)" << std::endl;

    free(ptr2);

    std::cout << "memory freed, no problemo" << std::endl;
}

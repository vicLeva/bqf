#include "rsqf.hpp"
#include "abstract_bqf.hpp" 
#include "bqf_ec.hpp"

#include <chrono>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>

void mock_resize(Bqf* bqf, int n){
    std::map<uint64_t, uint64_t> inserted_elements = bqf->enumerate();

    bqf->quotient_size += n;
    bqf->remainder_size -= n;
    
    uint64_t num_quots = 1ULL << bqf->quotient_size; 
    uint64_t num_of_words = num_quots * (MET_UNIT + bqf->remainder_size) / MEM_UNIT; 

    bqf->size_limit = num_quots * 0.95;

    // In machine words
    bqf->number_blocks = std::ceil(num_quots / BLOCK_SIZE);

    bqf->filter = std::vector<uint64_t>(num_of_words);

    bqf->elements_inside = 0;

    for (auto const& elem : inserted_elements){
        bqf->insert(elem.first, elem.second);
    }
}

// std::string generateRandomKMer(int kmer_size) {
//     const char alphabet[] = "ACGT";
//     const int alphabetSize = sizeof(alphabet) - 1;
//     static std::random_device rd;
//     static std::mt19937 gen(rd());
//     static std::uniform_int_distribution<> dis(0, alphabetSize - 1);

//     std::stringstream randomKmer;
//     for (int i = 0; i < kmer_size; ++i) {
//         randomKmer << alphabet[dis(gen)];
//     }

//     return randomKmer.str();
// }

// struct Perf
// {
//     uint64_t q_size;
//     uint64_t inserted_elements;
//     double mock_time;
//     double resize_time;
// };

// Perf time_for_n_insertions(Bqf_ec* mock, Bqf_ec* resize){
//     // almost full
//     // while(resize->elements_inside < resize->size_limit - 1) {
//     //     std::string kmer = generateRandomKMer(resize->kmer_size);
//     //     resize->insert(kmer, 1);
//     //     mock->insert(kmer, 1);
//     // }

//     // calculating times
//     auto mockStart = std::chrono::high_resolution_clock::now();
//     mock_resize(mock, 1);
//     double mockTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - mockStart).count();
    
//     auto resizeStart = std::chrono::high_resolution_clock::now();
//     resize->resize(1);
//     double resizeTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - resizeStart).count();

//     return Perf { resize->quotient_size, resize->elements_inside, mockTime, resizeTime};
// }

int main() {
    const uint64_t q = 8;
    const uint64_t c = 5;
    const uint64_t k = 32;
    const uint64_t z = 32 - 19;
    Bqf_ec mock = Bqf_ec(q, c, k, z, false);
    Bqf_ec resize = Bqf_ec(q, c, k, z, false);
    const std::string test_file = "/scratch/vlevallois/data/7_1_counted_k19.txt";

    try {
        std::ifstream infile(test_file);

        if (!infile) {
            throw std::runtime_error("File not found: " + test_file);
        }

        std::string smer; 
        uint64_t count;
        infile >> smer >> count;
        
        if (smer.size() == resize.smer_size){
            mock.insert(smer, count);
            resize.insert(smer, count);
        } else {
            std::cerr << "BQF has been configured to welcome " << resize.smer_size << "mers but trying to insert " << smer.size() << "mers, end of insertions" << std::endl;
            exit(1);
        }
        
        std::cout << "Quotient size, Number of inserted elements, Old time, new time" << std::endl;

        while (infile >> smer >> count) {
            mock.insert(smer, count);
            resize.insert(smer, count);
            if(resize.elements_inside == resize.size_limit - 1){
                // calculating times
                auto mockStart = std::chrono::high_resolution_clock::now();
                mock_resize(&mock, 1);
                double mockTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - mockStart).count();
                
                auto resizeStart = std::chrono::high_resolution_clock::now();
                resize.resize(1);
                double resizeTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - resizeStart).count();

                for (uint i = 0; i < mock.filter.size(); ++i){
                    if (mock.filter[i] != resize.filter[i]) {
                        std::cerr << "mock and resize have diferent values !!" << std::endl;
                        exit(1);
                    }
                }

                std::cout << (resize.quotient_size - 1) << ',' << (resize.elements_inside + 1) << ',' << mockTime << ',' << resizeTime << std::endl;
            }
            
        }

        infile.close();
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}
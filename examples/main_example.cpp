#include "bqf_ec.hpp"
#include "bqf_cf.hpp"

#include <chrono>

std::chrono::steady_clock::time_point begin;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ./bqf <command>" << std::endl;
        std::cout << "Commands:" << std::endl;
        std::cout << "./bqf build [-q <quotient size=8>] [-c <count size=5>] [-k <k=32>] [-z <z=5>] -i <counted_smers> -o <BQF_file>" << std::endl;
        std::cout << "./bqf query -b <bqf_file> -i <reads_to_query> -o <results>" << std::endl;
        std::cout << "./bqf filter [-q <quotient size=8>] [-k <k=32>] -i <fasta_file> -o <outfile>" << std::endl;
        std::cout << "./bqf help" << std::endl;
        return EXIT_FAILURE;
    }

    std::string command = argv[1];

    std::string input_file;
    std::string output_file;

    if (command == "build") {
        int q = 8;
        int c = 5;
        int k = 32;
        int z = 11;

        for (int i = 2; i < argc; i++) {
            if (std::string(argv[i]) == "-q") {
                if (i + 1 < argc) {
                    q = std::stoi(argv[++i]);
                } else {
                    std::cerr << "The -q option requires a value." << std::endl;
                    return EXIT_FAILURE;
                }
            } else if (std::string(argv[i]) == "-c") {
                if (i + 1 < argc) {
                    c = std::stoi(argv[++i]);
                } else {
                    std::cerr << "The -c option requires a value." << std::endl;
                    return EXIT_FAILURE;
                }
            } else if (std::string(argv[i]) == "-k") {
                if (i + 1 < argc) {
                    k = std::stoi(argv[++i]);
                } else {
                    std::cerr << "The -k option requires a value." << std::endl;
                    return EXIT_FAILURE;
                }
            } else if (std::string(argv[i]) == "-z") {
                if (i + 1 < argc) {
                    z = std::stoi(argv[++i]);
                } else {
                    std::cerr << "The -z option requires a value." << std::endl;
                    return EXIT_FAILURE;
                }
            } else if (std::string(argv[i]) == "-i") {
                if (i + 1 < argc) {
                    input_file = argv[++i];
                } else {
                    std::cerr << "The -i option requires an input file name." << std::endl;
                    return EXIT_FAILURE;
				}
        	} else if (std::string(argv[i]) == "-o") {
                if (i + 1 < argc) {
                    output_file = argv[++i];
                } else {
                    std::cerr << "The -o option requires an output file name." << std::endl;
                    return EXIT_FAILURE;
                }
            } else {
                std::cerr << "Invalid argument : " << argv[i] <<std::endl;
                return EXIT_FAILURE;
            }
        }

        if (q <= 7 || c <= 0) {
            std::cerr << "Values of q, and c must be greater than 7, and 0." << std::endl;
            return EXIT_FAILURE;
        }

        begin = std::chrono::steady_clock::now();

		Bqf_ec bqf = Bqf_ec(q, c, k, z, false);
		bqf.insert(input_file);
		bqf.save_on_disk(output_file);

        std::cout << "BQF constructed successfully and saved to " << output_file << std::endl;
        std::cout << "Construction time = " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - begin).count() << "ms" << std::endl;
    } 
	
	else if (command == "query") {

        std::string input_reads_file_to_query;
        for (int i = 2; i < argc; i++) {
            if (std::string(argv[i]) == "-b") {
                if (i + 1 < argc) {
                    input_file = argv[++i];
                } else {
                    std::cerr << "The -b option requires a value." << std::endl;
                    return EXIT_FAILURE;
                }
            } else if (std::string(argv[i]) == "-i") {
                if (i + 1 < argc) {
                    input_reads_file_to_query = argv[++i];
                } else {
                    std::cerr << "The -i option requires an input file name." << std::endl;
                    return EXIT_FAILURE;
                } 
            } else if (std::string(argv[i]) == "-o") {
                if (i + 1 < argc) {
                    output_file = argv[++i];
                } else {
                    std::cerr << "The -o option requires an output file name." << std::endl;
                    return EXIT_FAILURE;
                }
            } else {
                std::cerr << "Invalid argument : " << argv[i] <<std::endl;
                return EXIT_FAILURE;
            }
        }
        
        if (input_file.empty() || input_reads_file_to_query.empty()) {
            std::cerr << "Input file names are missing." << std::endl;
            std::cerr << "- input_file.empty()                = " << input_file.empty() << std::endl;
            std::cerr << "- input_file                        = " << input_file         << std::endl;
            std::cerr << "- input_reads_file_to_query.empty() = " << input_reads_file_to_query.empty() << std::endl;
            std::cerr << "- input_reads_file_to_query         = " << input_reads_file_to_query         << std::endl;
            std::cerr << "Error location : file = [" << __FILE__ << "] @ line " << __LINE__ << std::endl;
            return EXIT_FAILURE;
        }

        begin = std::chrono::steady_clock::now();
		Bqf_ec bqf = Bqf_ec::load_from_disk(input_file);

		try {
			std::ifstream infile(input_reads_file_to_query);
            std::ofstream outfile(output_file);

			if (!infile) {
                std::cerr << "Exception happened..." << std::endl;
                std::cerr << "Error location : file = [" << __FILE__ << "] @ line " << __LINE__ << std::endl;
				throw std::runtime_error("File not found: " + input_reads_file_to_query);
			}
            if (!outfile) {
                std::cerr << "Exception happened..." << std::endl;
                std::cerr << "Error location : file = [" << __FILE__ << "] @ line " << __LINE__ << std::endl;
				throw std::runtime_error("File can not be created: " + output_file);
			}

            bqf.query(infile, outfile);
			
		} catch (const std::exception &e) {
			std::cerr << "Error: " << e.what() << std::endl;
		}


        std::cout << "Queries executed successfully." << std::endl;
        std::cout << "Load + queries time = " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - begin).count() << "ms" << std::endl;
    
    } else if (command == "filter") {
        int q = 8;
        int k = 32;
        std::vector<std::string> filenames;

        for (int i = 2; i < argc; i++) {
            if (std::string(argv[i]) == "-q") {
                if (i + 1 < argc) {
                    q = std::stoi(argv[++i]);
                } else {
                    std::cerr << "The -q option requires a value." << std::endl;
                    return EXIT_FAILURE;
                }
            } else if (std::string(argv[i]) == "-k") {
                if (i + 1 < argc) {
                    k = std::stoi(argv[++i]);
                } else {
                    std::cerr << "The -k option requires a value." << std::endl;
                    return EXIT_FAILURE;
                }
            } else if (std::string(argv[i]) == "-i") {
                if (i + 1 < argc) {
                    input_file = argv[++i];
                    filenames.push_back(input_file);
                } else {
                    std::cerr << "The -i option requires an input file name." << std::endl;
                    return EXIT_FAILURE;
				}
        	} else if (std::string(argv[i]) == "-o") {
                if (i + 1 < argc) {
                    output_file = argv[++i];
                } else {
                    std::cerr << "The -o option requires an output file name." << std::endl;
                    return EXIT_FAILURE;
                }
            } else {
                std::cerr << "Invalid argument : " << argv[i] <<std::endl;
                return EXIT_FAILURE;
            }
        }

        if (q < 7) {
            std::cout << "Value of q must be greater than 7" << std::endl;
            return EXIT_FAILURE;
        }
        if (input_file.empty()) {
            std::cerr << "Input file name is missing" << std::endl;
            return EXIT_FAILURE;
        }
        if (output_file.empty()) {
            std::cerr << "Output file name is missing" << std::endl;
            return EXIT_FAILURE;
        }

        try {
            begin = std::chrono::steady_clock::now();

            Bqf_cf bqf = Bqf_cf(q, k);
            bqf.filter_fastx_file(filenames, output_file);

            std::cout << "File successfully filtered. " << k << "-mers present more than once are stored in " << output_file << std::endl;
            std::cout << "Total time = " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - begin).count() << "ms" << std::endl;
        }
        catch (const std::exception &e) {
            std::cerr << "Error : " << e.what() << std::endl;
        }
        

    } else if (command == "help" || command == "h") {
        std::cout << "Usage: ./bqf <command>" << std::endl;
        std::cout << "Commands:" << std::endl;
        std::cout << "./bqf build -q <quotient size> [-c <count size=5>] [-k <k=32>] [-z <z=5>] -i <counted_smers> -o <BQF_file>" << std::endl;
        std::cout << "./bqf query -b <bqf_file> -i <reads_to_query> -o <results_file>" << std::endl;
        std::cout << "./bqf filter [-q <quotient size=8>] [-k <k=32>] -i <fasta_file> -o <outfile>" << std::endl;
        std::cout << "./bqf help" << std::endl;

        std::cout << "-q is quotient size, it sets the filter size (there will be 2^q slots) so 2^(q-1) < nb_unique_elements < 2^q is higly recommanded" << std::endl;
        std::cout << "-c is the number of bits reserved for counters of each element. (2^c)-1 will be the maximum value" << std::endl;
        std::cout << "-k is the kmer size. The result of the query of a sequence S will be determined by the queries of all the kmers of S" << std::endl;
        std::cout << "-z is fimpera parameter. kmers are queried through the query of all their smers. s = k-z and smers are effectively inserted in the filter" << std::endl;
        std::cout << "-i is input_file, can be counted smers for \"build\" tool, sequences to query for \"query\" tool or a fasta/q file for \"filter\" tool" << std::endl;
        std::cout << "-o is the file on which the BQF is saved in binary form after building (weights around 2^q*(3+c+r) bits, r being 2s-q) for \"build\" tool, sequences results file for \"query\" tool or a list of kmer for \"filter\" tool" << std::endl;
        std::cout << "-b is the file from which the BQF is loaded" << std::endl;
    } else {
        std::cerr << "Invalid command or incorrect number of arguments." << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
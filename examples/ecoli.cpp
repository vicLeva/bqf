#include <iostream>

#include "bqf_ec.hpp" 
#include "additional_methods.hpp"

std::string cwd = DATA_DIR;

using namespace std;



int main (int argc, char * argv[]) {
        cout << 
        "\nExample of a user indexing genomic data from ecoli genome (examples/data/)\n"
        "Goal : query sequences of length 32\n"
        "=> index 28-mers instead (k=32, s=28), hash will be 2*28=56 bits long\n"
        "kmc -k28 gives a 28-mers counting of the genome, showing 182.332 distinct smers,\n"
        "at least 182.332 unique slots are needed\n"
        "=> 2^17=130.000 < 182.332 < 2^18=260.000 so we choose q=18 and r=56-18\n\n"
        "with 5 bits for counter, the quotient filter will weight 2^18*(3+5+38)bits = 1.5MB\n\n";

        Bqf_ec ecoli_cqf(18, 56-18, 5, false);
        ecoli_cqf.insert(cwd + "data/ecoli_count28.txt");

        cout << ecoli_cqf.query("CGTATGAAAAAACGTATCCCCACTCTCC", 28).minimum << endl; //28
        cout << ecoli_cqf.query("CTCGTTGATACCGGGCGCCCTGCTTGCC", 28).minimum << endl; //1
        cout << ecoli_cqf.query("AATCTGAAAGATCACAACGAGCAGGTCA", 28).minimum << endl; //25
        cout << ecoli_cqf.query("TTTTTTTTTATAAAAACGAGCAGAACAA", 28).minimum << endl; //0

        ecoli_cqf.remove("CGTATGAAAAAACGTATCCCCACTCTCC", 6);

        cout << ecoli_cqf.query("CGTATGAAAAAACGTATCCCCACTCTCC", 28).minimum << " removed 6 from 28" << endl;

        cout << ecoli_cqf.query("ACAACGTTTGCTCGATGATCGCCTGCTCATCG", 28).minimum << endl; //30
        cout << ecoli_cqf.query("ACCTTTGTTCAGGCGAGTCAGGGCGTCACGGA", 28).minimum << endl; //24 but 25 because of FalsePositive

        result_query bigQuery = ecoli_cqf.query("ACCTTTGTTCAGGCGAGTCAGGGCGTCACGGAG", 32);
        cout << bigQuery << endl;

        

        ecoli_cqf.save_on_disk("/udd/vlevallo/documents/ecoli_cqf");
        std::cout << ecoli_cqf.block2string(0) << "\n";

        Bqf_ec ecoli_cqf2 = Bqf_ec::load_from_disk("/udd/vlevallo/documents/ecoli_cqf");

        uint64_t val = (3ULL << 18) + 8; 
        ecoli_cqf2.insert(val, 1);

        cout << ecoli_cqf2.query(val) << endl;

        std::cout << ecoli_cqf2.block2string(0) << "\n";
        cout << ecoli_cqf2.query("CGTATGAAAAAACGTATCCCCACTCTCC", 28).minimum << endl; //28
        cout << ecoli_cqf2.query("CTCGTTGATACCGGGCGCCCTGCTTGCC", 28).minimum << endl; //1
        cout << ecoli_cqf2.query("AATCTGAAAGATCACAACGAGCAGGTCA", 28).minimum << endl; //25
        cout << ecoli_cqf2.query("TTTTTTTTTATAAAAACGAGCAGAACAA", 28).minimum << endl; //0

	
	return 0;
}
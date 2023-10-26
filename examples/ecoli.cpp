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

        Bqf_ec ecoli_cqf(18, 5, 32, 4, false);
        ecoli_cqf.insert(cwd + "data/ecoli_count28.txt");

        
        ecoli_cqf.save_on_disk("/udd/vlevallo/documents/ecoli_bqf");
        std::cout << ecoli_cqf.block2string(0) << "\n";

        Bqf_ec ecoli_cqf2 = Bqf_ec::load_from_disk("/udd/vlevallo/documents/ecoli_bqf");

        std::cout << ecoli_cqf2.block2string(2244) << "\n";
        cout << ecoli_cqf2.query("ACAACGTTTGCTCGATGATCGCCTGCTCATCG").minimum << endl; //30
        cout << ecoli_cqf2.query("ACCTTTGTTCAGGCGAGTCAGGGCGTCACGGA").minimum << endl; //25
        cout << ecoli_cqf2.query("ACCTTTGTTCAGGCGAGTCAGGGCGTCACGGAG").minimum << endl; //0

	
	return 0;
}
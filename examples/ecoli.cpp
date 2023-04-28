#include <iostream>

#include "bcqf_ec.hpp" 
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

        Bcqf_ec ecoli_cqf(18, 56-18, 5, false);
        ecoli_cqf.insert(cwd + "data/ecoli_count.txt");

        cout << ecoli_cqf.query("CGTATGAAAAAACGTATCCCCACTCTCC", 28) << endl; //28
        cout << ecoli_cqf.query("CTCGTTGATACCGGGCGCCCTGCTTGCC", 28) << endl; //1
        cout << ecoli_cqf.query("AATCTGAAAGATCACAACGAGCAGGTCA", 28) << endl; //25
        cout << ecoli_cqf.query("TTTTTTTTTATAAAAACGAGCAGAACAA", 28) << endl; //0

        ecoli_cqf.remove("CGTATGAAAAAACGTATCCCCACTCTCC", 6);

        cout << ecoli_cqf.query("CGTATGAAAAAACGTATCCCCACTCTCC", 28) << " removed 6 from 28" << endl;

        cout << ecoli_cqf.query("ACAACGTTTGCTCGATGATCGCCTGCTCATCG", 32) << endl; //30
        cout << ecoli_cqf.query("ACCTTTGTTCAGGCGAGTCAGGGCGTCACGGA", 32) << endl; //24 but 25 because of FalsePositive

	
	return 0;
}
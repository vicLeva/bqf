#include <iostream>


using namespace std;

int main (int argc, char * argv[]) {
	cout <<
        "I'm a user and my name is taro_oceans.\n"
        "I have 4.7 billions 32-mers to index (2.64Bil. are distinct, 1.98Bil. appear once)\n"
        "With z=4 (z=k-s), 4.97Bil. 28-mers will effectively be inserted into the BCQF\n"
        "inserting 28-mers means a 56bits hash function to make sure we have access to revertible hash function (resize)\n"
        "There are 2.62Bil. distinct 28-mers so 32 bits among the 56 of the hash are needed : q=32\n"
        "q=32 means 2^32 ~= 4.3Bil. thus we'll have 2.62/4.3 = 60% load\n"
        "q=32 and r = 56-32 = 24, with 5 bits for each counter, we expect\n"
        "2^32 * (3+24+5) bits = 17Gigabytes index\n";


	
	return 0;
}
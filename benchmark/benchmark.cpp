/*
	Generates n random integers and feed the cqf with them.
*/
#include "rsqf.hpp"

#include <fstream>
#include <random>
#include <ctime>
#include <getopt.h>
#include <chrono>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>

using namespace std;

std::string cwd = BENCHMARK_DIR;


void writeQueryTimes(vector<int> & vec){
	int len = vec.size();

	ofstream myfile;
  	myfile.open(cwd + "statsFiles/tmpQueryTime");

	vector<Rsqf> qfs;
	for (int i=0; i<len; i++){ 
		myfile << "q" << vec[i];
		qfs.push_back(Rsqf(vec[i], 64-vec[i], false)); 
	}
	myfile << "\n";

	// uint64 generators
	default_random_engine generator;
	generator.seed(time(NULL));
	uniform_int_distribution<uint64_t> distribution;
	uint64_t val;

	std::chrono::_V2::system_clock::time_point t1, t2;
	

	for (int i=0; i<len; i++){ //for each qf
		unordered_set<uint64_t> inserted;
		uint64_t max_load = ((1ULL << vec[i]) - 1);
		uint64_t load = 0;

		for (float perc=1; perc<96; perc++){ //for each load factor (1-95%)
			while (load < round(max_load * (perc/100))){ //fill 1% more
				load++;
				val = distribution(generator);  
				inserted.insert(val);    
				qfs[i].insert(val);
			}

			//100k positive queries
			t1 = std::chrono::high_resolution_clock::now();
			int cpt = 100000;
			while (true) { 
				for (auto itr = inserted.begin(); itr != inserted.end(); ++itr) {
					qfs[i].query(*itr);
					if (--cpt == 0){ goto out_while; } //break out of while loop 
				}
			}
			out_while:
			t2 = std::chrono::high_resolution_clock::now();
			myfile << to_string( std::chrono::duration<double, std::milli>( t2 - t1 ).count()) << "\n";
			

			//100k negative queries
			t1 = std::chrono::high_resolution_clock::now();
			for (int j=0; j<100000; j++){ 
				val = distribution(generator); //assumed not inserted yet
				qfs[i].query(val);
			}
			t2 = std::chrono::high_resolution_clock::now();
			myfile << to_string( std::chrono::duration<double, std::milli>( t2 - t1 ).count()) << "\n";
		}

	}

	myfile.close();
}

void writeInsertTimes(int q){
	//q >= 11 for having 100insert less than 0.5% load addition
	uint64_t max_load = ((1ULL << q) - 1);
	uint64_t half_percent_load = max_load/200;

	ofstream myfile;
  	myfile.open(cwd + "statsFiles/tmpInsertTime");
	myfile << "q" << q << ";unit" << half_percent_load << "\n";

	vector<Rsqf> qfs;
	for (int i=0; i<10; i++){ 
		qfs.push_back(Rsqf(q, 64-q, false)); 
	}

	// uint64 generators
	default_random_engine generator;
	generator.seed(time(NULL));
	uniform_int_distribution<uint64_t> distribution;
	uint64_t val;

	std::chrono::_V2::system_clock::time_point t1, t2;
	
	
	uint64_t load = 0;
	for (float perc=1; perc<96; perc++){ //for each load factor (1-95%)
		while (load < round(max_load * (perc/100))){ //fill 1% more
			load++;
			for (size_t i=0; i<qfs.size(); i++){ //for each qf
				val = distribution(generator);  
				qfs[i].insert(val);
			}
		}


		load = load + half_percent_load;
		t1 = std::chrono::high_resolution_clock::now();
		for (uint64_t j=0; j<half_percent_load; j++){ //some inserts in each qf
			for (size_t i=0; i<qfs.size(); i++){ //for each qf
				val = distribution(generator);  
				qfs[i].insert(val);
			}
		}
		t2 = std::chrono::high_resolution_clock::now();
		myfile << to_string( std::chrono::duration<double, std::milli>( t2 - t1 ).count()) << "\n";
	}
	

	myfile.close();
}

void printMaxMemUsage(int mem){
	cout << cwd + "statsFiles/tmpMaxMem" << endl;

	std::string cmd = "/usr/bin/time -v ";
	cmd = cmd + cwd + "/../build/bin/random_uints -s " + to_string(mem) + " -n " + to_string(500);
	cmd = cmd + " 2> " + cwd + "statsFiles/tmpMaxMem";
	system(cmd.c_str());

	ifstream myfile;
  	myfile.open(cwd + "statsFiles/tmpMaxMem");
	int c=0;
	for( std::string line; getline( myfile, line ); ) {
		if(++c == 10) cout << "Max memory usage : " << line.substr(37) << " KB while building " << mem << "MB QF" << endl;
	}
}

int main (int argc, char * argv[]) {
	
	vector<int> qsizes = {7, 17, 20};
	writeQueryTimes(qsizes);
	
	
	//q >= 11 for having 100insert less than 0.5% load addition
	writeInsertTimes(20);

	//mem = 4  ==>  q = 19
	printMaxMemUsage(4);
	
	return 0;
}

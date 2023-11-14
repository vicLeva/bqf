# cqf_implementation

Implementation of a Counting & Backpack Quotient Filter

## Compilation of the project
  
From the project root
```bash
cmake -B build
cd build && make 
```

## Tool

From build/bin/

### Usage  
```bash
./bqf <command> [parameters]
```

### Commands:  
```bash
./bqf build -q <quotient size> [-c <count size=5>] [-k <k=32>] [-z <z=5>] -i <counted_smers> -o <BQF_file>
./bqf query -b <bqf_file> -i <reads_to_query>
./bqf help
```

### Parameters  
  + `-q` is quotient size, it sets the filter size (there will be 2^q slots) so 2^(q-1) < nb_unique_elements < 2^q is needed
  + `-c` is the number of bits reserved for counters of each element. 2^c will be the maximum value
  + `-k` is the kmer size. The result of the query of a sequence S will be the minimum of the queries of all the kmers of S
  + `-z` is [fimpera](https://academic.oup.com/bioinformatics/article/39/5/btad305/7169157) parameter. kmers are queried through the query of all their smers. s = k-z and smers are effectively inserted in the filter
  + `-i` is input_file, can be counted smers for `build` command (usually from [KMC](https://github.com/refresh-bio/KMC) or equivalent) or sequences to query for `query` command (1 sequence / line)
  + `-o` is the file on which the BQF is saved in binary form after building (weights around 2^q*(3+c+r) bits, r being 2s-q)
  + `-b` is the file from which the BQF is loaded


### Experiments details

Protocol available in the [Wiki-protocol](https://github.com/vicLeva/bqf/wiki/Experiments-details-and-protocol) page of the repository

### Examples

  1. + `./bqf build -q 18 -z 4 -i data/ecoli_count28.txt -o /tmp/ecoli_bqf`
       - build a 2^18 slots filter with (32-4 = 28)-mers aiming to query 32-mers later. 5 bits for counters, max value =2^5=64  
     + `./bqf query -b /tmp/ecoli_bqf -i data/queries.fasta`
       - load bqf then query each line (=sequence) of the file given with `-i`

  2. + `./bqf build -q 31 -c 5 -k 32 -z 10 -i /scratch/vlevallois/data/AHX_ACXIOSF_6_1_22_all.txt -o /scratch/vlevallois/bqf_tmp`
     + `./bqf query -b /scratch/vlevallois/bqf_tmp -i ~/data/queries.fasta`


## Documentation

The documentation can be generated using doxygen with the following command in the root project
```bash
doxygen Doxyfile
```
Then you can find a html file (`index.html`) in the so-created html directory.


## Unitary tests

From build directory
```bash
ctest
```

## Benchmark

in build/
```bash
make && bin/benchmark && python ../benchmark/plot.py
```


## Flamegraph

From root directory  
```bash
flamegraph/flamegraph.sh [quicktests|random_uints|unit_tests]
```

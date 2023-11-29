# BQF - Backpack Quotient Filter

Implementation of a Counting & Backpack Quotient Filter

## About

The [docs](docs) directory contains slides used to present the BQF during [SeqBim2023](https://seqbim.cnrs.fr/seqbim-2023/).

### Overview

A variant of the original [counting quotient filter](https://github.com/splatlab/cqf). The quotient filter is a hash table-like data structure where part of the information inserted is stored implicitly within the address in the table where it is written.
At the price of a slight non-null false positive rate ($10^{-11}$ in our experiments), The BQF is more space-efficient than the CQF thanks to the way it handles the abundance of indexed elements. In the BQF, each slot stores a fingerprint plus a counter using itself *c* bits.
Having the counter attached to every slot, every fingerprint holds its own count. 

Using [fimpera](https://github.com/lrobidou/fimpera), the space required by the $c$ bits per stored element does not impact the final structure size. 
With fimpera, instead of storing *k-mers*, the structure indexes *s-mers* with $s < k$. At query time a *kmer* is considered as present if all its *smers* are present. This induces the tiny non-null false positive rate (as long as *smers* are big enough) but this frees $2 \times (k-s)$ bits per element, used for storing the count of elements.

### Comparison with CQF

If a kmer is present twice or more, then using only one slot becomes more space-efficient than CQF encoding scheme. The strategy results to a double benefit for the BQF, **less space is needed for every slots**, and **less slots are used** (for abundances > 2), meaning we postpone the doubling up of the structure when dynamically inserting elements. All of this is at the cost of a negligible false positive rate if parameterized correctly.   

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
  
(binaries in build/bin/)  
  
1. + `./bqf build -q 18 -z 4 -i examples/data/ecoli_count28.txt -o /tmp/ecoli_bqf`
     - build a 2^18 slots filter with (32-4 = 28)-mers aiming to query 32-mers later. 5 bits for counters, max value =2^5=64  
   + `./bqf query -b /tmp/ecoli_bqf -i examples/data/queries.fasta`
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

## BQF - Backpack Quotient Filter

## Overview

The Backpack Quotient Filter (BQF) is an indexing data structure with abundance. Although the data can be anything, it's been thought to index genomic datasets. 
The BQF is a dynamic structure, with a correct hash function it can add, delete and enumerate elements. Thus the structure can resize itself when needed. The main features are **indexing** (building the BQF over a dataset *D*) and **querying** (searching for a sequence in *D*)

The BQF is able to index metagenomics datasets (low redundancy, high complexity datasets) with an average of 25 bits per element. This value tends to lower as the datasets grow. Compared to the main variant, the [Counting Quotient Filter](https://github.com/splatlab/cqf) (CQF), the BQF is 4 to 5 times smaller according to our [experiments](#Experiments-results).

It relies on a hash-table-like structure called Quotient Filter. Part of the information inserted is stored implicitly within the address in the table where it is written

## ToC

 + [Installation](#Installation)
 + [Tool usage](#Tool-usage)
 + [Examples](#Examples)
 + [API Documentation](#Documentation)
 + [Unitary tests](#Unitary-tests)
 + [Slides presentation](#Slides-presentation)
 + [Paper experiments results](#Experiments-results)

## Installation

```bash
git clone git@github.com:vicLeva/bqf.git
cd bqf
mkdir build && cd build
cmake ..
make 
```

## Binary usage 

From `bqf/build/`

```
./bin/bqf [TOOL] [PARAMETERS]

[TOOL] : 
    build           Build a BQF from a counted s-mers file
    query           Query sequence(s) in a BQF
    help            Display commands and parameters


[PARAMETERS] : build
    -q  (default = 8)   Quotient size, defines BQF's size (2^q slots). For n distinct s-mers to index, it is advised to initialize q to ceil(log2(n)). In other terms : log2(n) < q is recommanded with smallest q possible

    -c  (default = 5)   Counters size, number of bits used to encode abundances. Max value is (2^c)-1. Incrementing c from x to x+1 will increase BQF's size of 2^q bits

    -k  (default = 32)  k-mer size. When querying a sequence S in BQF, all substrings of size k of S are queried.

    -z  (default = 11)  Fimpera parameter. BQF inserts s-mers, of size s=k-z. Abundances of k-mers are virtualized through abundances of s-mers. Increasing z by 1 lowers BQF size by 2(2^q) bits but increases false positive rate. Negligible at first, it increases exponentially around z=15 (with k=32).

    -i  [mandatory]     Input file path. Must be counted s-mers. Output of KMC (https://github.com/refresh-bio/KMC) execution is recommanded.

    -o  [mandatory]     Output file path. Binary BQF on disk.


[PARAMETERS] : query
    -b  [mandatory]     BQF file path.

    -i  [mandatory]     Input file path. Sequences to query in the BQF, 1 per line or .fasta format.

    -o  [mandatory]     Output file path. Results of queries, 1 per line. Results are formated this way (min:X, max:X, average:X, presence ratio:X). min, max and average are k-mers abundances statistics of the queried sequence S. presence ratio is the ratio of present k-mers over all k-mers of S.
```

## Examples 
  
1. + `./bqf build -q 18 -z 4 -i examples/data/ecoli_count28.txt -o /tmp/ecoli_bqf`
     - build a 2^18 slots filter with (32-4 = 28)-mers aiming to query 32-mers later. 5 bits for counters, max value =2^5=64  
   + `./bqf query -b /tmp/ecoli_bqf -i examples/data/queries.fasta -o ./results`
     - load bqf then query each line (=sequence) of the file given with `-i`

2. + `./bqf build -q 31 -c 5 -k 32 -z 10 -i /scratch/vlevallois/data/AHX_ACXIOSF_6_1_22_all.txt -o /scratch/vlevallois/bqf_tmp`
   + `./bqf query -b /scratch/vlevallois/bqf_tmp -i ~/data/queries.fasta -o ./results`


## Documentation

The documentation can be generated using doxygen with the following command in the root project

```bash
doxygen Doxyfile
```

Then you can find a html file (`index.html`) in the so-created html directory.

## Unitary tests

From `build/` directory

```bash
ctest
```

## Slides presentation

Available [here](https://vicleva.github.io/) (talks section) or through this [download link](https://vicleva.github.io/assets/slides/presentation_seqbim_2023.pdf) (from [SeqBim2023](https://seqbim.cnrs.fr/seqbim-2023/) conference).


## Experiments results

All datasets used and experiments done are detailed in the [wiki](https://github.com/vicLeva/bqf/wiki) or directly by this [link](https://github.com/vicLeva/bqf/wiki/Experiments-details-and-protocol-for-BQF-paper-results).
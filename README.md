## BQF - Backpack Quotient Filter

**Preprint** : https://www.biorxiv.org/content/10.1101/2024.02.15.580441v1

## Overview

The Backpack Quotient Filter (BQF) is an indexing data structure with abundance. Although the data can be anything, it's been thought to index genomic datasets. 
The BQF is a dynamic structure, with a correct hash function it can add, delete and enumerate elements. Thus the structure can resize itself when needed. The main features are **indexing** (building the BQF over a dataset *D*) and **querying** (searching for a sequence in *D*)

The BQF is able to index metagenomics datasets (low redundancy, high complexity datasets) with an average of 25 bits per element. This value tends to lower as the datasets grow. Compared to the main variant, the [Counting Quotient Filter](https://github.com/splatlab/cqf) (CQF), the BQF is 4 to 5 times smaller according to our [experiments](https://github.com/vicLeva/bqf/wiki/Experiments-details-and-protocol-for-BQF-paper-results).

It relies on a hash-table-like structure called Quotient Filter. Part of the information inserted is stored implicitly within the address in the table where it is written. You are going to read about *k*-mers and *s*-mers, both are words of size *k* or *s* with *k* $\geq$ *s*. BQF inserts and query *s*-mers but virtualizes the presence of *k*-mers at query time. In other words, a query sequence is broken down into *k*-mers, and each *k*-mer is virtually queried through all of its *s*-mers.  

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

(KMC step example with Dataset XX, count 32-mers)

``` bash 
kmc -k [32] -m24 -ci2 -v XX.fastq XX.res tmp_working_dir
kmc_tools dump XX.res XX_counted
rm XX.res.kmc_*
```

From `bqf/build/`
  
1. + `./bin/bqf build -q 18 -z 4 -i ../examples/data/ecoli_28_counted -o /tmp/ecoli_bqf`
     - build a 2^18 slots filter with (32-4 = 28)-mers aiming to query 32-mers later. 5 bits for counters, max value =2^5=64  
   + `./bin/bqf query -b /tmp/ecoli_bqf -i ../examples/data/queries.fasta -o ./results.out`
     - load bqf then query each sequence of the file given with `-i`

   + (`rm /tmp/ecoli_bqf`)


Real scale data example with [this dataset](ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR172/ERR1726642/AHX_ACXIOSF_6_1_C2FGHACXX.IND4_clean.fastq.gz) from Tara Oceans metagenomic project (7.7GB). Assuming 19-mers have been counted with KMC.
2. + `./bin/bqf build -q 29 -c 5 -k 32 -z 13 -i /path/to/6_1_19_counted -o /path/to/6_1_bqf`
   + `./bin/bqf query -b /path/to/6_1_bqf -i ../examples/data/queries.fasta -o ./results.out`


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

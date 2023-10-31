## Retrieve data
data @genouest : 
/groups/genscale/NGSdatasets/metagenomics/TARA/tara_set3/QQSS/11SUR1QQSS11/AHX_ACXIOSF_6_1_C2FGHACXX.IND4_clean.fastq.gz

`srun --mem 64G -c 4 -p genscale --pty bash`  
`source /local/env/envcmake-3.21.3.sh`  
`source /local/env/envgcc-11.2.0.sh`  

## File kmer count
`source /local/env/envkmc-2.3.0.sh`  
  + count for indexing  
`kmc -k23 -m24 -ci1 -v AHX_ACXIOSF_6_1_C2FGHACXX.IND4_clean.fastq /scratch/vlevallois/data/6_1_23.res /scratch/vlevallois/kmc_working/`  
`kmc_tools dump 6_1_23.res AHX_ACXIOSF_6_1_23_all.txt`  
  - 1.607.526.303 unique 23-mers
  - AHX_ACXIOSF_6_1_C2FGHACXX.IND4_clean.fastq *7.7GB*
  - AHX_ACXIOSF_6_1_23_all.txt *39GB*

  + count for verifying  
`kmc -k32 -m24 -ci1 -v AHX_ACXIOSF_6_1_C2FGHACXX.IND4_clean.fastq /scratch/vlevallois/data/6_1_32.res /scratch/vlevallois/kmc_working/`  
`kmc_tools dump 6_1_32.res AHX_ACXIOSF_6_1_32_all.txt`  

## BQF
https://github.com/frankandreace/cqf_implementation/tree/dev  
`cmake -B build`  
`cd build && make`  
`./bin/bqf build -q 31 -c 5 -k 32 -z 9 -i /scratch/vlevallois/data/AHX_ACXIOSF_6_1_23_all.txt -o /scratch/vlevallois/bqf_tmp`
  - 1385s (~23min)
`./bin/bqf query -b /scratch/vlevallois/bqf_tmp -i ~/data/queries.fasta`
  - 190s for 772.322.833 positive smer queries => 4.064.857 smer/s
  - 127s for 654.974.399 negative smer queries => 5.124.313 smer/s

## CQF 
https://github.com/splatlab/cqf/tree/68939f587b0f33f8c45be4f0e2a18e12bffa9fa6  
`. /local/env/envconda.sh`  
`conda create -p ~/my_env gcc cmake openssl`  
`conda create -p ~/my_env gcc cmake openssl`  
`conda activate ~/my_env`  
(in case openssl bug include `cp -r ~/my_env/include/openssl /home/genouest/genscale/vlevallois/my_env/bin/../lib/gcc/x86_64-conda-linux-gnu/12.2.0/include/`)  
  
`make test`  
`./test 31 33` (src/test.c)  
  - build (inserts) : 1186s (~20min)
  - 266s for 682.454.044 positive kmer queries => 2.565.616 kmer/s
  - 190s for 564.974.399 negative kmer queries => 2.973.549 kmer/s
  - 1.951.662.794 slots used (90%)

## cBF
`cd build && make && ./app/fimpera_index -k 32 -z 4 --canonical /scratch/vlevallois/data/AHX_ACXIOSF_6_1_32_all.txt out.out 65065438397`  
(app/src/index.cpp)  

## Retrieve data
data @genouest : 
/groups/genscale/NGSdatasets/metagenomics/TARA/tara_set3/QQSS/11SUR1QQSS11/AHX_ACXIOSF_6_1_C2FGHACXX.IND4_clean.fastq.gz

`srun --mem 64G -c 4 -p genscale --pty bash`  
`source /local/env/envcmake-3.21.3.sh`  
`source /local/env/envgcc-11.2.0.sh`  

## File kmer count
`source /local/env/envkmc-2.3.0.sh`  
count for indexing  
`kmc -k28 -m24 -ci1 -v AHX_ACXIOSF_6_1_C2FGHACXX.IND4_clean.fastq /scratch/vlevallois/data/6_1_28.res /scratch/vlevallois/kmc_working/`  
`kmc_tools dump 6_1_28.res AHX_ACXIOSF_6_1_28_all.txt`  
count for verifying  
`kmc -k32 -m24 -ci1 -v AHX_ACXIOSF_6_1_C2FGHACXX.IND4_clean.fastq /scratch/vlevallois/data/6_1_32.res /scratch/vlevallois/kmc_working/`  
`kmc_tools dump 6_1_32.res AHX_ACXIOSF_6_1_32_all.txt`  

## BQF
https://github.com/frankandreace/cqf_implementation/tree/dev  
`cmake -B build`  
`cd build`  
`./bin/quicktests` (test_8_GB() in test/test.cpp)  

## CQF 
https://github.com/splatlab/cqf/tree/68939f587b0f33f8c45be4f0e2a18e12bffa9fa6  
`. /local/env/envconda.sh`  
`conda create -p ~/my_env gcc cmake openssl`  
`conda create -p ~/my_env gcc cmake openssl`  
`conda activate ~/my_env`  
(in case openssl bug include `cp -r ~/my_env/include/openssl /home/genouest/genscale/vlevallois/my_env/bin/../lib/gcc/x86_64-conda-linux-gnu/12.2.0/include/`)  
  
`make test`  
`./test 31 33` (src/test.c)  

## cBF
`cd build && make && ./app/fimpera_index -k 32 -z 4 --canonical /scratch/vlevallois/data/AHX_ACXIOSF_6_1_32_all.txt out.out 65065438397`  
(app/src/index.cpp)  

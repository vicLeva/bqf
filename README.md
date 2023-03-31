# cqf_implementation

Implementation of a RSQF

## Compilation of the project
  
else

From the project root
```bash
cmake -B build
cmake --build build 
```
  
or   
  
From the project root
```bash
cmake -B build
cd build && make 
```

## Examples

**Binaries are in build/bin/**  

From project root  

```bash
./build/bin/random_uints -h
./build/bin/random_uints -n 500 -s 6 (insert 500 numbers in a 6MB RSQF)
```

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

# Benchmark

in build/
```bash
make && bin/benchmark && python ../benchmark/plot.py
```



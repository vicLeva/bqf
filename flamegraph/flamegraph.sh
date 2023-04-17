#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "Need 1 binary as arg (cqf_test|random_uints|unit_tests)"
    exit
fi


cd build
make
perf record -o ../flamegraph/perf.data --call-graph dwarf ./bin/"$1"
cd ../flamegraph
perf script | ./stackcollapse-perf.pl |./flamegraph.pl > perf.svg
firefox perf.svg

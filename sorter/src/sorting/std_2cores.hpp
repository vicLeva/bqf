#pragma once
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <omp.h>

#include "../tools/tools.hpp"

extern void my_sort(std::vector<uint64_t>& test);
extern void p_sort(std::vector<uint64_t>& test);

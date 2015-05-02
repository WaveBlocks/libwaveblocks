#ifndef PROFILING_HELPER_HPP
#define PROFILING_HELPER_HPP

#include <cstddef>


void incr_flag();
void decr_flag();

std::size_t total_sample_count();
std::size_t positive_sample_count();

void start_profiler();
void stop_profiler();

#endif
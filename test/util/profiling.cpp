#include "profiling.hpp"

#include <atomic>
#include <thread>
#include <chrono>

std::atomic<int> counter;
std::atomic<int> interrupt;
std::atomic_flag finished;

volatile std::size_t positive, total;

void incr_flag()
{
    counter.fetch_add(1);
}

void decr_flag()
{
    counter.fetch_sub(1);
}

std::size_t total_sample_count()
{
    return total;
}

std::size_t positive_sample_count()
{
    return positive;
}

void profile()
{
    while (interrupt.load() == 0) {
        ++total;
        int value = counter.load();
        if (value > 0)
            ++positive;
        std::this_thread::sleep_for (std::chrono::milliseconds(1));
    }
    finished.test_and_set();
}

void start_profiler()
{
    counter.store(0);
    interrupt.store(0);
    finished.clear();
    
    positive = total = 0;
    
    std::thread profiler(profile);
    profiler.detach();
}

void stop_profiler()
{
    interrupt.fetch_add(1);
    while (finished.test_and_set()) {}
}
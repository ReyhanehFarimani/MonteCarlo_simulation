/// @file rng.h
/// @brief Simple uniform random number generator wrapper around std::mt19937, for parallel system.
/// @author Reyhaneh
/// @date 2025-06-07
/// @location Ljubjana

#include "rng_parallel.h"

RNG_parallel::RNG_parallel(unsigned global_seed, int rank)
    : dist01_(0.0, 1.0)
{
    // Use seed_seq for stable, decorrelated streams
    std::seed_seq seq{global_seed, static_cast<unsigned>(rank)}; 
    engine_ = std::mt19937(seq);
}

double RNG_parallel::uniform01() {
    return dist01_(engine_);
}

double RNG_parallel::uniform(double a, double b) {
    return a + (b - a) * uniform01();
}

int RNG_parallel::uniformInt(int low, int high) {
    std::uniform_int_distribution<int> dist(low, high);
    return dist(engine_);
}

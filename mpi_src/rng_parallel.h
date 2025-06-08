/// @file rng_parallel.h
/// @brief Simple uniform random number generator wrapper around std::mt19937, for parallel systems.
/// @author Reyhaneh
/// @date 2025-06-07
/// @location Ljubljana

#ifndef RNG_PARALLEL_H
#define RNG_PARALLEL_H

#include <random>

/**
 * @class RNG_parallel
 * @brief MPI-aware, reproducible, decorrelated random number generator.
 *
 * Generates independent random streams per MPI rank using seed sequences.
 */
class RNG_parallel {
public:
    /**
     * @brief Construct RNG using global seed and rank to ensure unique streams.
     * @param global_seed Shared base seed.
     * @param rank MPI rank (ensures distinct RNG stream).
     */
    RNG_parallel(unsigned global_seed, int rank);

    /// @brief Uniform double in [0, 1)
    double uniform01();

    /// @brief Uniform double in [a, b)
    double uniform(double a, double b);

    /// @brief Uniform integer in [low, high]
    int uniformInt(int low, int high);

private:
    std::mt19937 engine_;
    std::uniform_real_distribution<double> dist01_;
};

#endif // RNG_PARALLEL_H


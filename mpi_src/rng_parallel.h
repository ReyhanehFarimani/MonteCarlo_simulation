#ifndef RNG_PARALLEL_H
#define RNG_PARALLEL_H

#include <random>
#include <cstdint>
#include <array>

/**
 * @brief Parallel-safe RNG with one independent, reproducible stream per rank.
 *
 * Design:
 *  - Fixed engine: std::mt19937_64 (deterministic across platforms).
 *  - Seeding via std::seed_seq with (global_seed, rank, fixed_mixers...) to
 *    decorrelate streams and ensure reproducibility.
 *  - uniform01(): returns double in [0,1) with full 53-bit mantissa precision.
 *
 * Usage:
 *    RNG_parallel rng(seed, rank);
 *    double x = rng.uniform01();
 */
class RNG_parallel {
public:
    RNG_parallel(unsigned global_seed, int rank)
    : dist01_(0.0, 1.0)
    {
        // Mix global seed and rank into a seed sequence to decorrelate streams.
        // Add a few fixed odd constants to improve diffusion.
        std::seed_seq seq{
            static_cast<unsigned>(global_seed),
            static_cast<unsigned>(rank),
            0x9E3779B9u, // golden ratio mix
            0x85EBCA6Bu, // from MurmurHash3
            0xC2B2AE35u  // from MurmurHash3
        };
        engine_.seed(seq);
    }

    /// Uniform in [0,1) with 53 bits of precision.
    inline double uniform01() {
        // Faster and well-defined precision vs std::uniform_real_distribution<double>
        return std::generate_canonical<double, 53>(engine_);
    }

    /// Uniform in [a,b)
    inline double uniform(double a, double b) {
        std::uniform_real_distribution<double> dist(a, b);
        return dist(engine_);
    }

    /// Normal with mean mu and stddev sigma
    inline double normal(double mu = 0.0, double sigma = 1.0) {
        std::normal_distribution<double> dist(mu, sigma);
        return dist(engine_);
    }

    /// Integer in [lo, hi] inclusive
    inline std::int64_t randint(std::int64_t lo, std::int64_t hi) {
        std::uniform_int_distribution<std::int64_t> dist(lo, hi);
        return dist(engine_);
    }

    /// Skip ahead by n draws (useful to partition work deterministically)
    inline void skip(std::uint64_t n) {
        engine_.discard(n);
    }

    /// Expose engine if needed for advanced use (e.g., STL algorithms)
    inline std::mt19937_64& engine() { return engine_; }
    inline const std::mt19937_64& engine() const { return engine_; }

private:
    std::mt19937_64 engine_;
    // Kept only to avoid re-allocation if you later want the distribution object form.
    std::uniform_real_distribution<double> dist01_;
};

#endif // RNG_PARALLEL_H

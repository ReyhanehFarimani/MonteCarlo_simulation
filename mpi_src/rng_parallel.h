#ifndef RNG_PARALLEL_H
#define RNG_PARALLEL_H

#include <random>
#include <cstdint>
#include <stdexcept>

/**
 * Parallel-safe RNG with one independent, reproducible stream per rank.
 * Engine: std::mt19937_64
 * Seeding: std::seed_seq(global_seed, rank, fixed mixers)
 * uniform01(): [0,1) with 53-bit precision.
 */
class RNG_parallel {
public:
    RNG_parallel(unsigned global_seed, int rank) {
        if (rank < 0) throw std::invalid_argument("RNG_parallel: rank must be >= 0");
        // Mix global seed and rank; add fixed odd constants to diffuse bits.
        std::seed_seq seq{
            static_cast<unsigned>(global_seed),
            static_cast<unsigned>(rank),
            0x9E3779B9u, // golden ratio
            0x85EBCA6Bu, // MurmurHash3
            0xC2B2AE35u  // MurmurHash3
        };
        engine_.seed(seq);
    }

    /// Uniform in [0,1) with 53 bits of precision.
    inline double uniform01() noexcept {
        return std::generate_canonical<double, 53>(engine_);
    }

    /// Uniform in [a,b)
    inline double uniform(double a, double b) {
        if (!(a < b)) throw std::invalid_argument("RNG_parallel::uniform: require a < b");
        std::uniform_real_distribution<double> dist(a, b);
        return dist(engine_);
    }

    /// Normal with mean mu and stddev sigma
    inline double normal(double mu = 0.0, double sigma = 1.0) {
        if (!(sigma > 0.0)) throw std::invalid_argument("RNG_parallel::normal: sigma must be > 0");
        std::normal_distribution<double> dist(mu, sigma);
        return dist(engine_);
    }

    /// Integer in [lo, hi] inclusive
    inline std::int64_t randint(std::int64_t lo, std::int64_t hi) {
        if (lo > hi) throw std::invalid_argument("RNG_parallel::randint: lo > hi");
        std::uniform_int_distribution<std::int64_t> dist(lo, hi);
        return dist(engine_);
    }

    /// Skip ahead by n draws (useful to partition work deterministically)
    inline void skip(std::uint64_t n) noexcept { engine_.discard(n); }

    /// Expose engine if needed (e.g., STL algorithms)
    inline std::mt19937_64&       engine()       noexcept { return engine_; }
    inline const std::mt19937_64& engine() const noexcept { return engine_; }

private:
    std::mt19937_64 engine_;
};

#endif // RNG_PARALLEL_H

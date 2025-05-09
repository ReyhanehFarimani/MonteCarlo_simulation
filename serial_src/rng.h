/// @file rng.h
/// @brief Simple uniform random number generator wrapper around std::mt19937.
/// @author Reyhaneh
/// @date 2025-05-09

#ifndef RNG_H
#define RNG_H

#include <random>

/**
 * @class RNG
 * @brief A simple random number generator with uniform distributions.
 *
 * Provides reproducible uniform random numbers in real and integer ranges.
 */
class RNG {
public:
    /**
     * @brief Construct the RNG with a fixed seed for reproducibility.
     *
     * @param seed The seed to initialize the underlying engine.
     */
    explicit RNG(unsigned seed);

    /**
     * @brief Generate a uniform real random number in [0, 1).
     *
     * @return A double in the interval [0, 1).
     */
    double uniform01();

    /**
     * @brief Generate a uniform real random number in [a, b).
     *
     * @param a Lower bound of the interval.
     * @param b Upper bound of the interval.
     * @return A double in the interval [a, b).
     */
    double uniform(double a, double b);

    /**
     * @brief Generate a uniform integer random number in [low, high].
     *
     * @param low Lower bound of the interval.
     * @param high Upper bound of the interval.
     * @return An integer in the interval [low, high].
     */
    int uniformInt(int low, int high);

private:
    std::mt19937 engine_; ///< Mersenne Twister engine
    std::uniform_real_distribution<double> dist01_; ///< Distribution for [0,1)
};

#endif // RNG_H

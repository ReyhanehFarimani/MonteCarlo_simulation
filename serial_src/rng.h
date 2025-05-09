// rng.h
#ifndef RNG_H
#define RNG_H

#include <random>

/// Simple uniform random number generator.
class RNG {
public:
    /// Construct with a fixed seed (for reproducibility).
    explicit RNG(unsigned seed);

    /// Uniform real in [0,1).
    double uniform01();

    /// Uniform real in [a,b).
    double uniform(double a, double b);

    /// Uniform integer in [low, high].
    int uniformInt(int low, int high);

private:
    std::mt19937 engine_;
    std::uniform_real_distribution<double> dist01_;
};

#endif // RNG_H

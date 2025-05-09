#include "rng.h"
#include <random>

RNG::RNG(unsigned seed)
    : engine_(seed), dist01_(0.0, 1.0) {
}

/// Returns a uniform real in [0,1).
double RNG::uniform01() {
    return dist01_(engine_);
}

/// Returns a uniform real in [a,b).
double RNG::uniform(double a, double b) {
    return a + (b - a) * uniform01();
}

/// Returns a uniform integer in [low, high].
int RNG::uniformInt(int low, int high) {
    std::uniform_int_distribution<int> dist(low, high);
    return dist(engine_);
}
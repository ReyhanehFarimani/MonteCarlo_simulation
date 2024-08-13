#include <vector>
#include <cmath>
#include "initial.h"

// Lennard-Jones potential energy calculation
double lennardJonesPotential(double r2) {
    const double epsilon = 1.0;
    const double sigma = 1.0;
    double r2inv = 1.0 / r2;
    double r6inv = r2inv * r2inv * r2inv;
    return 4.0 * epsilon * r6inv * (r6inv - 1.0);
}

// Weeks-Chandler-Andersen (WCA) potential energy calculation
double wcaPotential(double r2) {
    const double epsilon = 1.0;
    const double sigma = 1.0;
    const double r2cutoff = std::pow(2.0, 1.0 / 3.0) * sigma * sigma;

    if (r2 < r2cutoff) {
        double r2inv = sigma * sigma / r2;
        double r6inv = r2inv * r2inv * r2inv;
        return 4.0 * epsilon * r6inv * (r6inv - 1.0) + epsilon;
    }
    return 0.0;
}

// Yukawa potential energy calculation
double yukawaPotential(double r2) {
    const double epsilon = 1.0;
    const double kappa = 1.0;
    double r = std::sqrt(r2);
    return epsilon * exp(-kappa * r) / r;
}



enum class PotentialType {
    LennardJones,
    HardSphere,
    Yukawa
};

PotentialType selectPotentialType(const std::string &potentialName) {
    if (potentialName == "LennardJones") return PotentialType::LennardJones;
    if (potentialName == "WCA") return PotentialType::WCA;
    if (potentialName == "Yukawa") return PotentialType::Yukawa;
    throw std::invalid_argument("Unknown potential type: " + potentialName);
}

std::function<double(double)> getPotentialFunction(PotentialType type) {
    switch (type) {
        case PotentialType::LennardJones:
            return lennardJonesPotential;
        case PotentialType::WCA:
            return WCA;
        case PotentialType::Yukawa:
            return yukawaPotential;
        default:
            throw std::invalid_argument("Invalid potential type");
    }
}

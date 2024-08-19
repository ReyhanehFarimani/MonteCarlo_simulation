#include "potential.h"
#include <stdexcept>
#include <string>
#include <functional>
#include <iostream>
/**
 * @brief Calculates the Lennard-Jones potential energy.
 * 
 * @param r2 The squared distance between two particles.
 * @return The Lennard-Jones potential energy between two particles.
 */
double lennardJonesPotential(double r2) {
    const double epsilon = 1.0;
    const double sigma = 1.0;
    double r2inv = sigma * sigma / r2;
    double r6inv = r2inv * r2inv * r2inv;
    return 4.0 * epsilon * r6inv * (r6inv - 1.0);
}

/**
 * @brief Calculates the Weeks-Chandler-Andersen (WCA) potential energy.
 * 
 * @param r2 The squared distance between two particles.
 * @return The WCA potential energy between two particles.
 */
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

/**
 * @brief Calculates the Yukawa potential energy.
 * 
 * @param r2 The squared distance between two particles.
 * @return The Yukawa potential energy between two particles.
 */
double yukawaPotential(double r2) {
    const double epsilon = 1.0;
    const double kappa = 1.0;
    double r = std::sqrt(r2);
    return epsilon * exp(-kappa * r) / r;
}

/**
 * @brief Selects the potential type based on a string input.
 * 
 * @param potentialName The name of the potential type as a string.
 * @return The corresponding PotentialType enum value.
 */
PotentialType selectPotentialType(const std::string &potentialName) {
    if (potentialName == "LennardJones") return PotentialType::LennardJones;
    if (potentialName == "WCA") return PotentialType::WCA;
    if (potentialName == "Yukawa") return PotentialType::Yukawa;
    throw std::invalid_argument("Unknown potential type: " + potentialName);
}

/**
 * @brief Returns a function pointer to the appropriate potential function.
 * 
 * @param type The type of potential.
 * @return A function pointer to the selected potential function.
 */
std::function<double(double)> getPotentialFunction(PotentialType type) {
    switch (type) {
        case PotentialType::LennardJones:
            return lennardJonesPotential;
        case PotentialType::WCA:
            return wcaPotential;
        case PotentialType::Yukawa:
            return yukawaPotential;
        default:
            throw std::invalid_argument("Invalid potential type");
    }
}

#include "potential.h"
#include <stdexcept>
#include <cmath>


/**
 * @brief Calculates the Ideal potential energy.
 * 
 * @return zero!
 */
double idealPotential(){
    return 0.0;
}

/**
 * @brief Calculates the force between two ideal particles.
 * 
 * @return zero!
 */
double idealForceDotR(){
    return 0.0;
}

/**
 * @brief Calculates the Lennard-Jones potential energy.
 * 
 * @param r2 The squared distance between two particles.
 * @return The Lennard-Jones potential energy between two particles.
 */
double lennardJonesPotential(double r2) {
    const double epsilon = 1.0;
    const double sigma = 1.0;
    double r2inv = (sigma * sigma) / r2;
    double r6inv = r2inv * r2inv * r2inv;
    return 4.0 * epsilon * r6inv * (r6inv - 1.0);
}

/**
 * @brief Calculates the force between two particles using Lennard-Jones potential.
 * 
 * @param r2 The squared distance between two particles.
 * @return The magnitude of the Lennard-Jones force between two particles.
 */
double lennardJonesForceDotR(double r2) {
    const double epsilon = 1.0;
    const double sigma = 1.0;
    double r2inv = (sigma * sigma) / r2;
    double r6inv = r2inv * r2inv * r2inv;
    return 48.0 * epsilon * r6inv * (r6inv - 0.5);
}

/**
 * @brief Calculates the WCA potential energy.
 * 
 * @param r2 The squared distance between two particles.
 * @return The WCA potential energy between two particles.
 */
double wcaPotential(double r2) {
    const double epsilon = 1.0;
    const double sigma = 1.0;
    // const double r2cutoff = std::pow(2.0, 1.0 / 3.0) * sigma * sigma;
    if (r2 < 1.2599) {
        double r2inv = sigma * sigma / r2;
        double r6inv = r2inv * r2inv * r2inv;
        return 4.0 * epsilon * r6inv * (r6inv - 1.0) + epsilon;
    }
    return 0.0;
}

/**
 * @brief Calculates the force between two particles using WCA potential.
 * 
 * @param r2 The squared distance between two particles.
 * @return The magnitude of the WCA force between two particles.
 */
double wcaForceDotR(double r2) {
    const double epsilon = 1.0;
    const double sigma = 1.0;
    // const double r2cutoff = std::pow(2.0, 1.0 / 3.0) * sigma * sigma;
    if (r2 < 1.2599) {
        double r2inv = sigma * sigma / r2;
        double r6inv = r2inv * r2inv * r2inv;
        return 48.0 * epsilon * r6inv * (r6inv - 0.5);
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
 * @brief Calculates the force between two particles using Yukawa potential.
 * 
 * @param r2 The squared distance between two particles.
 * @return The magnitude of the Yukawa force between two particles.
 */
double yukawaForceDotR(double r2) {
    const double epsilon = 1.0;
    const double kappa = 1.0;
    double r = std::sqrt(r2);
    return epsilon * (kappa + 1.0 / r) * exp(-kappa * r) ;
}

/**
 * @brief Calculates the logarithmic potential energy.
 * 
 * @param f_dependant star polymer functtionality dependant part.
 * @param r2 The squared distance between two particles.
 * @return The potential energy between two star polymer.
 */
double athermalStarPotential(double r2, float f_dependant, float alpha){
    double e;
    if (r2<1) {
        double r = sqrt(r2);
        e = -log(r) + alpha;
    }
    else {
        e = alpha * exp((1 - r2)/ 2 / alpha);
    }
    return f_dependant * e;
}

/**
 * @brief Calculates the force between two star polymer cores.
 * 
 * @param f_dependant star polymer functtionality dependant part.
 * @param r2 The squared distance between two particles.
 * @return The magnitude of the force between two star polymer.
 */
double athermalStarForceDotR(double r2, float f_Dependant, float alpha){
    double f;
    if (r2<1) {
        return f_Dependant;
    }
    return f_Dependant * r2 * exp((1 - r2)/2/alpha);
}

/**
 * @brief Calculates the logarithmic potential energy.
 * 
 * @param f_dependant star polymer functtionality dependant part.
 * @param r2 The squared distance between two particles.
 * @return The potential energy between two star polymer.
 */
double thermalStarPotential(double r2, float f_dependant, float f_dependent_2, float kappa){
    double e1;
    double r = sqrt(r2);
    if (r2<1) {  
        e1 = -log(r) + 0.5;
    }
    else {
        e1 = 0.5 * exp(1 - r2);
    }
    double e2;

    e2 = f_dependent_2 * exp(-kappa * r);
    return f_dependant * e1 - e2;
}

/**
 * @brief Calculates the force between two star polymer cores.
 * 
 * @param f_dependant star polymer functtionality dependant part.
 * @param r2 The squared distance between two particles.
 * @return The magnitude of the force between two star polymer.
 */
double thermalStarForceDotR(double r2, float f_Dependant, float f_dependent_2, float kappa){
    double f;
    double r = sqrt(r2);
    f = f_dependent_2 * kappa * exp(-kappa * r) * r;
    if (r2<1) {
        return f_Dependant + f;
    }
    return f_Dependant * r2 * exp(1 - r2) - f;
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
    if (potentialName == "AthermalStar") return PotentialType::AthermalStar;
    if (potentialName == "ThermalStar") return PotentialType::ThermalStar;
    if (potentialName == "Ideal") return PotentialType::Ideal;
    throw std::invalid_argument("Unknown potential type: " + potentialName);
}

/**
 * @brief Returns the pair potential betwen two particles
 * 
 * @param r2 the squared distance between the two particles
 * @param type The type of potential
 * @param f_prime for the case of the ultrasoft potential.
 * @return The amount of potential.
 */
double computePairPotential(double r2, PotentialType potentialType, float f_prime, float f_d_prime, float kappa) {
    switch (potentialType) {
        case PotentialType::LennardJones:
            return lennardJonesPotential(r2);
        case PotentialType::WCA:
            return wcaPotential(r2);
        case PotentialType::Yukawa:
            return yukawaPotential(r2);
        case PotentialType::AthermalStar:
            return athermalStarPotential(r2, f_prime);
        case PotentialType::ThermalStar:
            return thermalStarPotential(r2, f_prime, f_d_prime, kappa);
        case PotentialType::Ideal:
            return idealPotential();
        default:
            throw std::invalid_argument("Unknown potential type");
    }
}

/**
 * @brief Returns the pair potential betwen two particles
 * 
 * @param r2 the squared distance between the two particles
 * @param type The type of potential
 * @param f_prime for the case of the ultrasoft potential.
 * @return The amount of force dot r.
 */
double computePairForce(double r2, PotentialType potentialType, float f_prime, float f_d_prime, float kappa) {
    switch (potentialType) {
        case PotentialType::LennardJones:
            return lennardJonesForceDotR(r2);
        case PotentialType::WCA:
            return wcaForceDotR(r2);
        case PotentialType::Yukawa:
            return yukawaForceDotR(r2);
        case PotentialType::AthermalStar:
            return athermalStarForceDotR(r2, f_prime);
        case PotentialType::ThermalStar:
            return thermalStarForceDotR(r2, f_prime, f_d_prime, kappa);
        case PotentialType::Ideal:
            return idealForceDotR();
        default:
            throw std::invalid_argument("Unknown potential type");
    }
}

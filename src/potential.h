#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <string>
#include <functional>
#include <cmath>

/**
 * @brief Enumeration of available potential types.
 */
enum class PotentialType {
    LennardJones,  ///< Lennard-Jones potential
    WCA,           ///< Weeks-Chandler-Andersen (WCA) potential
    Yukawa,         ///< Yukawa potential
    AthermalStar,   ///< potential between two cores of athermal star polymer.
    Ideal           ///< non interacting particls.
};

/**
 * @brief Calculates the Ideal potential energy.
 * 
 * @return zero!
 */
double idealPotential();

/**
 * @brief Calculates the force between two ideal particles.
 * 
 * @return zero!
 */
double idealForceDotR();

/**
 * @brief Calculates the Lennard-Jones potential energy.
 * 
 * @param r2 The squared distance between two particles.
 * @return The Lennard-Jones potential energy between two particles.
 */
double lennardJonesPotential(double r2);

/**
 * @brief Calculates the force between two particles using Lennard-Jones potential.
 * 
 * @param r2 The squared distance between two particles.
 * @return The magnitude of the Lennard-Jones force between two particles.
 */
double lennardJonesForceDotR(double r2);

/**
 * @brief Calculates the Weeks-Chandler-Andersen (WCA) potential energy.
 * 
 * The WCA potential is a truncated and shifted Lennard-Jones potential,
 * ensuring purely repulsive interactions.
 * 
 * @param r2 The squared distance between two particles.
 * @return The WCA potential energy between two particles.
 */
double wcaPotential(double r2);

/**
 * @brief Calculates the force between two particles using WCA potential.
 * 
 * @param r2 The squared distance between two particles.
 * @return The magnitude of the WCA force between two particles.
 */
double wcaForceDotR(double r2);

/**
 * @brief Calculates the Yukawa potential energy.
 * 
 * The Yukawa potential models screened Coulomb interactions, often used
 * in plasma physics and colloidal interactions.
 * 
 * @param r2 The squared distance between two particles.
 * @return The Yukawa potential energy between two particles.
 */
double yukawaPotential(double r2);

/**
 * @brief Calculates the force between two particles using Yukawa potential.
 * 
 * @param r2 The squared distance between two particles.
 * @return The magnitude of the Yukawa force between two particles.
 */
double yukawaForceDotR(double r2);

/**
 * @brief Calculates the logarithmic potential energy.
 * 
 * @param r2 The squared distance between two particles.
 * @param f_dependant star polymer functtionality dependant part.
 * @return The potential energy between two star polymer.
 */
double athermalStarPotential(double r2, float f_dependenat);

/**
 * @brief Calculates the force between two star polymer cores.
 * 
 * @param r2 The squared distance between two particles.
 * @param f_dependant star polymer functtionality dependant part.
 * @return The magnitude of the force between two star polymer.
 */
double athermalStarForceDotR(double r2, float f_dependant);

/**
 * @brief Selects the potential type based on a string input.
 * 
 * @param potentialName The name of the potential type as a string.
 * @return The corresponding PotentialType enum value.
 */
PotentialType selectPotentialType(const std::string &potentialName);

/**
 * @brief Returns the pair potential betwen two particles
 * 
 * @param r2 the squared distance between the two particles
 * @param type The type of potential
 * @param f_prime for the case of the ultrasoft potential.
 * @return The amount of potential.
 */
double computePairPotential(double r2, PotentialType potentialType, float f_prime);

/**
 * @brief Returns the pair potential betwen two particles
 * 
 * @param r2 the squared distance between the two particles
 * @param type The type of potential
 * @param f_prime for the case of the ultrasoft potential.
 * @return The amount of force dot r.
 */
double computePairForce(double r2, PotentialType potentialType, float f_prime);

#endif // POTENTIAL_H

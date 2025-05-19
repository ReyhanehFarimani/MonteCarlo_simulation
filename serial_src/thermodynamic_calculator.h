#ifndef THERMODYNAMIC_CALCULATOR_H
#define THERMODYNAMIC_CALCULATOR_H

#include <vector>
#include "initial.h"       ///< For Particle
#include "potential.h"     ///< For PotentialType and computePairPotential
#include "cell_list.h"

/**
 * @class ThermodynamicCalculator
 * @brief Computes state accessors and thermodynamic properties of the system.
 *
 * Provides methods to query system state (number of particles, temperature,
 * volume, density) and compute energetic and thermodynamic properties such as
 * total potential energy, local energy around a particle, virial, pressure,
 * and tail corrections in 2D.
 */
class ThermodynamicCalculator {
public:
    /**
     * @brief Construct with fixed simulation parameters.
     *
     * Internally computes r2cut_ = rcut * rcut.
     *
     * @param temperature System temperature.
     * @param potentialType Type of interaction potential.
     * @param rcut Cutoff distance for interactions.
     * @param f_prime Parameter for athermal star potential (optional).
     * @param alpha Parameter for tail strength in logarithmic potentials (optional).
     * @param f_d_prime Parameter for thermal star potential (optional).
     * @param kappa Parameter for thermal star potential (optional).
     */
    ThermodynamicCalculator(double temperature,
                             PotentialType potentialType,
                             double rcut,
                             double f_prime = 0.0,
                             double alpha   = 0.0,
                             double f_d_prime = 0.0,
                             double kappa    = 0.0);

    //── State Accessors ──────────────────────────────────────────────────
    /**
     * @brief Get the number of particles in the current state.
     * @param particles Vector of Particle objects.
     * @return Number of particles (size of vector).
     */
    size_t getNumParticles(const std::vector<Particle>& particles) const;

    /**
     * @brief Get the system temperature.
     * @return Temperature provided at construction.
     */
    double getTemperature() const;

    /**
     * @brief Compute the volume of the simulation box.
     * @param box SimulationBox describing system dimensions.
     * @return Volume of the box.
     */
    double getVolume(const SimulationBox& box) const;

    /**
     * @brief Compute the number density of the system.
     * @param particles Vector of Particle objects.
     * @param box SimulationBox object.
     * @return Density = N / V.
     */
    double computeDensity(const std::vector<Particle>& particles,
                          const SimulationBox& box) const;

    //── Energy & Virial ─────────────────────────────────────────────────
    /**
     * @brief Compute the total potential energy of all particle pairs.
     * @param particles Vector of Particle objects.
     * @param box SimulationBox for minimum-image calculations.
     * @return Total potential energy.
     */
    double computeTotalEnergy(const std::vector<Particle>& particles,
                              const SimulationBox& box) const;

    /**
     * @brief Compute the local potential energy around a given particle.
     * Uses a generic neighbor-list interface for efficiency.
     * @tparam NeighborList Any container providing iteration over neighbor indices.
     * @param particleIndex Index of the reference particle.
     * @param particles Vector of Particle objects.
     * @param box SimulationBox for distance calculations.
     * @param neighborList Iterable neighbor indices for that particle.
     * @return Local potential energy of the specified particle.
     */
    template <typename NeighborList>
    double computeLocalEnergy(size_t particleIndex,
                              const std::vector<Particle>& particles,
                              const SimulationBox& box,
                              const NeighborList& neighborList) const;

    /**
     * @brief Compute the total virial (sum of r_ij · f_ij) for pressure calculations.
     * @param particles Vector of Particle objects.
     * @param box SimulationBox for distance calculations.
     * @return Total virial term W defined via W = -1/2 Σ r_ij·f_ij.
     */
    double computeTotalVirial(const std::vector<Particle>& particles,
                              const SimulationBox& box) const;

    //── Pressure ────────────────────────────────────────────────────────
    /**
     * @brief Compute the system pressure using the 2D virial theorem.
     *
     * Implements P = ρ k_B T + W/V, where W is defined by
     * W = -1/2 Σ_{i<j} r_{ij}·f_{ij} (Allen & Tildesley, Eqn. 2.65, p. 64) and
     * PV = Nk_B T + ⟨W⟩ (Allen & Tildesley, Eqn. 2.60, p. 63) citeturn2file14.
     *
     * @param particles Vector of Particle objects.
     * @param box SimulationBox object.
     * @return System pressure.
     */
    double computePressure(const std::vector<Particle>& particles,
                           const SimulationBox& box) const;

    //── Tail Corrections (2D) ──────────────────────────────────────────
    /**
     * @brief Compute the tail correction to the energy in 2D.
     * @note Only valid for strictly two-dimensional systems.
     *       Derived from 
     *       \f$U_\text{tail} = \pi \rho \int_{r_\mathrm{cut}}^{\infty} r\,u(r)\,dr\f$ hmm 
     *       assuming \f$g(r)=1\f$ 
     *       (Allen & Tildesley, *Computer Simulation of Liquids*, Eqn. 2.74, p. 70) :contentReference[oaicite:0]{index=0}:contentReference[oaicite:1]{index=1}.
     *
     * @param particles Vector of Particle objects.
     * @param box        SimulationBox describing system dimensions.
     * @return           Energy tail correction.
     */
    double computeTailCorrectionEnergy2D(const std::vector<Particle>& particles,
        const SimulationBox& box) const;

    /**
    * @brief Compute the tail correction to the pressure in 2D.
    * @note Only valid for strictly two-dimensional systems.
    *       Derived from 
    *       \f$P_\text{tail} = \frac{\pi \rho^2}{2}\int_{r_\mathrm{cut}}^{\infty} r^2\,f(r)\,dr\f$ 
    *       assuming \f$g(r)=1\f$ 
    *       (Allen & Tildesley, *Computer Simulation of Liquids*, discussion around Eqns. 2.64–2.67) :contentReference[oaicite:2]{index=2}:contentReference[oaicite:3]{index=3}.
    *
    * @param particles Vector of Particle objects.
    * @param box        SimulationBox describing system dimensions.
    * @return           Pressure tail correction.
    */
    double computeTailCorrectionPressure2D(const std::vector<Particle>& particles,
            const SimulationBox& box) const;

    /**
     * @brief Compute total potential energy using a CellList for neighbor queries.
     * @param particles Particle positions.
     * @param box SimulationBox defining PBC.
     * @return Total potential energy (unique pairs).
     */
    double computeTotalEnergyCellList(const std::vector<Particle>& particles,
                                      const SimulationBox& box) const;

    /**
     * @brief Compute total virial using a CellList for neighbor queries.
     * @param particles Particle positions.
     * @param box SimulationBox defining PBC.
     * @return Total virial W = Σ_{i<j} r_ij·f_ij.
     */
    double computeTotalVirialCellList(const std::vector<Particle>& particles,
                                      const SimulationBox& box) const;

                                          /**
     * @brief Compute pressure using a CellList for neighbor queries.
     * @param particles Particle positions.
     * @param box SimulationBox defining PBC.
     * @return System pressure.
     */
    double computePressureCellList(const std::vector<Particle>& particles,
                                      const SimulationBox& box) const;
private:
    //── Thermodynamic parameters ───────────────────────────────────────
    double temperature_;          ///< System temperature (T)
    double rcut_;                 ///< Cutoff distance (r_cut)
    double r2cut_;                ///< Square of cutoff distance (r_cut^2)

    //── Potential-specific parameters ─────────────────────────────────
    PotentialType potentialType_; ///< Interaction potential type
    double f_prime_;              ///< Athermal star potential parameter
    double alpha_;                ///< Logarithmic/tail strength parameter
    double f_d_prime_;            ///< Thermal star potential parameter
    double kappa_;                ///< Thermal star potential range parameter
};

#endif // THERMODYNAMIC_CALCULATOR_H

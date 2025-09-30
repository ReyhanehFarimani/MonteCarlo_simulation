#ifndef THERMODYNAMIC_CALCULATOR_PARALLEL_H
#define THERMODYNAMIC_CALCULATOR_PARALLEL_H

#include <mpi.h>
#include <vector>
#include <cassert>
#include "particle.h"
#include "simulation_box.h"
#include "potential.h"
#include "cell_list_parallel.h"
#include "particle_exchange.h"

/**
 * @brief MPI/halo-aware thermodynamics using CellListParallel + ParticleExchange.
 */
class ThermodynamicCalculatorParallel {
public:
    ThermodynamicCalculatorParallel(MPI_Comm comm,
                                    double temperature,
                                    PotentialType potentialType,
                                    double rcut,
                                    double mu   = 0.0,
                                    double f    = 0.0,
                                    double alpha= 0.0,
                                    double A0   = 0.0,
                                    double kappa= 0.0);


    // ---- simple parameter getters (match serial style) ----
    double getTemperature() const { return temperature_; }
    double getActivity()   const { return mu_; }
    double getCutoff()     const { return rcut_; }
    double getCutoff2()    const { return r2cut_; }
    PotentialType getPotentialType() const { return potentialType_; }
    double getFPrime()     const { return f_prime_; }
    double getFPrimeAttraction() const { return f_d_prime_; }
    double getAlpha()      const { return alpha_; }
    double getKappa()      const { return kappa_; }

    // ---- global state accessors ----
    std::size_t numParticlesLocal(const std::vector<Particle>& owned) const;
    std::size_t numParticlesGlobal(const std::vector<Particle>& owned) const;
    double densityGlobal(const std::vector<Particle>& owned,
                         const SimulationBox& box) const;


    // ---- Energetics ----
    double totalEnergy(const std::vector<Particle>& owned,
                       const SimulationBox& box,
                       const CellListParallel& cl,
                       const ParticleExchange& px) const;

    double totalVirial(const std::vector<Particle>& owned,
                       const SimulationBox& box,
                       const CellListParallel& cl,
                       const ParticleExchange& px) const;

    double pressure(const std::vector<Particle>& owned,
                    const SimulationBox& box,
                    const CellListParallel& cl,
                    const ParticleExchange& px) const;

    // ---- Tail corrections ----
    double tailCorrectionEnergy2D(const std::vector<Particle>& owned,
                                  const SimulationBox& box) const;

    double tailCorrectionPressure2D(const std::vector<Particle>& owned,
                                    const SimulationBox& box) const;

private:
    MPI_Comm comm_;
    int rank_{0}, size_{1};

    // params
    double temperature_{0.0};
    double rcut_{0.0}, r2cut_{0.0};
    PotentialType potentialType_;
    double f_prime_{0.0}, alpha_{0.0}, f_d_prime_{0.0}, kappa_{0.0};
    double mu_{0.0};
};

#endif // THERMODYNAMIC_CALCULATOR_PARALLEL_H

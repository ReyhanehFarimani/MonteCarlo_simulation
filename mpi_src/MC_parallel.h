#ifndef MC_PARALLEL_H
#define MC_PARALLEL_H

#include <mpi.h>
#include <vector>
#include <cstdint>
#include "particle.h"
#include "simulation_box.h"
#include "cell_list_parallel.h"
#include "particle_exchange.h"
#include "thermodynamic_calculator_parallel.h"
#include "rng_parallel.h"
#include "logging_data_mpi.h"
#include "logging_traj_mpi.h"

/**
 * Minimal NVT Monte Carlo driver (MPI + halos).
 * - No global energy cache (avoid races).
 * - One random cell & one random particle per parity step.
 * - After each parity step: flush ghost updates (incremental halo).
 */
class MonteCarloNVT_MPI {
public:
    struct Params {
        double rcut;             // cutoff (must match cl/pex)
        double delta;            // max displacement magnitude
        int    out_every{0};     // frames per sweep for logging (0 = off)
        int    halo_every{0};    // optional ghost refresh every k attempts (0 = never)

        // NEW: how often to do a full rebuild+migrate+refresh.
        // Units = "attempts" (each attempt = 4 parity micro-steps).
        // 0 disables; e.g. 32 means do it every 32 attempts.
        int    rebuild_every_attempts{0};
    };


    MonteCarloNVT_MPI(MPI_Comm comm,
                      SimulationBox& box,
                      CellListParallel& cl,
                      ParticleExchange& pex,
                      ThermodynamicCalculatorParallel& thermo,
                      std::vector<Particle>& owned,
                      RNG_parallel& rng,
                      LoggingTrajMPI* traj,
                      LoggingDataMPI* data,
                      const Params& p);

    /**
     * Run nsweeps Monte Carlo sweeps. Returns global acceptance ratio.
     * A "sweep" targets ~|owned| attempted displacements per rank.
     * start_timestep is used only for logging (dump/data).
     */
    double run(std::size_t nsweeps, int start_timestep = 0);

private:
    // local energy contribution for an owned particle using unique-pair safe list
    double local_energy_of_(int i) const;

    // energy at a trial point (exclude self if it appears)
    double local_energy_of_point_(int i_moved, double x, double y) const;

    // one Metropolis attempt on owned[i]
    bool try_displacement_(int i);

    // one “parity micro-step”: pick a random cell with this parity, pick one random owned,
    // attempt move, flush ghost updates
    bool do_one_parity_step_(Parity par, long long& acc, long long& att);

    // collective seed so all ranks shuffle parities/cells consistently
    std::uint64_t bcast_seed_();

private:
    MPI_Comm   comm_;
    int        rank_{0}, size_{1};

    SimulationBox&                     box_;
    CellListParallel&                  cl_;
    ParticleExchange&                  pex_;
    ThermodynamicCalculatorParallel&   thermo_;
    std::vector<Particle>&             owned_;
    RNG_parallel&                      rng_;
    Params                             p_;

    LoggingTrajMPI*  traj_{nullptr};
    LoggingDataMPI*  data_{nullptr};

    double beta_{1.0};
};

#endif // MC_PARALLEL_H

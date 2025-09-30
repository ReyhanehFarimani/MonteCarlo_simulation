#ifndef MC_PARALLEL_H
#define MC_PARALLEL_H

#include <mpi.h>
#include <vector>
#include "particle.h"
#include "simulation_box.h"
#include "cell_list_parallel.h"
#include "particle_exchange.h"
#include "thermodynamic_calculator_parallel.h"
#include "logging_traj_mpi.h"
#include "logging_data_mpi.h"
#include "rng_parallel.h"

class MonteCarloNVT_MPI {
public:
    struct Params {
        double      delta      = 0.1;   // max component displacement in [-delta, delta]
        std::size_t halo_every = 1;     // full halo refresh every k parity blocks (>=1)
        std::size_t out_every  = 0;     // output every k sweeps (0 = off)
    };

    MonteCarloNVT_MPI(MPI_Comm comm,
                      SimulationBox& box,
                      CellListParallel& cl,
                      ParticleExchange& pex,
                      ThermodynamicCalculatorParallel& thermo,
                      std::vector<Particle>& owned,
                      RNG_parallel& rng,
                      LoggingTrajMPI* traj = nullptr,
                      LoggingDataMPI* data = nullptr,
                      const Params& p = Params{0.1, 1, 0});

    // Run 'nsweeps' MC sweeps (1 sweep = 4 parity sub-sweeps). Returns global acceptance ratio.
    double run(std::size_t nsweeps, int start_timestep = 0);

    // Optional loggers (can be changed after construction)
    void set_trajectory_logger(LoggingTrajMPI* t) noexcept { traj_ = t; }
    void set_data_logger(LoggingDataMPI* d) noexcept { data_ = d; }

private:
    // Attempt one Metropolis displacement for owned_[i]
    bool try_displacement_(int i);

    // Local energy of owned_[i] via neighbors (owned + ghosts)
    double local_energy_of_(int i) const;

    // Energy of a hypothetical point at (x,y) vs neighbors
    double local_energy_of_point_(double x, double y) const;

    // Broadcast a 64-bit seed from rank 0 to all
    std::uint64_t bcast_seed_();

private:
    MPI_Comm comm_;
    int rank_{0}, size_{1};

    SimulationBox&                     box_;
    CellListParallel&                  cl_;
    ParticleExchange&                  pex_;
    ThermodynamicCalculatorParallel&   thermo_;
    std::vector<Particle>&             owned_;
    RNG_parallel&                      rng_;
    Params                             p_;
    double                             beta_{1.0};

    LoggingTrajMPI*  traj_{nullptr};
    LoggingDataMPI*  data_{nullptr};
};

#endif // MC_PARALLEL_H

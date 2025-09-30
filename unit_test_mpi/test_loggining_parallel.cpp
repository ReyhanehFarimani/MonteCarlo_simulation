// unit_test_mpi/test_logging_parallel.cpp
#include "catch.hpp"
#include <mpi.h>
#include <vector>

#include "../mpi_src/particle.h"
#include "../mpi_src/simulation_box.h"
#include "../mpi_src/cell_list_parallel.h"
#include "../mpi_src/particle_exchange.h"
#include "../mpi_src/potential.h"
#include "../mpi_src/thermodynamic_calculator_parallel.h"
#include "../mpi_src/rng_parallel.h"
#include "../mpi_src/initial_mpi.h"
#include "../mpi_src/logging_traj_mpi.h"
#include "../mpi_src/logging_data_mpi.h"   // the scalar logger we added

TEST_CASE("Log 10 random configurations: trajectory + thermodynamics", "[MPI][logging][trajectory][thermo]") {
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // --- Geometry & interaction parameters (safe for any Px*Py decomposition) ---
    const double Lx = 160.0, Ly = 160.0;       // large enough so local Lloc >= 2*rcut
    const double rcut = 3.0;                   // halo = one cell, OK for our cell list
    const int    N_global = 2000;              // total particles across all ranks

    const double T = 1.0;
    const PotentialType pot = PotentialType::LennardJones;
    const double mu = 0.0, f = 0.0, alpha = 0.0, A0 = 0.0, kappa = 1.0; // unused for LJ

    SimulationBox box(Lx, Ly);
    auto decomp = box.bestDecomposition(size);

    // --- Loggers ---
    LoggingTrajMPI traj("demo_run", LoggingTrajMPI::Mode::MPIIO, MPI_COMM_WORLD, /*append=*/false, /*rank_as_type=*/true);
    LoggingDataMPI data("demo_run", MPI_COMM_WORLD, /*append=*/false);

    // --- Thermo calc (parallel) ---
    ThermodynamicCalculatorParallel thermo(MPI_COMM_WORLD, T, pot, rcut, mu, f, alpha, A0, kappa);

    // We’ll generate 10 independent random configurations (different seeds)
    const int nframes = 10;
    const unsigned base_seed = 424242;

    for (int frame = 0; frame < nframes; ++frame) {
        // 1) Local initialization (random uniform in each subdomain)
        std::vector<Particle> owned;
        initializeParticles_globalN(owned, box, decomp, rank, N_global,
                                    /*random=*/true, /*seed=*/base_seed + frame,
                                    /*id_stride=*/1000000LL);

        // 2) Cell list + ghosts
        CellListParallel::Params clp_params{rcut};
        CellListParallel clp(box, decomp, rank, clp_params);
        clp.buildInterior(owned);

        ParticleExchange::Params pex_params{/*enable_incremental=*/false};
        ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pex_params);
        pex.refreshGhosts(owned, clp);

        // 3) Log one trajectory frame (GLOBAL positions, z=0)
        traj.log_dump(owned, box, /*timestep=*/frame);   // LAMMPS dump
        traj.log_xyz (owned, box);                       // XYZ (optional)

        // 4) Log thermodynamic scalars (rank 0 CSV)
        data.log_step(owned, box, clp, pex, thermo, /*timestep=*/frame);

        // (Optional) small deterministic displacement to make the trajectory non-static.
        // Commented out because user asked for independent random initial configurations:
        // for (auto& p : owned) { p.updatePosition(0.05, -0.03); box.applyPBC(p); }
        // clp.buildInterior(owned); pex.refreshGhosts(owned, clp);
    }

    // Close logs (safe to call more than once)
    traj.close();
    data.close();

    // Minimal sanity: every rank did produce its portion; nothing to REQUIRE about files here.
    SUCCEED();
}

// sanity_check/main_parallel.cpp
#include <mpi.h>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>

#include "../mpi_src/particle.h"
#include "../mpi_src/simulation_box.h"
#include "../mpi_src/cell_list_parallel.h"
#include "../mpi_src/particle_exchange.h"
#include "../mpi_src/thermodynamic_calculator_parallel.h"
#include "../mpi_src/initial_mpi.h"
#include "../mpi_src/logging_traj_mpi.h"

// Simple CSV writer on rank 0
static void write_csv_header_if_new(std::ofstream& ofs) {
    ofs << "frame,N_global,rho,Energy,Virial,Pressure\n";
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // ---- Config ----
    const int    frames    = 200;            // number of initial configs
    const int    N_global  = 24000;           // total particles
    const double Lx        = 1000.0;
    const double Ly        = 800.0;      // temperature for pressure
    const unsigned seed0   = 12345;
    // --- choose a softer potential for sanity check ---
    const PotentialType pot = PotentialType::AthermalStar;
    const double T     = 1.0;
    const double rcut  = 5.0;
    const double mu    = 0.0;      // unused here
    const double f     = 8.0;      // star functionality factor
    const double alpha = 0.3;      // tail strength
    const double A0    = 0.0;      // (thermal term off)
    const double kappa = 0.0;      // (thermal term off)

    ThermodynamicCalculatorParallel thermo(
        MPI_COMM_WORLD, T, pot, rcut, mu, f, alpha, A0, kappa
    );

    // ---- Box & decomp ----
    SimulationBox box(Lx, Ly);
    auto decomp = box.bestDecomposition(size);

    // Choose rcut so the local extent condition holds: L_loc >= 2*rcut
    double x0,x1,y0,y1; decomp.localBounds(rank, Lx, Ly, x0, x1, y0, y1);
    const double Llocx = x1 - x0, Llocy = y1 - y0;
    // const double rcut  = 0.4 * std::min(Llocx, Llocy);  // conservative
    // (This is large for a demo; for real sims you’ll pick 2.5 etc. and ensure decomp fits)

    // ---- Logging: gather everything to rank 0 so we get a single XYZ ----
    LoggingTrajMPI traj("sanity_parallel_traj", LoggingTrajMPI::Mode::GatherRoot,
                        MPI_COMM_WORLD, /*append=*/false, /*rank_as_type=*/true);

    // Rank 0 data CSV (thermo)
    std::ofstream csv;
    if (rank == 0) {
        csv.open("sanity_parallel_data.csv", std::ios::out | std::ios::trunc);
        write_csv_header_if_new(csv);
    }


    // ---- Loop over independent initial configurations ----
    for (int frame = 0; frame < frames; ++frame) {
        // 1) Local initialization (random, reproducible stream per rank)
        std::vector<Particle> owned;
        initializeParticles_globalN(owned, box, decomp, rank,
                                    N_global, /*random=*/true,
                                    /*seed=*/seed0 + frame, /*id_stride=*/1000000LL);

        // 2) Build cell list + halos
        CellListParallel::Params clp_params{rcut};
        CellListParallel clp(box, decomp, rank, clp_params);
        clp.buildInterior(owned);

        ParticleExchange::Params pex_params{/*enable_incremental=*/false};
        ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pex_params);
        pex.refreshGhosts(owned, clp);

        // 3) Thermo (global)
        const double U = thermo.totalEnergy(owned, box, clp, pex);
        const double W = thermo.totalVirial(owned, box, clp, pex);
        const double P = thermo.pressure(owned, box, clp, pex);
        const double rho = thermo.densityGlobal(owned, box);
        const std::size_t Ng = thermo.numParticlesGlobal(owned);

        // 4) Log positions (single XYZ file with 10 frames) — GLOBAL coordinates
        traj.log_xyz(owned, box);

        // 5) Rank-0 writes CSV
        if (rank == 0) {
            csv << frame << "," << Ng << "," << std::setprecision(17) << rho << ","
                << U << "," << W << "," << P << "\n";
        }
    }

    if (rank == 0) csv.close();
    traj.close();
    MPI_Finalize();
    return 0;
}

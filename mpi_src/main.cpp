// file: mpi_src/main_mc_parallel.cpp
#include <mpi.h>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>

#include "input.h"
#include "initial_mpi.h"  // initializeParticles / initializeParticles_from_file
#include "particle.h"
#include "simulation_box.h"
#include "cell_list_parallel.h"
#include "particle_exchange.h"
#include "potential.h"
#include "thermodynamic_calculator_parallel.h"
#include "logging_traj_mpi.h"
#include "logging_data_mpi.h"
#include "rng_parallel.h"
#include "MC_parallel.h"

// Helper: same behavior as serial getOrDefault (Input already returns 0.0 if missing)
static inline double getOrDefault(const Input& input, const std::string& key) {
    return input.getConstant(key);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <input_file>\n";
        }
        MPI_Finalize();
        return 1;
    }

    std::string inputFile = argv[1];

    try {
        // ---------------- Parse input ----------------
        Input input(inputFile);

        // Required
        double Lx = input.getConstant("Lx");
        double Ly = input.getConstant("Ly");
        int    Nglob = static_cast<int>(input.getConstant("N"));
        double rcut  = input.getConstant("rcut");
        double T     = input.getConstant("T");
        double nSteps        = input.getConstant("nSteps");
        double eSteps        = input.getConstant("eSteps");
        double outputFreq    = input.getConstant("outputFreq");
        double cellUpdateFreq= input.getConstant("cellUpdateFreq"); // not used directly here; MC handles halo
        std::string potentialName = input.getFilename("potential");
        std::string out_xyz  = input.getFilename("out_xyz");   // used to derive basename
        std::string out_data = input.getFilename("out_data");  // ditto
        double delta = input.getConstant("delta");

        // Optional
        double mu    = getOrDefault(input, "mu");     // unused in NVT, accepted for completeness
        double f     = getOrDefault(input, "f");
        double alpha = getOrDefault(input, "alpha");
        double A0    = getOrDefault(input, "A_0");
        double kappa = getOrDefault(input, "kappa");
        int    seed  = static_cast<int>(getOrDefault(input, "seed"));
        std::string position_file = input.getFilename("position_file");
        std::string ensemble_name = input.getFilename("ensemble");
        double pressure = getOrDefault(input, "P");   // unused in NVT
        double delta_V  = getOrDefault(input, "delta_V"); // unused in NVT

        // Ensemble gate: MPI impl supports only NVT here
        if (ensemble_name != "NVT" && !ensemble_name.empty()) {
            if (rank == 0) {
                std::cerr << "Error: MPI driver currently supports only NVT. Got '" << ensemble_name << "'.\n";
            }
            MPI_Finalize();
            return 1;
        }

        // -------------- Summarize (rank 0) --------------
        if (rank == 0) {
            std::cout << "\n[Input Summary / MPI]\n";
            std::cout << "Ranks = " << size << "\n";
            std::cout << "Lx = " << Lx << ", Ly = " << Ly << ", N_global = " << Nglob << "\n";
            std::cout << "rcut = " << rcut << ", T = " << T << "\n";
            std::cout << "eSteps = " << eSteps << ", nSteps = " << nSteps << "\n";
            std::cout << "outputFreq = " << outputFreq << ", cellUpdateFreq = " << cellUpdateFreq << "\n";
            std::cout << "potential = " << potentialName << "\n";
            std::cout << "mu = " << mu << ", f = " << f << ", alpha = " << alpha
                      << ", A_0 = " << A0 << ", kappa = " << kappa << "\n";
            std::cout << "seed = " << seed << "\n";
            std::cout << "position_file = " << (position_file.empty() ? "[none]" : position_file) << "\n";
            std::cout << "ensemble = " << (ensemble_name.empty() ? "NVT (default)" : ensemble_name) << "\n";
            std::cout << "Output XYZ = " << out_xyz << ", Output Data = " << out_data << "\n\n";
        }

        // -------------- Geometry & decomposition --------------
        SimulationBox box(Lx, Ly);
        auto decomp = box.bestDecomposition(size);

        // derive this rank's bounds
        double x0,x1,y0,y1;
        decomp.localBounds(rank, Lx, Ly, x0, x1, y0, y1);
        const double Llocx = x1 - x0, Llocy = y1 - y0;

        // Cap rcut so each tile can build >= 2 cells per direction (matches your tests)
        const double cap = 0.49 * std::min(Llocx, Llocy);
        rcut = std::max(1e-8, std::min(rcut, cap));

        // -------------- Build global config (deterministic) and filter owned --------------
        std::vector<Particle> global;
        if (!position_file.empty()) {
            // Use the same helper as serial (reads global coords)
            std::cerr<<"system yet have not this option!"<<std::endl;
            return 1;
        } else {
            initializeParticles_globalN(
                global,           // std::vector<Particle>&
                box,                 // const SimulationBox&
                decomp,              // const SimulationBox::Decomposition&
                rank,                // this rank
                /*N_global=*/Nglob,      // global particle count
                /*random=*/true,     // or false if you have a deterministic layout
                static_cast<unsigned>(seed), // RNG seed
                /*id_stride=*/size   // ensures unique IDs across ranks (id = base + k*size)
            );

        }

        std::vector<Particle> owned;
        owned.reserve(global.size());
        for (const auto& p : global) {
            if (p.x >= x0 && p.x < x1 && p.y >= y0 && p.y < y1) owned.push_back(p);
        }

        // -------------- Parallel infra: cell list + ghosts --------------
        CellListParallel::Params clp_params{rcut};
        CellListParallel cl(box, decomp, rank, clp_params);
        cl.buildInterior(owned);

        ParticleExchange::Params pex_params{/*enable_incremental=*/true};
        ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, cl, pex_params);
        pex.refreshGhosts(owned, cl);

        // -------------- Potential & Thermodynamics --------------
        PotentialType pot = selectPotentialType(potentialName);
        ThermodynamicCalculatorParallel thermo(MPI_COMM_WORLD, T, pot, rcut,
                                               /*mu*/mu, /*f*/f, /*alpha*/alpha, /*A0*/A0, /*kappa*/kappa);

        // -------------- Logging --------------
        // Use out_xyz/data to form a single basename (strip extension if any)
        auto strip_ext = [](const std::string& s)->std::string {
            auto pos = s.rfind('.');
            return (pos == std::string::npos) ? s : s.substr(0, pos);
        };
        std::string basename = !out_xyz.empty() ? strip_ext(out_xyz)
                              : (!out_data.empty() ? strip_ext(out_data)
                                                   : std::string("run"));
        // Per-rank trajectory files by default (robust and simple)
        LoggingTrajMPI traj(basename, LoggingTrajMPI::Mode::MPIIO,
                            MPI_COMM_WORLD, /*append=*/false, /*rank_as_type=*/true);
        LoggingDataMPI data(basename, MPI_COMM_WORLD, /*append=*/false);

        // -------------- RNG (per-rank independent stream) --------------
        RNG_parallel rng(static_cast<unsigned>(seed), rank);

        // -------------- MC driver params & construction --------------
        MonteCarloNVT_MPI::Params mp;
        mp.delta      = delta;
        mp.halo_every = static_cast<int>(1/delta);                               // keep halos fresh (uses incremental queue + periodic flush)
        mp.out_every  = static_cast<int>(outputFreq);    // write thermo/xyz every 'outputFreq' sweeps (0 disables)
        mp.rebuild_every_attempts = static_cast<int>(cellUpdateFreq);
        MonteCarloNVT_MPI mc(MPI_COMM_WORLD, box, cl, pex, thermo, owned, rng,
                             &traj, &data, mp);

        // -------------- Run: equilibration then production --------------
        if (eSteps > 0) {
            if (rank == 0) std::cout << "Equilibration (NVT, MPI)…\n";
            (void)mc.run(static_cast<std::size_t>(eSteps));
        }

        if (nSteps > 0) {
            if (rank == 0) std::cout << "Production (NVT, MPI)…\n";
            (void)mc.run(static_cast<std::size_t>(nSteps));
        }

        // -------------- Close logs --------------
        traj.close();
        data.close();

        if (rank == 0) std::cout << "Simulation completed.\n";
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << "Error: " << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 2);
        return 2;
    }

    MPI_Finalize();
    return 0;
}

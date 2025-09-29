#include "catch.hpp"
#include <mpi.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <numeric>
#include <cstdio>
#include <iomanip>

#include "../mpi_src/particle.h"
#include "../mpi_src/simulation_box.h"
#include "../mpi_src/initial_mpi.h"
#include "../mpi_src/logging_traj_mpi.h"

// ---------------- Helpers ----------------

static bool file_exists(const std::string& f) {
    std::ifstream in(f);
    return in.good();
}

static void remove_if_exists(const std::string& f) {
    if (file_exists(f)) std::remove(f.c_str());
}

static size_t file_count_lines(const std::string& f) {
    std::ifstream in(f);
    size_t n = 0;
    std::string s;
    while (std::getline(in, s)) ++n;
    return n;
}

// ---------------- Tests ----------------

TEST_CASE("LoggingTrajMPI: GatherRoot dump rank-as-type; global coords within box", "[MPI][log][dump]") {
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    SimulationBox box(110.0, 75.5);
    auto decomp = box.bestDecomposition(size);

    const int N_global = 25700 + size;
    std::vector<Particle> local;
    initializeParticles_globalN(local, box, decomp, rank, N_global, /*random=*/true, /*seed=*/2025);

    const std::string base = "test_traj_gather";
    if (rank == 0) remove_if_exists(base + ".dump");
    MPI_Barrier(MPI_COMM_WORLD);

    {
        LoggingTrajMPI logger(base, LoggingTrajMPI::Mode::GatherRoot, MPI_COMM_WORLD,
                              /*append=*/false, /*rank_as_type=*/true);
        logger.log_dump(local, box, /*timestep=*/42);
        logger.close();
    }

    if (rank == 0) {
        const std::string fname = base + ".dump";
        REQUIRE(file_exists(fname));

        std::ifstream in(fname);
        REQUIRE(in.good());

        std::string line;
        std::getline(in, line); REQUIRE(line == "ITEM: TIMESTEP");
        std::getline(in, line); REQUIRE(std::stoi(line) == 42);
        std::getline(in, line); REQUIRE(line == "ITEM: NUMBER OF ATOMS");
        std::getline(in, line); const int Ntot = std::stoi(line);
        REQUIRE(Ntot == N_global);
        std::getline(in, line); REQUIRE(line.find("ITEM: BOX BOUNDS") == 0);
        std::getline(in, line); // x bounds
        std::getline(in, line); // y bounds
        std::getline(in, line); // z bounds
        std::getline(in, line); REQUIRE(line == "ITEM: ATOMS id type x y z");

        std::vector<int> type_seen(size, 0);
        int count = 0;
        const double Lx = box.getLx(), Ly = box.getLy();
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            int id, type; double x, y, z;
            REQUIRE(static_cast<bool>(iss >> id >> type >> x >> y >> z));
            REQUIRE(type >= 1);
            REQUIRE(type <= size);
            type_seen[type-1] += 1;
            REQUIRE(x >= 0.0 - 1e-9);
            REQUIRE(x <= Lx + 1e-9);
            REQUIRE(y >= 0.0 - 1e-9);
            REQUIRE(y <= Ly + 1e-9);
            REQUIRE(z == Approx(0.0));
            ++count;
        }
        REQUIRE(count == N_global);
        for (int r = 0; r < size; ++r) REQUIRE(type_seen[r] >= 1);

        // optional cleanup:
        // std::remove(fname.c_str());
    }
}

TEST_CASE("LoggingTrajMPI: PerRank xyz & dump basic write", "[MPI][log][perrank]") {
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    SimulationBox box(200.0, 160.0);
    auto decomp = box.bestDecomposition(size);

    std::vector<Particle> local;
    initializeParticles_globalN(local, box, decomp, rank, /*N_global=*/11*size + 3, /*random=*/true);

    std::string base = "test_traj_rank";

    std::ostringstream xyz;  xyz << base << ".rank" << std::setfill('0') << std::setw(4) << rank << ".xyz";
    std::ostringstream dump; dump << base << ".rank" << std::setfill('0') << std::setw(4) << rank << ".dump";

    // Each rank removes its own output files first
    remove_if_exists(xyz.str());
    remove_if_exists(dump.str());
    MPI_Barrier(MPI_COMM_WORLD);

    {
        LoggingTrajMPI logger(base, LoggingTrajMPI::Mode::PerRank, MPI_COMM_WORLD,
                              /*append=*/false, /*rank_as_type=*/true);
        logger.log_xyz(local, box);
        logger.log_dump(local, box, /*timestep=*/7);
        logger.close();
    }

    REQUIRE(file_exists(xyz.str()));
    REQUIRE(file_exists(dump.str()));

    const size_t xyz_lines  = file_count_lines(xyz.str());
    REQUIRE(xyz_lines == 2 + local.size());

    const size_t dump_lines = file_count_lines(dump.str());
    REQUIRE(dump_lines == 9 + local.size()); // dump header(8+1) + atoms

    // optional cleanup:
    // std::remove(xyz.str().c_str());
    // std::remove(dump.str().c_str());
}

TEST_CASE("LoggingTrajMPI: MPI-IO dump rank-as-type; structure valid", "[MPI][log][mpiio][dump]") {
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);

    SimulationBox box(180.0, 123.0);
    auto decomp = box.bestDecomposition(size);

    const int N_global = 3330 + size; // ensure non-uniform remainder
    std::vector<Particle> local;
    initializeParticles_globalN(local, box, decomp, rank, N_global, /*random=*/true, /*seed=*/2026);

    const std::string base = "test_traj_mpiio_dump";
    if (rank == 0) remove_if_exists(base + ".dump");
    MPI_Barrier(MPI_COMM_WORLD);

    {
        LoggingTrajMPI logger(base, LoggingTrajMPI::Mode::MPIIO, MPI_COMM_WORLD,
                              /*append=*/false, /*rank_as_type=*/true);
        logger.log_dump(local, box, /*timestep=*/99);
        logger.close();
    }
    MPI_Barrier(MPI_COMM_WORLD); // ensure file is closed before reading

    if (rank == 0) {
        const std::string fname = base + ".dump";
        REQUIRE(file_exists(fname));

        std::ifstream in(fname);
        REQUIRE(in.good());

        std::string line;
        std::getline(in, line); REQUIRE(line == "ITEM: TIMESTEP");
        std::getline(in, line); REQUIRE(std::stoi(line) == 99);
        std::getline(in, line); REQUIRE(line == "ITEM: NUMBER OF ATOMS");
        std::getline(in, line); const int Ntot = std::stoi(line);
        REQUIRE(Ntot == N_global);
        std::getline(in, line); REQUIRE(line.find("ITEM: BOX BOUNDS") == 0);
        std::getline(in, line); // x bounds
        std::getline(in, line); // y bounds
        std::getline(in, line); // z bounds
        std::getline(in, line); REQUIRE(line == "ITEM: ATOMS id type x y z");

        std::vector<int> type_seen(size, 0);
        int count = 0;
        const double Lx = box.getLx(), Ly = box.getLy();
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            int id, type; double x, y, z;
            REQUIRE(static_cast<bool>(iss >> id >> type >> x >> y >> z));
            REQUIRE(type >= 1);
            REQUIRE(type <= size);
            type_seen[type-1] += 1;
            REQUIRE(x >= 0.0 - 1e-9);
            REQUIRE(x <= Lx + 1e-9);
            REQUIRE(y >= 0.0 - 1e-9);
            REQUIRE(y <= Ly + 1e-9);
            REQUIRE(z == Approx(0.0));
            ++count;
        }
        REQUIRE(count == N_global);
        for (int r = 0; r < size; ++r) REQUIRE(type_seen[r] >= 1);
    }
}

TEST_CASE("LoggingTrajMPI: MPI-IO xyz single frame", "[MPI][log][mpiio][xyz]") {
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);

    SimulationBox box(140.0, 140.0);
    auto decomp = box.bestDecomposition(size);

    const int N_global = 200 + 2*size;
    std::vector<Particle> local;
    initializeParticles_globalN(local, box, decomp, rank, N_global, /*random=*/false);

    const std::string base = "test_traj_mpiio_xyz";
    if (rank == 0) remove_if_exists(base + ".xyz");
    MPI_Barrier(MPI_COMM_WORLD);

    {
        LoggingTrajMPI logger(base, LoggingTrajMPI::Mode::MPIIO, MPI_COMM_WORLD,
                              /*append=*/false, /*rank_as_type=*/false);
        logger.log_xyz(local, box);
        logger.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        const std::string fname = base + ".xyz";
        REQUIRE(file_exists(fname));

        std::ifstream in(fname);
        REQUIRE(in.good());

        int Nread = -1;
        in >> Nread;
        REQUIRE(Nread == N_global);
        std::string comment;
        std::getline(in, comment); // rest of first line
        std::getline(in, comment); // second line (comment)
        REQUIRE(comment.find("Lx=") != std::string::npos);

        size_t nlines = 0;
        std::string line;
        while (std::getline(in, line)) if (!line.empty()) ++nlines;
        REQUIRE(nlines == static_cast<size_t>(N_global));
    }
}

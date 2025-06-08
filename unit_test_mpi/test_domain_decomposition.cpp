#include "catch.hpp"
#include "../mpi_src/domain_decomposition.h"
#include "../mpi_src/domain_decomposition.h"
#include <mpi.h>
#include <vector>
#include <iostream>

TEST_CASE("DomainDecomposition: grid dimensions and rank mapping", "[mpi]") {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    double Lx = 8.0;
    double Ly = 6.0;
    SimulationBox globalBox(Lx, Ly);
    double rcut = 1.0;

    DomainDecomposition dd(globalBox, rcut, comm);

    int Px = dd.getPx();
    int Py = dd.getPy();
    REQUIRE(Px * Py == size);

    int coordX = dd.getCoordX();
    int coordY = dd.getCoordY();

    int coords[2] = {coordX, coordY};
    int expectedRank;
    MPI_Cart_rank(dd.getCartComm(), coords, &expectedRank);
    REQUIRE(rank == expectedRank);
}

TEST_CASE("DomainDecomposition: local box size", "[mpi]") {
    MPI_Comm comm = MPI_COMM_WORLD;
    double Lx = 8.0;
    double Ly = 6.0;
    SimulationBox globalBox(Lx, Ly);
    double rcut = 1.0;

    DomainDecomposition dd(globalBox, rcut, comm);
    const SimulationBox& localBox = dd.getLocalBox();

    REQUIRE(localBox.getLx() == Approx(Lx / dd.getPx()));
    REQUIRE(localBox.getLy() == Approx(Ly / dd.getPy()));
}

TEST_CASE("DomainDecomposition: neighbor ranks valid and complete", "[mpi]") {
    MPI_Comm comm = MPI_COMM_WORLD;
    DomainDecomposition dd(SimulationBox(8.0, 6.0), 1.0, comm);
    int size = dd.getSize();

    std::vector<int> neighbors = dd.getNeighborRanks();
    REQUIRE(neighbors.size() == 8);

    for (int r : neighbors) {
        REQUIRE(r >= 0);
        REQUIRE(r < size);
    }
}

TEST_CASE("DomainDecomposition: debug rank info", "[mpi]") {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    DomainDecomposition dd(SimulationBox(8.0, 6.0), 1.0, comm);

    if (rank == 0) {
        std::cout << "[Debug] Domain decomposition " << dd.getPx() << " x " << dd.getPy() << "\n";
    }
    MPI_Barrier(comm);
    std::cout << "Rank " << rank
              << " coord = (" << dd.getCoordX() << ", " << dd.getCoordY() << ")"
              << ", neighbors = ";
    for (int n : dd.getNeighborRanks()) std::cout << n << " ";
    std::cout << std::endl;
    MPI_Barrier(comm);
}
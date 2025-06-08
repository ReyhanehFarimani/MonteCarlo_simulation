#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include <mpi.h>
#include <iostream>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "\n[INFO] Running RNG_parallel MPI tests with Catch2...\n" << std::endl;
    }

    int result = Catch::Session().run(argc, argv);
    std::cout << "[Rank " << rank << "] Starting test run\n" << std::flush;


    MPI_Finalize();
    return result;
}

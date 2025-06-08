#include "catch.hpp"
#include "../mpi_src/rng_parallel.h"
#include <mpi.h>
#include <vector>
#include <iostream>

/**
 * @brief This test ensures:
 * - Each rank gets a different random stream (first values differ).
 * - Reproducibility: same stream on each rank when re-initialized.
 */
TEST_CASE("RNG_parallel reproducibility and uniqueness", "[RNG_parallel][MPI]") {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const unsigned seed = 2025;

    RNG_parallel rng1(seed, rank);
    RNG_parallel rng2(seed, rank);  // same seed and rank = same stream

    std::vector<double> seq1(2025);
    std::vector<double> seq2(2025);
    for (int i = 0; i < 2025; ++i) {
        seq1[i] = rng1.uniform01();
        seq2[i] = rng2.uniform01();
    }

    // ✅ Reproducibility
    REQUIRE(seq1 == seq2);

    // ✅ Uniqueness across ranks
    double my_first = seq1[0];
    std::vector<double> all_first(size);
    MPI_Gather(&my_first, 1, MPI_DOUBLE, all_first.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                REQUIRE(all_first[i] != Approx(all_first[j]));
            }
        }
    }
}

TEST_CASE("RNG_parallel inter-rank correlation is low", "[RNG_parallel][correlation][MPI]") {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const unsigned seed = 2025;
    RNG_parallel rng(seed, rank);

    const int N = 10000;
    std::vector<double> my_seq(N);
    for (int i = 0; i < N; ++i) {
        my_seq[i] = rng.uniform01();
    }

    // Gather all sequences to rank 0
    std::vector<double> all_seqs;
    if (rank == 0) {
        all_seqs.resize(N * size);
    }

    MPI_Gather(my_seq.data(), N, MPI_DOUBLE,
               all_seqs.data(), N, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    if (rank == 0) {
        auto get_seq = [&](int r) -> std::vector<double> {
            return std::vector<double>(all_seqs.begin() + r * N, all_seqs.begin() + (r + 1) * N);
        };

        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                auto xi = get_seq(i);
                auto xj = get_seq(j);

                // Compute Pearson correlation coefficient
                double mean_x = std::accumulate(xi.begin(), xi.end(), 0.0) / N;
                double mean_y = std::accumulate(xj.begin(), xj.end(), 0.0) / N;

                double num = 0.0, den_x = 0.0, den_y = 0.0;
                for (int k = 0; k < N; ++k) {
                    double dx = xi[k] - mean_x;
                    double dy = xj[k] - mean_y;
                    num += dx * dy;
                    den_x += dx * dx;
                    den_y += dy * dy;
                }

                double corr = num / std::sqrt(den_x * den_y);
                INFO("Correlation between rank " << i << " and rank " << j << " = " << corr);
                std::cout << "[INFO] Correlation between rank " << i << " and rank " << j << " = " << corr << std::endl;
                REQUIRE(std::fabs(corr) < 0.05);  // Allow weak noise
            }
        }
    }
}


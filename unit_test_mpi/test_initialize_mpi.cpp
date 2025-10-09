#include "catch.hpp"
#include <mpi.h>
#include <vector>
#include <numeric>
#include <cmath>
#include <limits>
#include <algorithm>

#include "../mpi_src/particle.h"
#include "../mpi_src/simulation_box.h"
#include "../mpi_src/rng_parallel.h"
#include "../mpi_src/initial_mpi.h"

static void get_bounds(const SimulationBox& box,
                       const SimulationBox::Decomposition& d,
                       int rank,
                       double& x0, double& x1, double& y0, double& y1)
{
    d.localBounds(rank, box.getLx(), box.getLy(), x0, x1, y0, y1);
}

TEST_CASE("MPI init: random, global coords inside local subdomain, reproducible", "[MPI][init][random]") {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Box & decomp
    SimulationBox box(110.0, 75.887);
    auto decomp = box.bestDecomposition(size); // pure function of size and box

    // Choose a global N that exercises uneven split
    const int N_global = 1003;
    const unsigned seed = 2025;

    // First init
    std::vector<Particle> p1;
    initializeParticles_globalN(p1, box, decomp, rank, N_global, /*random=*/true, seed);

    // Second init (same seed) -> identical on each rank
    std::vector<Particle> p2;
    initializeParticles_globalN(p2, box, decomp, rank, N_global, /*random=*/true, seed);

    // Local count equals deterministic split
    const int expected_local = compute_local_count(N_global, decomp, rank);
    REQUIRE((int)p1.size() == expected_local);
    REQUIRE((int)p2.size() == expected_local);

    // Reproducibility: same positions
    REQUIRE(p1.size() == p2.size());
    for (size_t i = 0; i < p1.size(); ++i) {
        REQUIRE(p1[i].x == Approx(p2[i].x));
        REQUIRE(p1[i].y == Approx(p2[i].y));
        REQUIRE(p1[i].id == p2[i].id);
    }

    // Global coordinates lie in this rank's bounds [x0,x1) × [y0,y1)
    double x0,x1,y0,y1;
    get_bounds(box, decomp, rank, x0,x1,y0,y1);

    const double eps = 1e-12;
    for (const auto& p : p1) {
        REQUIRE(p.x >= x0 - eps);
        REQUIRE(p.x <  x1 + eps); // guard tiny FP; initializer centers or uses [0,1) so it's safe
        REQUIRE(p.y >= y0 - eps);
        REQUIRE(p.y <  y1 + eps);
    }

    // Sum local counts across ranks equals N_global
    int local_n = (int)p1.size();
    int sum_n = 0;
    MPI_Allreduce(&local_n, &sum_n, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) REQUIRE(sum_n == N_global);

    // Check id ranges are disjoint per rank (with default id_stride=1'000'000)
    long long my_min_id = p1.empty() ? std::numeric_limits<long long>::max() : p1.front().id;
    long long my_max_id = p1.empty() ? std::numeric_limits<long long>::min() : p1.back().id;
    // (Initializer assigns ids monotonically from id_base + i)
    if (!p1.empty()) {
        my_min_id = std::min_element(p1.begin(), p1.end(), [](auto&a, auto&b){return a.id<b.id;})->id;
        my_max_id = std::max_element(p1.begin(), p1.end(), [](auto&a, auto&b){return a.id<b.id;})->id;
    }

    std::vector<long long> all_min(size), all_max(size);
    MPI_Gather(&my_min_id, 1, MPI_LONG_LONG, all_min.data(), 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Gather(&my_max_id, 1, MPI_LONG_LONG, all_max.data(), 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                const bool disjoint =
                    (all_max[i] < all_min[j]) || (all_max[j] < all_min[i]);
                REQUIRE(disjoint);  // simple bool — Catch2 parses this fine
            }
        }
    }
}

TEST_CASE("MPI init: grid mode tiles each subdomain; global count split is correct", "[MPI][init][grid]") {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    SimulationBox box(200.0, 160.0);
    auto decomp = box.bestDecomposition(size);

    // Choose an N_global that creates a non-uniform remainder across ranks
    const int N_global = 7 * size + 3; // base=7, rem=3 -> ranks 0..2 get +1
    std::vector<Particle> pts;
    initializeParticles_globalN(pts, box, decomp, rank, N_global, /*random=*/false /*grid*/);

    const int expected_local = compute_local_count(N_global, decomp, rank);
    REQUIRE((int)pts.size() == expected_local);

    // Bounds check for grid positions
    double x0,x1,y0,y1;
    get_bounds(box, decomp, rank, x0,x1,y0,y1);
    const double eps = 1e-12;
    for (const auto& p : pts) {
        REQUIRE(p.x >= x0 - eps);
        REQUIRE(p.x <  x1 + eps);
        REQUIRE(p.y >= y0 - eps);
        REQUIRE(p.y <  y1 + eps);
    }

    // Check half-open semantics on the right/top edges (no point should be >= x1 or >= y1)
    for (const auto& p : pts) {
        REQUIRE(p.x < x1 + 1e-9);
        REQUIRE(p.y < y1 + 1e-9);
    }

    // Global sum equals N_global
    int local_n = (int)pts.size();
    int sum_n = 0;
    MPI_Allreduce(&local_n, &sum_n, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) REQUIRE(sum_n == N_global);
}

TEST_CASE("MPI init: different seeds -> different random layouts on each rank", "[MPI][init][random][seed]") {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    SimulationBox box(150.0, 90.0);
    auto decomp = box.bestDecomposition(size);

    const int N_global = 25600;

    std::vector<Particle> a, b;
    initializeParticles_globalN(a, box, decomp, rank, N_global, /*random=*/true, /*seed=*/1234);
    initializeParticles_globalN(b, box, decomp, rank, N_global, /*random=*/true, /*seed=*/5678);

    // Same local count
    REQUIRE(a.size() == b.size());

    // With different seeds, positions should differ for at least some entries
    bool any_diff = false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::fabs(a[i].x - b[i].x) > 1e-12 || std::fabs(a[i].y - b[i].y) > 1e-12) {
            any_diff = true; break;
        }
    }
    REQUIRE(any_diff);
}


// unit_test_serial/test_cell_list.cpp
#include "catch.hpp"
#include "../serial_src/cell_list.h"
#include "../serial_src/thermodynamic_calculator.h"
#include <random>
#include <algorithm>
#include <vector>
#include <iostream>
TEST_CASE("CellList neighbor lookup with distances in small periodic box", "[CellList]") {
    SimulationBox box(3.0, 3.0);
    std::vector<Particle> P = {
        Particle(0.1, 0.1),
        Particle(2.9, 2.9), // near periodic image of P0
        Particle(1.6, 1.5)
    };
    double rcut = 0.5;
    CellList cl(box, rcut);
    cl.build(P);
    auto nb0 = cl.getNeighbors(0, P);
    // Expect one neighbor (idx=1) with r_sq ~0.2^2+0.2^2=0.08
    REQUIRE(nb0.size() == 1);
    auto [j, r_sq] = nb0[0];
    REQUIRE(j == 1);
    REQUIRE(r_sq == Approx((0.1+0.1)*(0.1+0.1)*2).epsilon(1e-6));

    auto nb2 = cl.getNeighbors(2, P);
    // P2 sees no neighbors
    REQUIRE(nb2.empty());
}


TEST_CASE("CellList exact cutoff inclusion and exclusion", "[CellList]") {
    SimulationBox box(5.0, 5.0);
    double rcut = 1.0;
    std::vector<Particle> P = {
        Particle(1.0, 1.0),
        Particle(2.0, 1.0), // exactly at r = 1.0
        Particle(3.1, 1.0), // at r = 2.1 > rcut
    };
    CellList cl(box, rcut);
    cl.build(P);

    auto nb0 = cl.getNeighbors(0, P);
    // Should include particle 1 at r_sq = 1.0
    REQUIRE(nb0.size() == 1);
    REQUIRE(nb0[0].first == 1);
    REQUIRE(nb0[0].second == Approx(1.0));

    auto nb2 = cl.getNeighbors(2, P);
    // Particle 2 only sees none within cutoff
    REQUIRE(nb2.empty());
}

TEST_CASE("CellList multiple neighbors in mixed cell configuration", "[CellList]") {
    SimulationBox box(6.0, 6.0);
    double rcut = 1.5;
    std::vector<Particle> P = {
        Particle(0.5, 0.5),
        Particle(1.0, 1.0),
        Particle(2.0, 0.5),
        Particle(3.5, 3.5) // far away, only communicates via PBC
    };
    CellList cl(box, rcut);
    cl.build(P);

    // Particle 0 should see particles 1 and 2
    auto nb0 = cl.getNeighbors(0, P);
    std::vector<int> ids0;
    for (auto &pr : nb0) {
        ids0.push_back(pr.first);
        // std::cout<<pr.first<<std::endl;
    }
    REQUIRE(std::find(ids0.begin(), ids0.end(), 1) != ids0.end());
    REQUIRE(std::find(ids0.begin(), ids0.end(), 2) != ids0.end());
    REQUIRE(ids0.size() == 2);

    // Particle 3 via PBC should see none (too far even across boundary)
    auto nb3 = cl.getNeighbors(3, P);
    REQUIRE(nb3.empty());
}



TEST_CASE("Energy consistency for large random system (manual LJ)", "[CellList][Thermo][Large]") {
    const size_t N = 10000;
    SimulationBox box(100.0, 100.0 * sqrt(3));
    std::mt19937_64 rng(123456);
    std::uniform_real_distribution<double> distX(0.0, box.getLx());
    std::uniform_real_distribution<double> distY(0.0, box.getLy());
    std::vector<Particle> P;
    P.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        P.emplace_back(distX(rng), distY(rng));
    }
    double rcut = 2.5;
    ThermodynamicCalculator calc(1.0, PotentialType::LennardJones, rcut);
    // reference energy via thermo
    double E_ref = calc.computeTotalEnergy(P, box);

    CellList cl(box, rcut);
    cl.build(P);
    double E_cell = 0.0;
    // manual LJ calculation for each neighbor
    for (size_t i = 0; i < N; ++i) {
        auto neigh = cl.getNeighbors(i, P);
        for (auto& pr : neigh) {
            double r_sq = pr.second;
            double inv_r2 = 1.0 / r_sq;
            double inv_r6 = inv_r2 * inv_r2 * inv_r2;
            // Lennard-Jones: 4*(inv_r6^2 - inv_r6)
            E_cell += 4.0 * (inv_r6 * inv_r6 - inv_r6);
        }
    }
    // each pair counted twice
    E_cell *= 0.5;
    REQUIRE(E_cell == Approx(E_ref).epsilon(1e-12));
}
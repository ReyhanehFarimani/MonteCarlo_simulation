#include "../serial_src/verlet_list.h"
#include "../serial_src/initial.h"
#include "../serial_src/thermodynamic_calculator.h"
#include "../serial_src/potential.h"
#include "catch.hpp"

TEST_CASE("build includes close particles in neighbor list", "VerletList") {
    double rc = 1.0;
    double skin = 0.2;
    double cellSize = rc + skin;
    VerletList vl(rc, skin);

    // Two particles separated by 0.7 < rc
    std::vector<Particle> positions = { {0.1, 0.1}, {0.8, 0.1} };
    SimulationBox box(2.0, 2.0);

    vl.build(positions, box);
    const auto& nbrs = vl.neighbors();

    REQUIRE(nbrs.size() == 2);
    // each should list the other
    REQUIRE(nbrs[0].size() == 1);
    REQUIRE(nbrs[0][0] == 1);
    REQUIRE(nbrs[1].size() == 1);
    REQUIRE(nbrs[1][0] == 0);
}


TEST_CASE("build excludes distant particles", "VerletList") {
    double rc = 1.0;
    double skin = 0.1;
    VerletList vl(rc, skin);

    // Two particles separated by 1.5 > rc+skin
    std::vector<Particle> positions = { {0.0, 0.0}, {1.5, 0.0} };
    SimulationBox box(3.0, 3.0);

    vl.build(positions, box);
    const auto& nbrs = vl.neighbors();

    REQUIRE(nbrs.size() == 2);
    REQUIRE(nbrs[0].empty());
    REQUIRE(nbrs[1].empty());
}


TEST_CASE("needsRebuild triggers after displacement > skin/2", "VerletList") {
    double rc = 1.0;
    double skin = 0.5;
    VerletList vl(rc, skin);

    std::vector<Particle> positions = { {0.0, 0.0} };
    SimulationBox box(10.0, 10.0);

    vl.build(positions, box);
    REQUIRE(!vl.needsRebuild());

    // Move particle by 0.3 > skin/2 = 0.25
    positions[0] = {0.3, 0.0};
    vl.computeMaxDisplacement(positions, box);
    REQUIRE(vl.needsRebuild());
}


TEST_CASE("getNeighborsForTest returns correct neighbors for insertion", "VerletList") {
    double rc = 1.0;
    double skin = 0.2;
    VerletList vl(rc, skin);

    std::vector<Particle> positions = { {0.1, 0.1}, {0.9, 0.1}, {2.0, 6.5} };
    SimulationBox box(4.0, 6.0);

    vl.build(positions, box);

    Particle p_new{0.5, 0.1};  // should be neighbors with particles 0 and 1
    std::vector<std::size_t> out;
    vl.getNeighborsForTest(p_new, positions, box, out);

    // Direct build uses rc only, so both within rc
    REQUIRE(out.size() == 2);
    // check presence of both indices without using || in REQUIRE
    bool has0 = false, has1 = false;
    for (auto idx : out) {
        if (idx == 0) has0 = true;
        if (idx == 1) has1 = true;
    }
    REQUIRE(has0);
    REQUIRE(has1);
}

TEST_CASE("Total energy consistency: VerletList vs ThermodynamicCalculator on random 100-particle system", "VerletList[Thermo]") {
    // Parameters
    double rc = 2.5;
    double skin = 0.2;
    VerletList vl(rc, skin);
    ThermodynamicCalculator calc(1.0, PotentialType::LennardJones, rc);

    // Generate 100 random particles in a 10x5 box with fixed seed
    const size_t N = 10000;
    SimulationBox box(10.0, 50.0);
    std::vector<Particle> positions;
    positions.reserve(N);
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> distX(0.0, box.getLx());
    std::uniform_real_distribution<double> distY(0.0, box.getLy());
    for (size_t i = 0; i < N; ++i) {
        positions.emplace_back(distX(rng), distY(rng));
    }

    // Build Verlet list for extended cutoff (rc + skin)
    vl.build(positions, box);
    const auto& nbrs = vl.neighbors();

    // Compute potential energy via VerletList pairs
    double energyFromList = 0.0;
    for (size_t i = 0; i < N; ++i) {
        for (auto j : nbrs[i]) {
            if (j > i) {
                double d2 = box.minimumImageDistanceSquared(positions[i], positions[j]);
                double u = 0.0;
                if (d2<rc*rc)
                    u = computePairPotential(d2, PotentialType::LennardJones, 0, 0, 0, 0);
                energyFromList += u;
            }
        }
    }

    // Compute total energy via ThermodynamicCalculator
    double energyCalc = calc.computeTotalEnergy(positions, box);

    // Compare energies
    REQUIRE(energyFromList == Approx(energyCalc).margin(1e-8));
}

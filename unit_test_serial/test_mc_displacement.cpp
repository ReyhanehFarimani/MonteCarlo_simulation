#include "catch.hpp"
#include "../serial_src/initial.h"
#include "../serial_src/thermodynamic_calculator.h"
#include "../serial_src/rng.h"
#include "../serial_src/logging.h"
#include "../serial_src/MC.h"

TEST_CASE("Single displacement consistency over many steps: cell list vs brute force", "[MC]") {
    const unsigned SEED = 4922;
    const double BOX_L = 10.0;
    const size_t N = 100;
    const double RCUT = 2.5;
    const double TEMP = 1.0;

    // Three RNGs with same seed: warm-up, cell-list driver, brute-force driver
    RNG rngWarm(SEED);
    RNG rngCell(SEED);
    RNG rngBF(SEED);

    // Simulation box
    SimulationBox box(BOX_L, BOX_L);

    // Initialize particle set and warm up via cell-list moves
    std::vector<Particle> pWarm;
    initializeParticles(pWarm, box, static_cast<int>(N), true, SEED);
    ThermodynamicCalculator calc(TEMP, 0.0,  PotentialType::LennardJones, RCUT, 0.0);
    Logging dummyLogger("tmp1", "tmp2");
    MonteCarlo mcWarm(box,  pWarm, calc, RCUT, 0.1, 0.1, Ensemble::NVT, dummyLogger, rngWarm);
    // Warm up: N*1000 displacement moves
    for (int t = 0; t < static_cast<int>(N) * 1000; ++t) {
        mcWarm.stepCellList();
        if (t % 10 == 0) mcWarm.updateCellList();
    }

    // Copy warmed configuration for both methods
    std::vector<Particle> p1 = pWarm;
    std::vector<Particle> p2 = pWarm;

    // Create two MC drivers on the warmed state
    MonteCarlo mcCell(box, p1, calc, RCUT,0.1, 0.1, Ensemble::NVT, dummyLogger, rngCell);
    MonteCarlo mcBF(box, p2, calc, RCUT, 0.1, 0.1, Ensemble::NVT, dummyLogger, rngBF);

    // Check initial energies
    double E0_cell = mcCell.getEnergy();
    double E0_bf   = mcBF.getEnergy();
    REQUIRE(E0_cell == Approx(E0_bf));

    // Perform 10,000 displacement steps and compare
    for (int t = 0; t < 10000; ++t) {
        bool accCell = mcCell.stepCellList();
        bool accBF   = mcBF.stepBruteForce();
        // Periodically ensure cell-list is rebuilt exactly
        if (t % 100 == 0) mcCell.updateCellList();
        
        // Energies should remain close
        double Ec = mcCell.getEnergy();
        double Eb = mcBF.getEnergy();
        REQUIRE(std::abs(Ec - Eb) / std::max(1.0, std::abs(Ec)) < 5e-4);
        // Accept/reject decisions must match
        REQUIRE(accCell == accBF);
    }

    // Final positions match exactly
    auto const &rCell = mcCell.getParticles();
    auto const &rBF   = mcBF.getParticles();
    REQUIRE(rCell.size() == rBF.size());
    for (size_t i = 0; i < rCell.size(); ++i) {
        REQUIRE(rCell[i].x == Approx(rBF[i].x));
        REQUIRE(rCell[i].y == Approx(rBF[i].y));
    }
}
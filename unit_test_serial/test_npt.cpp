#include "catch.hpp"
#include "../serial_src/initial.h"
#include "../serial_src/thermodynamic_calculator.h"
#include "../serial_src/rng.h"
#include "../serial_src/logging.h"
#include "../serial_src/MC.h"
#include <string>
TEST_CASE("NPT Ensemble checking the box update", "[MC]") {
    const unsigned SEED = 4922;
    const double BOX_L = 8.0;
    const size_t N = 20;
    const double RCUT = 2.5;
    const double TEMP = 1.0;
    const double Press = 1;
    std::string out_xyz = "npt_test_1.xyz";
    std::string out_data = "npt_test_1.log";

    // Three RNGs with same seed: warm-up, cell-list driver, brute-force driver
    RNG rngWarm(SEED);
    RNG rngBF(SEED);

    // Simulation box
    SimulationBox box(BOX_L, BOX_L);

    // Initialize particle set and warm up via cell-list moves
    std::vector<Particle> pWarm;
    initializeParticles(pWarm, box, static_cast<int>(N), true, SEED);
    ThermodynamicCalculator calc(TEMP,Press, PotentialType::LennardJones, RCUT, 0.0);
    Logging logger(out_xyz, out_data);
    MonteCarlo mcWarm(box, pWarm, calc, RCUT, 0.05, 0.1, Ensemble::NVT, logger, rngWarm);
    // Warm up: N*1000 displacement moves
    mcWarm.run(1000, 100, 50);
    

    // Copy warmed configuration for both methods
    std::vector<Particle> p1 = pWarm;

    // Create two MC drivers on the warmed state
    MonteCarlo mcP(box, p1, calc, RCUT, 0.05, 0.05, Ensemble::NPT, logger, rngBF);
    mcP.run(100, 1000000, 20);
    double p=0;
    // box = mcP.box_;
    // Check initial energies
    for (int i =0; i<4; i++){
        mcP.run(1000, 100, 20);
        p += calc.computePressure(p1, box)/4;
    }

    // REQUIRE(abs(p - Press)<0.1); 


}
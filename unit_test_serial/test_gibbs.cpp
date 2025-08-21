#include "catch.hpp"
#include "../serial_src/initial.h"
#include "../serial_src/thermodynamic_calculator.h"
#include "../serial_src/rng.h"
#include "../serial_src/logging.h"
#include "../serial_src/GibbsMC.h"
#include <string>
TEST_CASE("Gibbs Ensmble testing the simple two box scenario", "[MC]") {
    const unsigned SEED = 4922;
    const double BOX_L_1 = 8.0;
    const size_t N_1 = 20;
    const double BOX_L_2 = 10.0;
    const size_t N_2 = 200;
    const double RCUT = 2.5;
    const double TEMP = 1.0;
    
    std::string out_xyz_1 = "gibbs_test_1.xyz";
    std::string out_data_1 = "gibbs_test_1.log";
    std::string out_xyz_2 = "gibbs_test_2.xyz";
    std::string out_data_2 = "gibbs_test_2.log";
    RNG rngWarm(SEED);
    SimulationBox box_1(BOX_L_1, BOX_L_1);
    SimulationBox box_2(BOX_L_2, BOX_L_2);

    // Initialize particle set and warm up via cell-list moves
    std::vector<Particle> p1;
    initializeParticles(p1, box_1, static_cast<int>(N_1), true, SEED);

    std::vector<Particle> p2;
    initializeParticles(p2, box_2, static_cast<int>(N_2), true, SEED);

    ThermodynamicCalculator calc_1(TEMP,0, PotentialType::LennardJones, RCUT, 0.0);
    Logging logger_1(out_xyz_1, out_data_1);

    ThermodynamicCalculator calc_2(TEMP,0, PotentialType::LennardJones, RCUT, 0.0);
    Logging logger_2(out_xyz_2, out_data_2);

    GibbsMonteCarlo gmc(box_1, box_2, p1, p2, calc_1, calc_2, RCUT, 0.05, 0.0, logger_1, logger_2, rngWarm);
   
    gmc.run(10000, 100, 10);

}
#include "catch.hpp"
#include "../serial_src/initial.h"
#include "../serial_src/thermodynamic_calculator.h"
#include "../serial_src/rng.h"
#include "../serial_src/logging.h"
#include "../serial_src/GibbsMC.h"
#include <string>
#include <iostream>
TEST_CASE("Gibbs Ensmble testing the simple two box scenario", "[MC]") {
    const unsigned SEED = 4922;
    const double BOX_L_1 = 30.0;
    const size_t N_1 = 100;
    const double BOX_L_2 = 10.0;
    const size_t N_2 = 200;
    const double RCUT = 2.0;
    const double TEMP = 1.0;
    
    std::string out_xyz_1 = "gibbs_test_1.xyz";
    std::string out_data_1 = "gibbs_test_1.log";
    std::string out_xyz_2 = "gibbs_test_2.xyz";
    std::string out_data_2 = "gibbs_test_2.log";
    RNG rngWarm(SEED);
    SimulationBox box_1(BOX_L_1, BOX_L_1);
    SimulationBox box_2(BOX_L_2, BOX_L_2);

    double total_v = box_1.getV() + box_2.getV();
    double total_n = N_1 + N_2;
    // Initialize particle set and warm up via cell-list moves
    std::vector<Particle> p1;
    initializeParticles(p1, box_1, static_cast<int>(N_1), true, SEED);

    std::vector<Particle> p2;
    initializeParticles(p2, box_2, static_cast<int>(N_2), true, SEED);

    ThermodynamicCalculator calc_1(TEMP,0, PotentialType::WCA, RCUT, 0.0);
    Logging logger_1(out_xyz_1, out_data_1);

    ThermodynamicCalculator calc_2(TEMP,0, PotentialType::WCA, RCUT, 0.0);
    Logging logger_2(out_xyz_2, out_data_2);


    

    GibbsMonteCarlo gmc2(box_1, box_2, p1, p2, calc_1, calc_2, RCUT, 0.05, 0.01, logger_1, logger_2, rngWarm);
   
    gmc2.run(10000, 100, 100);

    double current_v = box_1.getV() + box_2.getV();
    std::cout<<(box_1.getV())<<std::endl; 
    double current_n = calc_1.getNumParticles(p1) + calc_2.getNumParticles(p2);

    REQUIRE(total_v == Approx(current_v));
    REQUIRE(total_n == current_n);
    

}
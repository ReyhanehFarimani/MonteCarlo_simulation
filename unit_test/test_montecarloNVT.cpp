#include "../src/initial.h"
#include "../src/input.h"
#include "../src/logging.h"
#include "../src/potential.h"
#include "../src/simulation.h"
#include <iostream>
#include "test.h"

void testMonteCarloEnergyReduction() {
    SimulationBox box(10.0, 10.0);
    Simulation sim(box, PotentialType::LennardJones, 0.1, 2, 0.1, 2.5, 1, 1, 1);
    sim.initializeParticles(true); // Random placement
    double initialEnergy = sim.getEnergy();
    // Create logger (you can also use a dummy logger if you don't want to log to files)
    Logging logger("positions_test.xyz", "simulation_data_test.dat");
    sim.run(10000000, 0, 10000, logger, SimulationType::MonteCarloNVT);

    double finalEnergy = sim.getEnergy();
    
    if (finalEnergy <= initialEnergy) {
        std::cout << "Test Passed: Final energy is less than or equal to initial energy." << std::endl;
    } else {
        std::cout << "Test Failed: Final energy is greater than initial energy." << std::endl;
    }
}
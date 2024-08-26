#include "../src/initial.h"
#include "../src/input.h"
#include "../src/logging.h"
#include "../src/potential.h"
#include "../src/simulation.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "test.h"

void ideal_gas_GCMC(){
    const double temperature = 1.0; // Temperature set to 10.0
    const double r_cutoff = 2.5; // Cutoff radius
    const double equilibration_steps = 5000000;
    const double production_steps = 1000000;
    const double mu_increment = 0.1;
    const double seed = 12345;
    double box_length_x = 20;
    double box_length_y = 20;
    SimulationBox box(box_length_x, box_length_y);
    for (double mu = -2.0; mu < 2.1; mu += mu_increment){

        Simulation sim(box, PotentialType::Ideal, SimulationType::GCMC, temperature, 100, 0.1, 10, 0 ,mu, seed, 1, 100000);
        sim.initializeParticles(0);

        // Log file names based on density
        std::string mu_str = std::to_string(mu).substr(0, 4); // Limiting to 3 decimal places
        Logging logger("position_" + mu_str + ".xyz", "data_" + mu_str + ".txt");


        // Run the simulation
        sim.run(production_steps, equilibration_steps, 10000, logger);

        std::system("python compare_theory_numerical.py");
    }

}
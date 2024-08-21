// pressure_NVT_LJ_test.cpp
#include "../src/initial.h"
#include "../src/input.h"
#include "../src/logging.h"
#include "../src/potential.h"
#include "../src/simulation.h"
#include <iostream>
#include <cmath>
#include "test.h"

void pressure_test_NVT_LJ()
{
    const int N = 32 * 32; // Reduced number of particles to 2500
    const double temperature = 10.0; // Temperature set to 10.0
    const double r_cutoff = 2.5; // Cutoff radius
    const double equilibration_steps = 5000000;
    const double production_steps = 1000000;
    const double density_increment = 0.0025;
    
    for (double density = 1.19; density < 1.23; density += density_increment) {
        double L2 = N / density;
        double box_length_x = sqrt(L2 / (sqrt(3.0) / 2.0));
        double box_length_y = sqrt(L2 * (sqrt(3.0) / 2.0));
        SimulationBox box(box_length_x, box_length_y);
        Simulation sim(box, PotentialType::LennardJones, temperature, N, 0.1, r_cutoff, 0, 1, 100);
        sim.initializeParticles(0);

        // Log file names based on density
        std::string density_str = std::to_string(density).substr(0, 4); // Limiting to 3 decimal places
        Logging logger("position_" + density_str + ".xyz", "data_" + density_str + ".txt");

        // Run the simulation
        sim.run(production_steps, equilibration_steps, 10000, logger, SimulationType::MonteCarloNVT);

        // Compute and print the pressure for this density
        double pressure = sim.calculatePressure();
        std::cout << "Density: " << density << " Pressure: " << pressure << std::endl;
    }
}

#include <iostream>
#include <cmath>
#include "../serial_src/initial.h"
#include "../serial_src/MC.h"
#include "../serial_src/thermodynamic_calculator.h"
#include "../serial_src/logging.h"
#include "../serial_src/rng.h"

int main() {
    // --- Parameters ---
    int N = 100;                    // number of particles
    double L= 40;                // box side length
    double temperature = 0.450;
    double rcut = 2.5;
    double z = 0.0242;
    int steps = 1000000;
    int fOut = 1000;
    int fCell = 1000;
    // --- Setup ---
    SimulationBox box(L * sqrt(sqrt(3)/2), L/ sqrt(sqrt(3)/2));
    std::vector<Particle> particles;
    initializeParticles(particles, box, N, true, 42); // seed = 42

    // --- Thermodynamics and logger ---
    ThermodynamicCalculator calc(temperature, PotentialType::LennardJones, rcut, z);
    Logging logger("positions_T0.45_z_0.0242_L40.xyz", "data_T0.45_z_0.0242_L40.txt");
    RNG rng(123); // reproducible RNG

    // --- Monte Carlo driver ---
    MonteCarlo mc(box, particles, calc, rcut, Ensemble::GCMC, logger, rng);
    // --- Equillibrate ---
    std::cout<<"Equilibration:"<<std::endl; 
    mc.run(100000, 10000000, 1000);  
    std::cout << "Equillibrated energy: " << mc.getEnergy() << "\n";
    // --- Run ---
    std::cout<<"Running:"<<std::endl;
    
    mc.run(steps, fOut, fCell);

    // --- Output final energy ---
    std::cout << "Final energy: " << mc.getEnergy() << "\n";
    logger.close();
    return 0;
}

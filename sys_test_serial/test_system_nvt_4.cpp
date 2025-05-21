#include <iostream>
#include <cmath>
#include "../serial_src/initial.h"
#include "../serial_src/MC.h"
#include "../serial_src/thermodynamic_calculator.h"
#include "../serial_src/logging.h"
#include "../serial_src/rng.h"

int main() {
    // --- Parameters ---
    int N = 32 * 32;                    // number of particles
    double rho = 1.215;
    double L= sqrt(N/rho);                // box side length
    double temperature = 10.0;
    double rcut = 2.5;
    int steps = 1000000;
    int fOut = 100000;
    int fCell = 100;
    // --- Setup ---
    SimulationBox box(L * sqrt(sqrt(3)/2), L/ sqrt(sqrt(3)/2));
    std::vector<Particle> particles;
    initializeParticles(particles, box, N, true, 90); // seed = 92

    // --- Thermodynamics and logger ---
    ThermodynamicCalculator calc(temperature, PotentialType::LennardJones, rcut);
    Logging logger("positions_T10.0_rho_1.215_N32_2.xyz", "data_T10.0_rho_1.215_N32_2.txt");
    RNG rng(129); // reproducible RNG

    // --- Monte Carlo driver ---
    MonteCarlo mc(box, particles, calc, rcut, Ensemble::NVT, logger, rng);
    // --- Equillibrate ---
    std::cout<<"Equilibration:"<<std::endl;
    mc.run(1000000, 1000000000, 10);    
    std::cout << "Equillibrated energy: " << mc.getEnergy() << "\n";
    // --- Run ---
    std::cout<<"Running:"<<std::endl;
    
    mc.run(steps, fOut, fCell);

    // --- Output final energy ---
    std::cout << "Final energy: " << mc.getEnergy() << "\n";
    logger.close();
    return 0;
}

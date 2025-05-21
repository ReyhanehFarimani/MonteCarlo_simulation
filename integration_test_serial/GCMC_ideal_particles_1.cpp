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
    double L= 10;                // box side length
    double temperature = 1.0;
    double rcut = 2.5;
    int steps = 100000;
    int fOut = 10000;
    int fCell = 10000;
    double z = 0.1;
    // --- Setup ---
    SimulationBox box(L * sqrt(sqrt(3)/2), L/ sqrt(sqrt(3)/2));
    std::vector<Particle> particles;
    initializeParticles(particles, box, N, true, 42); // seed = 42

    // --- Thermodynamics and logger ---
    ThermodynamicCalculator calc(temperature, PotentialType::Ideal, rcut, z=z);
    Logging logger("positions_z_0.1_V_100_T1.0.xyz", "data__z_0.1_V_100_T1.0.txt");
    RNG rng(123); // reproducible RNG

    // --- Monte Carlo driver ---
    MonteCarlo mc(box, particles, calc, rcut, Ensemble::GCMC, logger, rng);

    // --- Run ---
    std::cout<<"z:"<<z<<std::endl;
    mc.run(steps, fOut, fCell);

    logger.close();

    z = 1.0;
    // --- Setup ---
    initializeParticles(particles, box, N, true, 420); // seed = 42

    // --- Thermodynamics and logger ---
    ThermodynamicCalculator calc1(temperature, PotentialType::Ideal, rcut, z=z);
    Logging logger1("positions_z_1.0_V_100_T1.0.xyz", "data__z_1.0_V_100_T1.0.txt");

    // --- Monte Carlo driver ---
    MonteCarlo mc1(box, particles, calc1, rcut, Ensemble::GCMC, logger1, rng);

    // --- Run ---
    std::cout<<"z:"<<z<<std::endl;
    mc1.run(steps, fOut, fCell);

    logger1.close();

    z = 10.0;
    // --- Setup ---
    initializeParticles(particles, box, N, true, 42); // seed = 42

    // --- Thermodynamics and logger ---
    ThermodynamicCalculator calc2(temperature, PotentialType::Ideal, rcut, z=z);
    Logging logger2("positions_z_10.0_V_100_T1.0.xyz", "data__z_10.0_V_100_T1.0.txt");

    // --- Monte Carlo driver ---
    MonteCarlo mc2(box, particles, calc2, rcut, Ensemble::GCMC, logger2, rng);

    // --- Run ---
    std::cout<<"z:"<<z<<std::endl;
    mc2.run(steps, fOut, fCell);

    logger2.close();

    z = 0.5;
    // --- Setup ---
    initializeParticles(particles, box, N, true, 42); // seed = 42

    // --- Thermodynamics and logger ---
    ThermodynamicCalculator calc3(temperature, PotentialType::Ideal, rcut, z=z);
    Logging logger3("positions_z_0.5_V_100_T1.0.xyz", "data__z_0.5_V_100_T1.0.txt");

    // --- Monte Carlo driver ---
    MonteCarlo mc3(box, particles, calc3, rcut, Ensemble::GCMC, logger3, rng);

    // --- Run ---
    std::cout<<"z:"<<z<<std::endl;
    mc3.run(steps, fOut, fCell);

    logger3.close();

    z = 5.0;
    // --- Setup ---
    initializeParticles(particles, box, N, true, 42); // seed = 42

    // --- Thermodynamics and logger ---
    ThermodynamicCalculator calc4(temperature, PotentialType::Ideal, rcut, z=z);
    Logging logger4("positions_z_5.0_V_100_T1.0.xyz", "data__z_5.0_V_100_T1.0.txt");

    // --- Monte Carlo driver ---
    MonteCarlo mc4(box, particles, calc4, rcut, Ensemble::GCMC, logger4, rng);

    // --- Run ---
    std::cout<<"z:"<<z<<std::endl;
    mc4.run(steps, fOut, fCell);

    logger4.close();


    return 0;
}

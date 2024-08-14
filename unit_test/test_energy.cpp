// test_energy.cpp
#include "../src/initial.h"
#include "../src/input.h"
#include "../src/logging.h"
#include "../src/potential.h"
#include "../src/simulation.h"
#include <iostream>
#include "test.h"

void testBoundaryEnergy() {
    SimulationBox box(10.0, 10.0);
    Simulation sim(box, PotentialType::LennardJones, 1.0, 2, 0.01, 2.5);

    std::vector<Particle> particles = {
        Particle(0.1, 0.1),
        Particle(9.9, 9.9)
    };
    sim.initializeParticles(false);

    double expectedEnergy = lennardJonesPotential(box.minimumImageDistanceSquared(particles[0], particles[1]));
    double calculatedEnergy = sim.getEnergy();

    if (fabs(expectedEnergy - calculatedEnergy) < 1e-5) {
        std::cout << "Boundary Energy Test Passed." << std::endl;
    } else {
        std::cout << "Boundary Energy Test Failed." << std::endl;
    }
}


void testEnergyLJ() {
    SimulationBox box(10.0, 10.0);
    Simulation sim(box, PotentialType::LennardJones, 1.0, 2, 0.01, 2.5);

    std::vector<Particle> particles = {
        Particle(1.0, 1.0),
        Particle(2.0, 2.0)
    };
    sim.initializeParticles(false);

    double expectedEnergy = lennardJonesPotential(box.minimumImageDistanceSquared(particles[0], particles[1]));
    double calculatedEnergy = sim.getEnergy();

    if (fabs(expectedEnergy - calculatedEnergy) < 1e-5) {
        std::cout << "LJ Energy Test Passed." << std::endl;
    } else {
        std::cout << "LJ Energy Test Failed." << std::endl;
    }
}
void testEnergyYukawa() {
    SimulationBox box(10.0, 10.0);
    Simulation sim(box, PotentialType::Yukawa, 1.0, 2, 0.01, 2.5);

    std::vector<Particle> particles = {
        Particle(1.0, 1.0),
        Particle(2.0, 2.0)
    };
    sim.initializeParticles(false);

    double expectedEnergy = yukawaPotential(box.minimumImageDistanceSquared(particles[0], particles[1]));
    double calculatedEnergy = sim.getEnergy();

    if (fabs(expectedEnergy - calculatedEnergy) < 1e-5) {
        std::cout << "Yukawa Energy Test Passed." << std::endl;
    } else {
        std::cout << "Yukawa Energy Test Failed." << std::endl;
    }
}

void testEnergyWCA() {
    SimulationBox box(10.0, 10.0); // Define a box of size 10x10
    Simulation sim(box, PotentialType::WCA, 1.0, 2, 0.01, 2.5);
    
    std::vector<Particle> particles = {
        Particle(1.0, 1.0),
        Particle(2.0, 2.0)
    };
    sim.initializeParticles(false); // Initialize particles at defined positions

    double expectedEnergy = wcaPotential(box.minimumImageDistanceSquared(particles[0], particles[1]));
    double calculatedEnergy = sim.getEnergy();

    if (fabs(expectedEnergy - calculatedEnergy) < 1e-5) {
        std::cout << "WCA Energy Test Passed." << std::endl;
    } else {
        std::cout << "WCA Energy Test Failed." << std::endl;
    }
}
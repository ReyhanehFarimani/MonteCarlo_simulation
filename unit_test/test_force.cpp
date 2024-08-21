// test_energy.cpp
#include "../src/initial.h"
#include "../src/input.h"
#include "../src/logging.h"
#include "../src/potential.h"
#include "../src/simulation.h"
#include <iostream>
#include "test.h"

void testBoundaryForce() {
    SimulationBox box(10.0, 10.0);
    Simulation sim(box, PotentialType::LennardJones, 1.0, 2, 0.1, 2.5, 0, 0, 0);

    std::vector<Particle> particles = {
        Particle(9.0, 9.0),
        Particle(0.0, 0.0)
    };
    sim.initializeParticles(false);
    sim.setParticlePosition(0, particles[0].x, particles[0].y);
    sim.setParticlePosition(1, particles[1].x, particles[1].y);

    double expectedForce = -6.0 * 3.0 / 8.0;
    double calculatedForce = sim.computeTotalForce();
    if (fabs(expectedForce - calculatedForce) < 1e-5) {
        std::cout << "Boundary Force Test Passed." << std::endl;
    } else {
        std::cout << "Boundary Force Test Failed." << std::endl;
    }
}


// void testEnergyYukawa() {
//     SimulationBox box(10.0, 10.0);
//     Simulation sim(box, PotentialType::Yukawa, 1.0, 2, 0.1, 2.5, 0, 0, 0);

//     std::vector<Particle> particles = {
//         Particle(1.0, 1.0),
//         Particle(2.0, 2.0)
//     };
//     sim.initializeParticles(false);
//     sim.setParticlePosition(0, particles[0].x, particles[0].y);
//     sim.setParticlePosition(1, particles[1].x, particles[1].y);

//     double expectedEnergy = yukawaPotential(box.minimumImageDistanceSquared(particles[0], particles[1]));
//     double calculatedEnergy = sim.getEnergy();

//     if (fabs(expectedEnergy - calculatedEnergy) < 1e-5) {
//         std::cout << "Yukawa Energy Test Passed." << std::endl;
//     } else {
//         std::cout << "Yukawa Energy Test Failed." << std::endl;
//     }
// }

void testForceWCA() {
    SimulationBox box(10.0, 10.0); // Define a box of size 10x10
    Simulation sim(box, PotentialType::WCA, 1.0, 2, 0.1, 2.5, 0, 0, 0);
    
    std::vector<Particle> particles = {
        Particle(1.0, 1.0),
        Particle(1.0, 1.25)
    };
    sim.initializeParticles(false); // Initialize particles at defined positions
    sim.setParticlePosition(0, particles[0].x, particles[0].y);
    sim.setParticlePosition(1, particles[1].x, particles[1].y);

    double expectedForce = 48 * (pow(4, 12) - pow(4, 6)/2);
    double calculatedForce = sim.computeTotalForce();
    if (fabs(expectedForce - calculatedForce) < 1e-5) {
        std::cout << "WCA Force Test Passed." << std::endl;
    } else {
        std::cout << "WCA Force Test Failed." << std::endl;
    }
}
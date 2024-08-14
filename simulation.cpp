#include "simulation.h"

/**
 * @brief Constructs a Simulation object with the specified parameters.
 * 
 * @param box The simulation box.
 * @param potentialType The type of potential used in the simulation.
 * @param temperature The temperature of the simulation.
 * @param numParticles The number of particles in the simulation.
 * @param timeStep The time step for the simulation.
 * @param r2cut The squared distance cutoff for potential calculations.
 */
Simulation::Simulation(const SimulationBox &box, PotentialType potentialType, double temperature, int numParticles, double timeStep, double r2cut)
    : box(box), potentialType(potentialType), temperature(temperature), numParticles(numParticles), timeStep(timeStep), r2cut(r2cut) {
    particles.resize(numParticles);
}

/**
 * @brief Initializes particles in the simulation box.
 * 
 * @param randomPlacement If true, particles are placed randomly; otherwise, in a grid.
 */
void Simulation::initializeParticles(bool randomPlacement) {
    ::initializeParticles(particles, box, numParticles, randomPlacement);
}

/**
 * @brief Calculates the total energy of the system.
 * 
 * @return The total energy of the system.
 */
double Simulation::calculateEnergy() const {
    double energy = 0.0;
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            double r2 = box.minimumImageDistanceSquared(particles[i], particles[j]);
            if (r2 < r2cut) {
                switch (potentialType) {
                    case PotentialType::LennardJones:
                        energy += lennardJonesPotential(r2);
                        break;
                    case PotentialType::WCA:
                        energy += wcaPotential(r2);
                        break;
                    case PotentialType::Yukawa:
                        energy += yukawaPotential(r2);
                        break;
                }
            }
        }
    }
    return energy;
}

/**
 * @brief Runs the simulation for a specified number of steps.
 * 
 * @param numSteps The number of steps to run the simulation.
 * @param equilibrationTime The number of steps for equilibration before logging.
 * @param outputFrequency How often to log the results.
 * @param logger The logging object for output.
 */
void Simulation::run(int numSteps, int equilibrationTime, int outputFrequency, Logging &logger) {
    for (int timestep = 0; timestep < numSteps; ++timestep) {
        // Perform necessary simulation steps (e.g., position update, force calculation)

        // Skip logging during equilibration
        if (timestep >= equilibrationTime) {
            if (timestep % outputFrequency == 0) {
                logger.logPositions_xyz(particles);
            }
        }
    }
}

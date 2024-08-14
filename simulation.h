#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include "initial.h"
#include "potential.h"
#include "logging.h"

/**
 * @brief Manages the entire simulation process, including initialization,
 * running the simulation, and logging results.
 */
class Simulation {
private:
    SimulationBox box;                   ///< The simulation box
    std::vector<Particle> particles;     ///< Particles in the simulation
    PotentialType potentialType;         ///< The type of potential used in the simulation
    double temperature;                  ///< Temperature of the simulation
    int numParticles;                    ///< Number of particles in the simulation
    double timeStep;                     ///< Time step for the simulation
    double r2cut;                        ///< Squared distance cutoff for potential calculations

public:
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
    Simulation(const SimulationBox &box, PotentialType potentialType, double temperature, int numParticles, double timeStep, double r2cut);

    /**
     * @brief Initializes particles in the simulation box.
     * 
     * @param randomPlacement If true, particles are placed randomly; otherwise, in a grid.
     */
    void initializeParticles(bool randomPlacement);

    /**
     * @brief Calculates the total energy of the system.
     * 
     * @return The total energy of the system.
     */
    double calculateEnergy() const;

    /**
     * @brief Runs the simulation for a specified number of steps.
     * 
     * @param numSteps The number of steps to run the simulation.
     * @param equilibrationTime The number of steps for equilibration before logging.
     * @param outputFrequency How often to log the results.
     * @param logger The logging object for output.
     */
    void run(int numSteps, int equilibrationTime, int outputFrequency, Logging &logger);
};

#endif // SIMULATION_H

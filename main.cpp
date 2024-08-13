#include <iostream>
#include "initial.h"
#include "logging.h"

int main() {
    // Define the dimensions of the simulation box
    double boxLengthX = 10.0;
    double boxLengthY = 10.0;

    // Create a SimulationBox object
    SimulationBox simBox(boxLengthX, boxLengthY);

    // Number of particles
    int N = 100;

    // Create a vector to hold the particles
    std::vector<Particle> particles;

    // Initialize particles (randomly or in a grid)
    bool randomPlacement = true;  // Set to false for grid placement
    initializeParticles(particles, simBox, N, randomPlacement);

    // Create a Logging object
    Logging logger("particle_positions.xyz", "simulation_data.dat");

    // Apply periodic boundary conditions and log positions
    for (int timestep = 0; timestep < 100; ++timestep) {
        // (Simulation steps would go here)

        // Apply PBC to all particles
        for (auto &particle : particles) {
            simBox.applyPBC(particle);
        }

        // Log positions in XYZ format
        logger.logPositions_xyz(particles);

        // Optionally, log in dump format as well
        // logger.logPositions_dump(particles);
    }

    // Close the log files when done
    logger.close();

    return 0;
}

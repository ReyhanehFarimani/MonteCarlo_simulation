#include <iostream>
#include "initial.h"
#include "logging.h"
#include "input.h"
#include "potential.h"
int main(int argc, char *argv[]) {
    // Parse input from file and command line
    Input input("input.txt");
    input.parseCommandLineArgs(argc, argv);

    // Retrieve constants from the input
    double boxLengthX = input.getConstant("boxLengthX");
    double boxLengthY = input.getConstant("boxLengthY");
    int numParticles = static_cast<int>(input.getConstant("numParticles"));
    bool randomPlacement = static_cast<bool>(input.getConstant("randomPlacement"));
    double timeStep = input.getConstant("timeStep");
    int outputFrequency = static_cast<int>(input.getConstant("outputFrequency"));
    int equilibrationTime = static_cast<int>(input.getConstant("equilibrationTime"));
    
    // Retrieve the potential type from the input
    std::string potentialName = input.getFilename("potentialType");
    PotentialType potentialType = selectPotentialType(potentialName);
    auto potentialFunction = getPotentialFunction(potentialType);

    // Retrieve filenames from the input
    std::string positionFile = input.getFilename("positionFile");
    std::string dataFile = input.getFilename("dataFile");
    // Create a SimulationBox object
    SimulationBox simBox(boxLengthX, boxLengthY);

    // Create a vector to hold the particles
    std::vector<Particle> particles;

    // Initialize particles (randomly or in a grid)
    initializeParticles(particles, simBox, numParticles, randomPlacement);

    // Create a Logging object
    Logging logger("particle_positions.xyz", "simulation_data.dat");

    // Simulation loop
    for (int timestep = 0; timestep < simulationTime + equilibrationTime; ++timestep) {
        // (Simulation steps would go here, updating particle positions etc.)

        // Apply PBC to all particles
        for (auto &particle : particles) {
            simBox.applyPBC(particle);
        }

        // Skip logging during equilibration
        if (timestep < equilibrationTime) continue;

        // Log positions based on output frequency
        else if (timestep % outputFrequency == 0) {
            logger.logPositions_xyz(particles);
        }
    }

    // Close the log files when done
    logger.close();

    return 0;
}

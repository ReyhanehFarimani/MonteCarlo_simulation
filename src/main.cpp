#include <iostream>
#include "input.h"
#include "logging.h"
#include "potential.h"
#include "simulation.h"
#include "initial.h"

int main(int argc, char *argv[]) {
    // Parse input from file and command line
    Input input("input.txt");

    // Retrieve constants from the input
    double boxLengthX = input.getConstant("boxLengthX");
    double boxLengthY = input.getConstant("boxLengthY");
    int numParticles = static_cast<int>(input.getConstant("numParticles"));
    bool randomPlacement = static_cast<bool>(input.getConstant("randomPlacement"));
    double temperature = input.getConstant("temperature");
    int numSteps = static_cast<int>(input.getConstant("numSteps"));
    int outputFrequency = static_cast<int>(input.getConstant("outputFrequency"));
    int equilibrationTime = static_cast<int>(input.getConstant("equilibrationTime"));
    double timeStep = input.getConstant("timeStep");
    double r2cut = input.getConstant("r2cut");
    unsigned int seed = input.getConstant("seed");
    const bool useCellList = input.getConstant("useCellList");
    const int cellListUpdateFrequency = input.getConstant("cellListUpdateFrequency");


    // Retrieve the potential type from the input
    std::string potentialName = input.getFilename("potentialType");
    PotentialType potentialType = selectPotentialType(potentialName);

    // Retrieve the simulation type from the input
    std::string simName = input.getFilename("simulationType");
    SimulationType simType = selectSimulationType(simName);

    // Create a SimulationBox object
    SimulationBox simBox(boxLengthX, boxLengthY);

    // Create a Logging object
    std::string positionFile = input.getFilename("positionFile");
    std::string dataFile = input.getFilename("dataFile");
    Logging logger(positionFile, dataFile);

    // Create and run the simulation
    Simulation simulation(simBox, potentialType, temperature, numParticles, timeStep, r2cut, seed, useCellList, cellListUpdateFrequency);
    simulation.initializeParticles(randomPlacement);
    simulation.run(numSteps, equilibrationTime, outputFrequency, logger, simType);

    return 0;
}

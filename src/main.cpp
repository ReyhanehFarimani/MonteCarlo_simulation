#include <mpi.h>
#include <iostream>
#include "input.h"
#include "logging.h"
#include "potential.h"
#include "simulation.h"
#include "initial.h"

int main(int argc, char *argv[]) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    // Getting number of preccesors and their rank:
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    double boxLengthX, boxLengthY;
    int numParticles,numSteps, outputFrequency, equilibrationTime, cellListUpdateFrequency;
    bool randomPlacement;
    double temperature;
    double timeStep, r2cut;
    unsigned int seed;
    bool useCellList;
    std::string potentialName;
    PotentialType potentialType;
    float f_prime = 0.0;
    float f_d_prime = 0.0;
    float kappa = 0.0;
    float A_0 = 0.0;
    std::string simName;
    SimulationType simType;
    double mu = 0.0;
    std::string startFile, positionFile, dataFile;
    

    int potentialTypeInt, simTypeInt;  // To hold enum values as int for MPI_Bcast

    if (world_rank == 0){
        // Parse input from file and command line
        Input input("input.txt");

        // Retrieve constants from the input
        boxLengthX = input.getConstant("boxLengthX");
        boxLengthY = input.getConstant("boxLengthY");
        numParticles = static_cast<int>(input.getConstant("numParticles"));
        randomPlacement = static_cast<bool>(input.getConstant("randomPlacement"));
        temperature = input.getConstant("temperature");
        numSteps = static_cast<int>(input.getConstant("numSteps"));
        outputFrequency = static_cast<int>(input.getConstant("outputFrequency"));
        equilibrationTime = static_cast<int>(input.getConstant("equilibrationTime"));
        timeStep = input.getConstant("timeStep");
        r2cut = input.getConstant("r2cut");
        seed = input.getConstant("seed");
        useCellList = input.getConstant("useCellList");
        cellListUpdateFrequency = input.getConstant("cellListUpdateFrequency");
        // input trajectory file:
        startFile = input.getFilename("startingPositionFile");

        // Retrieve the potential type from the input
        potentialName = input.getFilename("potentialType");
        potentialType = selectPotentialType(potentialName);

        if (potentialType == PotentialType::AthermalStar){
            int f = input.getConstant("f");
            f_prime = (2.0 + 9.0 * f * f) / 24.0;
        }

        if (potentialType == PotentialType::ThermalStar){
            int f = input.getConstant("f");
            f_prime = (2.0 + 9.0 * f * f) / 24.0;
            A_0 = input.getConstant("A_0");
            kappa = input.getConstant("kappa");
            f_d_prime = A_0 * f * f / kappa; 
        }
        // Retrieve the simulation type from the input
        simName = input.getFilename("simulationType");
        simType = selectSimulationType(simName);
        mu = 0.0;
        if (simType == SimulationType::GCMC){
            mu = input.getConstant("mu");
        }
        

        // Create a Logging object
        positionFile = input.getFilename("positionFile");
        dataFile = input.getFilename("dataFile");

        potentialTypeInt = static_cast<int>(potentialType);  // Convert enum to int for the broadcasting
        simTypeInt = static_cast<int>(simType);  // Convert enum to int


    }

    MPI_Bcast(&boxLengthX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&boxLengthY, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numParticles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&randomPlacement, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&temperature, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&outputFrequency, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&equilibrationTime, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&timeStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&r2cut, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&useCellList, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cellListUpdateFrequency, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&f_prime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&f_d_prime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kappa, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mu, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    // Broadcast the enums as integers
    MPI_Bcast(&potentialTypeInt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&simTypeInt, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Convert back from int to enum
    potentialType = static_cast<PotentialType>(potentialTypeInt);
    simType = static_cast<SimulationType>(simTypeInt);


    // Create a SimulationBox object in each and every cpu!
    SimulationBox simBox(boxLengthX, boxLengthY);
    
    // Create Logging object only on rank 0
    Logging *logger = nullptr;
    if (world_rank == 0) {
        logger = new Logging(positionFile, dataFile);  // Rank 0 creates the logger
    }
    // Logging logger(positionFile, dataFile);

    // // Create the simulation
    Simulation simulation(simBox, potentialType, simType, temperature, numParticles, timeStep, r2cut, f_prime, f_d_prime, kappa, mu, seed, useCellList, cellListUpdateFrequency);
    simulation.initializeParticles(randomPlacement, startFile);
    simulation.run(numSteps, equilibrationTime, outputFrequency, logger);

    // Finalize the MPI environment
    MPI_Finalize();

    return 0;
}

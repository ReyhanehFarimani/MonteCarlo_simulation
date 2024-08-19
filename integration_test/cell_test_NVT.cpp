#include <cstdlib>  // For system() to run Python script
#include <iostream>
#include "../src/initial.h"
#include "../src/simulation.h"
#include "../src/input.h"
#include "../src/logging.h"
#include "../src/potential.h"
void cell_test_NVT() {
    // First Simulation
    Input input1("input_cell_list_NVT_1.txt");

    SimulationBox box1(10.0, 10.0);
    Simulation sim1(box1, PotentialType::LennardJones, 
                    input1.getConstant("temperature"), 
                    static_cast<int>(input1.getConstant("numParticles")), 
                    input1.getConstant("timeStep"), 
                    input1.getConstant("r2cut"), 
                    static_cast<unsigned int>(input1.getConstant("seed")),
                    1, //input1.getConstant("useCellList"),
                    1);//input1.getConstant("cellListUpdateFrequency"));

    sim1.initializeParticles(0); // Random placement
    Logging logger1("positions1.xyz", "simulation_data1.dat");

    sim1.run(static_cast<int>(input1.getConstant("numSteps")),
             static_cast<int>(input1.getConstant("equilibrationTime")),
             static_cast<int>(input1.getConstant("outputFrequency")),
             logger1,
             SimulationType::MonteCarloNVT);

    // Second Simulation (same seed)
    std::cout<<"second test running..."<<std::endl;
    Input input2("input_cell_list_NVT_2.txt");

    SimulationBox box2(10.0, 10.0);
    Simulation sim2(box2, PotentialType::LennardJones, 
                    input2.getConstant("temperature"), 
                    static_cast<int>(input2.getConstant("numParticles")), 
                    input2.getConstant("timeStep"), 
                    input2.getConstant("r2cut"), 
                    static_cast<unsigned int>(input2.getConstant("seed")),
                    0, //input2.getConstant("useCellList"),
                    0);//input2.getConstant("cellListUpdateFrequency"));

    sim2.initializeParticles(0); // Random placement
    Logging logger2("positions2.xyz", "simulation_data2.dat");

    sim2.run(static_cast<int>(input2.getConstant("numSteps")),
             static_cast<int>(input2.getConstant("equilibrationTime")),
             static_cast<int>(input2.getConstant("outputFrequency")),
             logger2,
             SimulationType::MonteCarloNVT);


    // Run the Python script to compare outputs
    int result = std::system("python compare_outputs.py");

    if (result == 0) {
        std::cout << "Test Passed: The simulation outputs are identical." << std::endl;
    } else {
        std::cerr << "Test Failed: The simulation outputs differ." << std::endl;
    }
    // std::system("rm position*");
}

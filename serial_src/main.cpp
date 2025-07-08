// main.cpp - Entry point for running the Grand Canonical Monte Carlo simulation
#include "input.h"
#include "initial.h"
#include "MC.h"
#include "logging.h"
#include "thermodynamic_calculator.h"
#include <iostream>
#include <string>

/**
 * @brief Main driver for the GCMC simulation.
 *
 * Reads parameters from input file, initializes the system,
 * sets up logging and Monte Carlo engine, and runs the simulation.
 */

/**
 * @brief Utility to get a constant from input file with default fallback to 0.
 *
 * @param input The Input object.
 * @param key   The name of the constant.
 * @return Value of the constant or 0.0 if not found.
 */
double getOrDefault(const Input& input, const std::string& key) {
    return input.getConstant(key);  // Already returns 0.0 if not found
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];

    try {
        // Parse input parameters
        Input input(inputFile);

        // Required parameters
        double Lx = input.getConstant("Lx");
        double Ly = input.getConstant("Ly");
        int N = static_cast<int>(input.getConstant("N"));
        double rcut = input.getConstant("rcut");
        double temperature = input.getConstant("T");
        double nSteps = input.getConstant("nSteps");
        double eSteps = input.getConstant("eSteps");
        double outputFreq = input.getConstant("outputFreq");
        double cellUpdateFreq = input.getConstant("cellUpdateFreq");
        std::string potentialName = input.getFilename("potential");
        std::string out_xyz = input.getFilename("out_xyz");
        std::string out_data = input.getFilename("out_data");

        // Optional parameters
        double mu = getOrDefault(input, "mu");
        double f = getOrDefault(input, "f");
        double alpha = getOrDefault(input, "alpha");
        double A_0 = getOrDefault(input, "A_0");
        double kappa = getOrDefault(input, "kappa");
        int seed = static_cast<int>(getOrDefault(input, "seed"));
        std::string position_file = input.getFilename("position_file");
        std::string ensemble_name = input.getFilename("ensemble");

        // Diagnostic: print parsed input values
        std::cout << "\n[Input Summary]" << std::endl;
        std::cout << "Lx = " << Lx << ", Ly = " << Ly << ", N = " << N << std::endl;
        std::cout << "rcut = " << rcut << ", T = " << temperature << std::endl;
        std::cout << "nSteps = " << nSteps << ", eSteps = " << eSteps << std::endl;
        std::cout << "outputFreq = " << outputFreq << ", cellUpdateFreq = " << cellUpdateFreq << std::endl;
        std::cout << "potential = " << potentialName << std::endl;
        std::cout << "mu = " << mu << ", f = " << f << ", alpha = " << alpha
                  << ", A_0 = " << A_0 << ", kappa = " << kappa << std::endl;
        std::cout << "seed = " << seed << std::endl;
        std::cout << "position_file = " << (position_file.empty() ? "[none]" : position_file) << std::endl;
        std::cout << "ensemble = " << (ensemble_name.empty() ? "[default GCMC]" : ensemble_name) << std::endl;
        std::cout << "Output XYZ = " << out_xyz << ", Output Data = " << out_data << std::endl << std::endl;

        // Choose ensemble
        Ensemble ensemble = Ensemble::GCMC; // Default
        if (ensemble_name == "NVT") {
            ensemble = Ensemble::NVT;
        } else if (ensemble_name == "GCMC") {
            ensemble = Ensemble::GCMC;
        }

        // Set up simulation box and particle container
        SimulationBox box(Lx, Ly);
        std::vector<Particle> particles;

        if (!position_file.empty()) {
            initializeParticles_from_file(particles, box, N, position_file);
        } else {
            initializeParticles(particles, box, N, true, seed * 10);
        }

        // Set up logging and RNG
        Logging logger(out_xyz, out_data);
        RNG rng(seed);

        // Set up potential
        PotentialType potentialType = selectPotentialType(potentialName);

        // Set up thermodynamic calculator
        ThermodynamicCalculator calc(
            temperature,
            potentialType,
            rcut,
            mu,
            f,
            alpha,
            A_0,
            kappa
        );

        // Create MC object
        MonteCarlo mc(box, particles, calc, rcut, ensemble, logger, rng);

        // Equilibration simulation
        std::cout << "Equilibration:" << std::endl;
        mc.run(static_cast<size_t>(eSteps),
               static_cast<size_t>(eSteps * 10),
               static_cast<size_t>(cellUpdateFreq));

        // Run simulation
        std::cout << "Running:" << std::endl;
        mc.run(static_cast<size_t>(nSteps),
               static_cast<size_t>(outputFreq),
               static_cast<size_t>(cellUpdateFreq));

        std::cout << "Simulation completed." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

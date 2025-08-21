// main.cpp - Entry point for running the Grand Canonical Monte Carlo simulation
#include "input.h"
#include "initial.h"
#include "MC.h"
#include "logging.h"
#include "thermodynamic_calculator.h"
#include "GibbsMC.h"
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
        double delta = input.getConstant("delta");
        

        // Optional parameters
        double mu = getOrDefault(input, "mu");
        double f = getOrDefault(input, "f");
        double alpha = getOrDefault(input, "alpha");
        double A_0 = getOrDefault(input, "A_0");
        double kappa = getOrDefault(input, "kappa");
        int seed = static_cast<int>(getOrDefault(input, "seed"));
        std::string position_file = input.getFilename("position_file");
        std::string ensemble_name = input.getFilename("ensemble");
        double pressure = getOrDefault(input, "P");
        double delta_V = getOrDefault(input, "delta_V");


        // Choose ensemble
        Ensemble ensemble = Ensemble::GCMC; // Default
        if (ensemble_name == "NVT") {
            ensemble = Ensemble::NVT;
        } else if (ensemble_name == "GCMC") {
            ensemble = Ensemble::GCMC;
        } else if (ensemble_name == "NPT") {
            ensemble = Ensemble::NPT;
        } else if (ensemble_name == "Gibbs") {
            ensemble == Ensemble::Gibbs;
        } else {
            std::cerr<<"Error: Unrecognised ensemble"<<std::endl;
            return 1;
        }
        
        if (ensemble != Ensemble::Gibbs)
        {
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
                pressure,
                potentialType,
                rcut,
                mu,
                f,
                alpha,
                A_0,
                kappa
            );

            // Create MC object
            MonteCarlo mc(box, particles, calc, rcut,delta, delta_V, ensemble, logger, rng);

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
    }else{
        // Required parameters
        double Lx2 = input.getConstant("Lx2");
        double Ly2 = input.getConstant("Ly2");
        int N2 = static_cast<int>(input.getConstant("N2"));
        std::string out_xyz2 = input.getFilename("out_xyz2");
        std::string out_data2 = input.getFilename("out_data2");
        std::string position_file2 = input.getFilename("position_file2");
        // Diagnostic: print parsed input values
        std::cout << "\n[Input Summary]" << std::endl;
        std::cout << "Lx_1 = " << Lx << ", Ly_1 = " << Ly << ", N_1 = " << N << std::endl;
        std::cout << "Lx_2 = " << Lx2 << ", Ly_2 = " << Ly2 << ", N_2 = " << N2 << std::endl;
        std::cout << "rcut = " << rcut << ", T = " << temperature << std::endl;
        std::cout << "nSteps = " << nSteps << ", eSteps = " << eSteps << std::endl;
        std::cout << "outputFreq = " << outputFreq << ", cellUpdateFreq = " << cellUpdateFreq << std::endl;
        std::cout << "potential = " << potentialName << std::endl;
        std::cout << "mu = " << mu << ", f = " << f << ", alpha = " << alpha
                << ", A_0 = " << A_0 << ", kappa = " << kappa << std::endl;
        std::cout << "seed = " << seed << std::endl;
        std::cout << "position_file_1 = " << (position_file.empty() ? "[none]" : position_file) << std::endl;
        std::cout << "position_file_2 = " << (position_file2.empty() ? "[none]" : position_file) << std::endl;
        std::cout << "ensemble = " << (ensemble_name.empty() ? "[default GCMC]" : ensemble_name) << std::endl;
        std::cout << "Output XYZ 1 = " << out_xyz << ", Output Data 1 = " << out_data << std::endl << std::endl;
        std::cout << "Output XYZ 2 = " << out_xyz2 << ", Output Data 2 = " << out_data2 << std::endl << std::endl;

        // Set up simulation box and particle container
        SimulationBox box_1(Lx, Ly);
        std::vector<Particle> particles_1;
        SimulationBox box_2(Lx2, Ly2);
        std::vector<Particle> particles_2;

        if (!position_file.empty()) {
            initializeParticles_from_file(particles_1, box_1, N, position_file);
        } else {
            initializeParticles(particles_1, box_1, N, true, seed * 10);
        }

        if (!position_file2.empty()) {
            initializeParticles_from_file(particles_2, box_2, N2, position_file2);
        } else {
            initializeParticles(particles_2, box_2, N2, true, seed * 5);
        }
        // Set up logging and RNG
        Logging logger_1(out_xyz, out_data);
        Logging logger_2(out_xyz2, out_data2);
        RNG rng(seed);

        // Set up potential
        PotentialType potentialType = selectPotentialType(potentialName);

        // Set up thermodynamic calculator
        ThermodynamicCalculator calc_1(
            temperature,
            pressure,
            potentialType,
            rcut,
            mu,
            f,
            alpha,
            A_0,
            kappa
        );

        ThermodynamicCalculator calc_2(
            temperature,
            pressure,
            potentialType,
            rcut,
            mu,
            f,
            alpha,
            A_0,
            kappa
        );

        GibbsMonteCarlo gmc(box_1, box_2, particles_1, particles_2, calc_1, calc_2, rcut, delta, delta_V, logger_1, logger_2, rng);
        
        // Equilibration simulation
        std::cout << "Equilibration:" << std::endl;
        gmc.run(static_cast<size_t>(eSteps),
            static_cast<size_t>(eSteps * 10),
            static_cast<size_t>(cellUpdateFreq));

        // Run simulation
        std::cout << "Running:" << std::endl;
        gmc.run(static_cast<size_t>(nSteps),
            static_cast<size_t>(outputFreq),
            static_cast<size_t>(cellUpdateFreq));

        std::cout << "Simulation completed." << std::endl;


    }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;

}

    return 0;
}

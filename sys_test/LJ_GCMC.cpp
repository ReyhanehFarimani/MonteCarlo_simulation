#include "../src/initial.h"
#include "../src/input.h"
#include "../src/logging.h"
#include "../src/potential.h"
#include "../src/simulation.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <mpi.h>
#include "test.h"

void LJ_GCMC() {
    // MPI initialization
    MPI_Init(NULL, NULL);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Get the rank of the process
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // Get the number of processes

    const double r_cutoff = 2.5; // Cutoff radius
    const double equilibration_steps = 10000000;
    const double production_steps = 10000000;
    const double mu_increment = 0.1;
    const double seed = 12345;
    double box_length_x = 40 * sqrt(3)/2;
    double box_length_y = 40;
    SimulationBox box(box_length_x, box_length_y);

    int total_tasks = 0;
    for (double T = 10.0; T < 10.1; T += 9.0) {
        for (double mu = -2.0; mu < 2.1; mu += mu_increment) {
            total_tasks++;
        }
    }

    int tasks_per_process = total_tasks / size;
    int remainder = total_tasks % size;

    int start_task = rank * tasks_per_process + std::min(rank, remainder);
    int end_task = start_task + tasks_per_process + (rank < remainder ? 1 : 0);

    int task_index = 0;
    for (double T = 10.0; T < 10.1; T += 1.0) {
        for (double mu = -2.0; mu < 3.1; mu += mu_increment) {
            if (task_index >= start_task && task_index < end_task) {
                Simulation sim(box, PotentialType::LennardJones, SimulationType::GCMC, T, 100, 0.1, 2.5, 0, mu, seed, 1, 100);
                sim.initializeParticles(0);

                std::string mu_str = std::to_string(mu).substr(0, 4); // Limiting to 3 decimal places
                std::string T_str = std::to_string(T).substr(0, 4); // Limiting to 3 decimal places
                Logging logger("LJ_position_" + mu_str + "_" + T_str + ".xyz", "LJ_data_" + mu_str + "_" + T_str + ".txt");

                // Run the simulation
                sim.run(production_steps, equilibration_steps, 1000, logger);
            }
            task_index++;
        }
    }

    // Ensure all processes finish before running the Python script
    MPI_Barrier(MPI_COMM_WORLD);

    // Only the root process runs the Python script
    // if (rank == 0) {
    //     std::system("python ideal_gas_GMCM.py");
    // }

    MPI_Finalize();
}

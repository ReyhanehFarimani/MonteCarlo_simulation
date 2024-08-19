#ifndef LOGGING_H
#define LOGGING_H

#include <fstream>
#include <vector>
#include "initial.h"  // Include the Particle class definition
// Forward declaration of the Simulation class
class Simulation;

/**
 * @brief Manages logging of simulation data to files.
 */
class Logging {
private:
    std::ofstream outFile_position; ///< Output file stream for logging position
    std::ofstream outFile_data;     ///< Output file stream for logging data
    std::string filename_position;  ///< Name of the file to log position
    std::string filename_data;      ///< Name of the file to log data

public:
    /**
     * @brief Constructs a Logging object and opens the file for writing.
     * @param filename The name of the file where data will be logged.
     */
    Logging(const std::string &filename_position, const std::string &filename_data);

    /**
     * @brief Destructor to close the file when the Logging object is destroyed.
     */
    ~Logging();

    /**
     * @brief Logs the positions of all particles to the file.
     * @param particles The vector of particles whose positions will be logged in the xyz format.
     */
    void logPositions_xyz(const std::vector<Particle> &particles, const SimulationBox &box, double r2cut);

    /**
     * @brief Logs the positions of all particles to the file.
     * @param particles The vector of particles whose positions will be logged in the dump format.
     */
    void logPositions_dump(const std::vector<Particle> &particles);

    /**
     * @brief Logs the simularion data like energy pressure, number of particles of the system to the data file.
     * @param energy The total energy of the system.
     * @param timestep The current timestep in the simulation.
     */
    void logSimulationData(const Simulation &sim, int timestep);

    /**
     * @brief Closes the log file.
     */
    void close();
};

#endif // LOGGING_H

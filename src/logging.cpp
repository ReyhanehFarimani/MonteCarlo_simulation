#include "logging.h"
#include "simulation.h"
#include <iomanip>  // For setting precision in output

Logging::Logging(const std::string &filename_position, const std::string &filename_data)
    : filename_position(filename_position), filename_data(filename_data) {
    outFile_position.open(filename_position);  // Open the position file for writing
    outFile_data.open(filename_data);          // Open the data file for writing

    if (!outFile_position) {
        throw std::runtime_error("Could not open file for logging positions: " + filename_position);
    }
    if (!outFile_data) {
        throw std::runtime_error("Could not open file for logging data: " + filename_data);
    }
}

Logging::~Logging() {
    if (outFile_position.is_open()) {
        outFile_position.close();  // Ensure the position file is closed when the object is destroyed
    }
    if (outFile_data.is_open()) {
        outFile_data.close();  // Ensure the data file is closed when the object is destroyed
    }
}

void Logging::logPositions_xyz(const std::vector<Particle> &particles, const SimulationBox &box, double r2cut) {
    outFile_position << particles.size() << "\n";
    outFile_position << "Cell index for each particle\n";

    int numCellsX = static_cast<int>(box.getLx() / r2cut);
    int numCellsY = static_cast<int>(box.getLy() / r2cut);

    for (const auto &particle : particles) {
        int cellX = static_cast<int>(particle.x / r2cut);
        int cellY = static_cast<int>(particle.y / r2cut);
        int cellIndex = cellY * numCellsX + cellX;

        outFile_position << cellIndex << " " << particle.x << " " << particle.y << "\n";
    }
}

void Logging::logPositions_dump(const std::vector<Particle> &particles) {
    if (!outFile_position.is_open()) {
        throw std::runtime_error("Position log file is not open");
    }

    for (const auto &particle : particles) {
        outFile_position << std::fixed << std::setprecision(5) 
                         << particle.x << " " << particle.y << "\n";
    }
    outFile_position << "\n";  // Separate each time step with a blank line
}
void Logging::logSimulationData(const Simulation &sim, int timestep){
    outFile_data << "Timestep:" << timestep << " " << std::fixed << std::setprecision(5)
            << ",\tEnergy:"<< sim.getEnergy() << ",\tTempreture:" << sim.getTemperature() << 
            ",\tPressure:" << sim.getPressure() <<
            ",\tNumberofParticles:"<<sim.getNumParticles()<< "\n";
}

void Logging::close() {
    if (outFile_position.is_open()) {
        outFile_position.close();  // Close the position file when logging is done
    }
    if (outFile_data.is_open()) {
        outFile_data.close();  // Close the data file when logging is done
    }
}
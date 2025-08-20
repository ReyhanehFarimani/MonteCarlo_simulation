#include "logging.h"
#include <iomanip>  // For setting precision in output
#include <cmath>
#include <unistd.h>     // for fsync
#include <fcntl.h>      // for file descriptor


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

    for (const auto &particle : particles) {


        outFile_position << "1" << " " << particle.x << " " << particle.y << "\n";
    }

    outFile_position.flush();
    // fsync(fileno(outFile_position));  // POSIX safe flush to disk

}

void Logging::logPositions_dump(const std::vector<Particle> &particles, const SimulationBox &box, int timestep) {
    if (!outFile_position.is_open()) {
        throw std::runtime_error("Position log file is not open");
    }

    // 1) Header
    outFile_position << "ITEM: TIMESTEP\n"
        << timestep << "\n"
        << "ITEM: NUMBER OF ATOMS\n"
        << particles.size() << "\n";

    // 2) Box bounds (2D → use a dummy z-range [0,0])
    double xlo = 0.0, xhi = box.getLx();
    double ylo = 0.0, yhi = box.getLy();
    double zlo = 0.0, zhi = 0.0;
    outFile_position << "ITEM: BOX BOUNDS pp pp ff\n"
    << xlo << " " << xhi << "\n"
    << ylo << " " << yhi << "\n"
    << zlo << " " << zhi << "\n";

    // 3) Per‐atom lines: id, x, y, z
    outFile_position << "ITEM: ATOMS id x y z\n";
    int id = 1;
    for (const auto &p : particles) {
        outFile_position 
            << id++      << " "
            << p.x       << " "
            << p.y       << " "
            << 0.0       // z‐coordinate dummy
            << "\n";
        }

    outFile_position.flush();
    // fsync(fileno(outFile_position));

}
void Logging::logSimulationData(const std::vector<Particle> &particles, SimulationBox &box, const ThermodynamicCalculator &cal, int timestep){
    outFile_data << "Timestep:" << timestep << " " << std::fixed << std::setprecision(5)
            << ",\tEnergy:"<< cal.computeTotalEnergyCellList(particles, box) << ",\tTempreture:" << cal.getTemperature() << 
            ",\tPressure:" << cal.computePressureCellList(particles, box) <<
            ",\tVolume:" << box.getV() <<
            ",\tdensity of particles:"<<cal.getNumParticles(particles)/ box.getV()<<
            ",\tNumberofParticles:"<<cal.getNumParticles(particles)<< "\n";

    outFile_data.flush();
    // fsync(fileno(outFile_data));
}

void Logging::close() {
    if (outFile_position.is_open()) {
        outFile_position.flush();
        // fsync(fileno(outFile_position));
        outFile_position.close();  // Close the position file when logging is done
    }
    if (outFile_data.is_open()) {
        outFile_position.flush();
        // fsync(fileno(outFile_data));
        outFile_data.close();  // Close the data file when logging is done
    }
}
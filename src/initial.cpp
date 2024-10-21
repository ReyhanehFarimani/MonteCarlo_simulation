#include "initial.h"
#include <cstdlib>  // For rand() and srand()
#include <ctime>    // For time()
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
// Particle class methods
Particle::Particle(double x_init, double y_init)
    : x(x_init), y(y_init) {}

void Particle::updatePosition(double dx, double dy) {
    x += dx;
    y += dy;
}

// Constructor for SimulationBox
SimulationBox::SimulationBox(double Lx, double Ly)
    : Lx(Lx), Ly(Ly), invLx(1.0 / Lx), invLy(1.0 / Ly), volume(Lx * Ly) {}

// Apply periodic boundary conditions to a particle
void SimulationBox::applyPBC(Particle &p) const {
    p.x -= Lx * std::floor(p.x * invLx);
    p.y -= Ly * std::floor(p.y * invLy);
}

// Calculate the minimum image distance between two particles considering PBC
double SimulationBox::minimumImageDistance(const Particle &p1, const Particle &p2) const {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;

    dx -= Lx * std::round(dx * invLx);
    dy -= Ly * std::round(dy * invLy);

    return std::sqrt(dx * dx + dy * dy);
}

// Calculate the minimum image distance squared between two particles considering PBC
double SimulationBox::minimumImageDistanceSquared(const Particle &p1, const Particle &p2) const {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;

    dx -= Lx * std::round(dx * invLx);
    dy -= Ly * std::round(dy * invLy);

    return (dx * dx + dy * dy);
}

// Get the length of the simulation box in the x-direction
double SimulationBox::getLx() const {
    return Lx;
}

// Get the length of the simulation box in the y-direction
double SimulationBox::getLy() const {
    return Ly;
}

// Get the simulation box volume
double SimulationBox::getV() const {
    return volume;
}



void initializeParticles(std::vector<Particle> &particles, const SimulationBox &box, int N, bool random, unsigned int seed, int rank, int world_size) {
    particles.clear();  // Clear the vector before initializing new particles


    // Each rank gets a portion of the total particles
    int local_N = N / world_size;
    int remainder = N % world_size;
    if (rank < remainder) local_N++;

    particles.reserve(local_N);  // Reserve memory for N particles

    // Compute the subdomain bounds for each rank
    double subdomain_x_min = rank * box.getLx() / world_size;
    double subdomain_x_max = (rank + 1) * box.getLx() / world_size;
    double subdomain_y_min = 0;
    double subdomain_y_max = box.getLy();

    if (random) {
        if (seed == 0){
        srand(static_cast<unsigned>(time(0)) + rank * 13 + 71);  // Seed the random number generator
        }
        else{
            srand(seed + 10 + rank * 13 + 71);
        }
        for (int i = 0; i < local_N; ++i) {
            double x = subdomain_x_min + static_cast<double>(rand()) / RAND_MAX * (subdomain_x_max - subdomain_x_min);
            double y = subdomain_y_min + static_cast<double>(rand()) / RAND_MAX * (subdomain_y_max - subdomain_y_min);
            particles.emplace_back(x, y);
        }
    } else {
        // Grid-based particle initialization within the subdomain
        int gridSizeX = static_cast<int>(std::sqrt(local_N));
        int gridSizeY = static_cast<int>(std::sqrt(local_N));
        double spacingX = (subdomain_x_max - subdomain_x_min) / gridSizeX;
        double spacingY = (subdomain_y_max - subdomain_y_min) / gridSizeY;

        for (int i = 0; i < gridSizeX; ++i) {
            for (int j = 0; j < gridSizeY; ++j) {
                if (particles.size() < local_N) {
                    double x = subdomain_x_min + i * spacingX;
                    double y = subdomain_y_min + j * spacingY;
                    particles.emplace_back(x, y);
                }
            }
        }
    }
}

void initializeParticles_from_file(std::vector<Particle> &particles, const SimulationBox &box, int N, const std::string &filename_data, int rank, int world_size) {
    particles.clear();  // Clear the vector before initializing new particles

    // Calculate the number of particles each rank will handle
    int local_N = N / world_size;
    int remainder = N % world_size;
    if (rank < remainder) local_N++;  // Handle uneven distribution

    particles.reserve(local_N);  // Reserve memory only for the local particles

    double subdomain_x_min = rank * box.getLx() / world_size;
    double subdomain_x_max = (rank + 1) * box.getLx() / world_size;

    std::ifstream infile(filename_data);  // Open the file for reading
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << filename_data << std::endl;
        return;
    }

    std::string line;
    int particle_count = 0;

    std::getline(infile, line);  // Skip the first line (header)
    std::getline(infile, line);  // Skip the second line (comment line)

    while (std::getline(infile, line) && particle_count < local_N) {
        std::istringstream iss(line);
        double x, y;
        std::string atom_type;  // For XYZ files, the first column is usually atom type
        if (!(iss >> atom_type >> x >> y)) {
            std::cerr << "Error reading line: " << line << std::endl;
            continue;
        }

        // Only add particles that belong to the rank's subdomain
        if (x >= subdomain_x_min && x < subdomain_x_max) {
            particles.emplace_back(x, y);
            particle_count++;
        }
    }
}

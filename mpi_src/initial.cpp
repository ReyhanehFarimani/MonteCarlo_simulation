#include "initial.h"
#include <cstdlib>  // For rand() and srand()
#include <ctime>    // For time()
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
// Particle class methods
Particle::Particle(double x_init, double y_init, int id_init)
    : x(x_init), y(y_init), id(id_init) {}

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



void initializeParticles(std::vector<Particle> &particles, const SimulationBox &box, int N, bool random, unsigned int seed) {
    // particles.clear();  // Clear the vector before initializing new particles
    // particles.reserve(N);  // Reserve memory for N particles

    // if (random) {
    //     if (seed == 0){
    //     srand(static_cast<unsigned>(time(0)));  // Seed the random number generator
    //     }
    //     else{
    //         srand(seed + 10);
    //     }
    //     for (int i = 0; i < N; ++i) {
    //         double x = static_cast<double>(rand()) / RAND_MAX * box.getLx();
    //         double y = static_cast<double>(rand()) / RAND_MAX * box.getLy();
    //         particles.emplace_back(x, y, i);
    //     }
    // } else {
    //     // Place particles in a simple grid pattern
    //     int gridSize = static_cast<int>(std::sqrt(N));
    //     double spacingX = box.getLx() / gridSize;
    //     double spacingY = box.getLy() / gridSize;

    //     for (int i = 0; i < gridSize ; ++i) {
    //         for (int j = 0; j < gridSize ; ++j) {
    //             if (particles.size() < gridSize * gridSize) {
                    
    //                 double x = i * spacingX;
    //                 double y = j * spacingY;
    //                 // std::cout<<x<<","<<y<<std::endl;
    //                 particles.emplace_back(x, y, i);
    //             }
    //         }
    //     }
    // }
}

void initializeParticles_from_file(std::vector<Particle> &particles, const SimulationBox &box, int N, const std::string &filename_data){
    // particles.clear();  // Clear the vector before initializing new particles
    // particles.reserve(N);  // Reserve memory for N particles

    // double Lx = box.getLx();
    // double Ly = box.getLy();

    // std::ifstream infile(filename_data);  // Open the file for reading
    // if (!infile.is_open()) {
    //     std::cerr << "Error: Could not open file " << filename_data << std::endl;
    //     return;
    // }
        
    // std::string line;
    // int particle_count = 0;

    // // Skip header lines (if any) depending on file format
    // // For an XYZ file, you may want to skip the first 2 lines
    // std::getline(infile, line);  // Skip the first line
    // std::getline(infile, line);  // Skip the second line (comment line)

    // // Read particle data
    // while (std::getline(infile, line) && particle_count < N) {
    //     std::istringstream iss(line);
    //     double x, y;

    //     // For XYZ file format: expect x, y, z as columns 2, 3, 4 (ignoring the particle type in the first column)
        
    //     std::string atom_type;  // For XYZ files, first column might be atom type (string)
    //     if (!(iss >> atom_type >> x >> y)) {
    //         std::cerr << "Error reading line: " << line << std::endl;
    //         continue;
    //     }

    //     // Check if the particle coordinates are inside the simulation box (optional)
    //     if (x >= 0 && x <= Lx && y >= 0 && y <= Ly) {
    //         particles.emplace_back(x, y);  // Add the particle to the list
    //         particle_count++;
    //     } 
    //     else {
    //         std::cerr << "Warning: Particle out of bounds (" << x << ", " << y << ")" << std::endl;
    //     }

    // }
}
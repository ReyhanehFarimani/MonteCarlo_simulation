#include "initial.h"
#include <cstdlib>  // For rand() and srand()
#include <ctime>    // For time()
#include <iostream>

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



void initializeParticles(std::vector<Particle> &particles, const SimulationBox &box, int N, bool random, unsigned int seed) {
    particles.clear();  // Clear the vector before initializing new particles
    particles.reserve(N);  // Reserve memory for N particles

    if (random) {
        if (seed == 0){
        srand(static_cast<unsigned>(time(0)));  // Seed the random number generator
        }
        else{
            srand(seed + 10);
        }
        for (int i = 0; i < N; ++i) {
            double x = static_cast<double>(rand()) / RAND_MAX * box.getLx();
            double y = static_cast<double>(rand()) / RAND_MAX * box.getLy();
            particles.emplace_back(x, y);
        }
    } else {
        // Place particles in a simple grid pattern
        int gridSize = static_cast<int>(std::sqrt(N));
        double spacingX = box.getLx() / gridSize;
        double spacingY = box.getLy() / gridSize;

        for (int i = 0; i < gridSize ; ++i) {
            for (int j = 0; j < gridSize ; ++j) {
                if (particles.size() < gridSize * gridSize) {
                    
                    double x = i * spacingX;
                    double y = j * spacingY;
                    // std::cout<<x<<","<<y<<std::endl;
                    particles.emplace_back(x, y);
                }
            }
        }
    }
}
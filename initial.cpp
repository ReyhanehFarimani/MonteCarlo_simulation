#include "initial.h"
#include <cstdlib>  // For rand() and srand()
#include <ctime>    // For time()


// Particle class methods
Particle::Particle(double x_init, double y_init)
    : x(x_init), y(y_init) {}

void Particle::updatePosition(double dx, double dy) {
    x += dx;
    y += dy;
}

// SimulationBox class methods
SimulationBox::SimulationBox(double Lx, double Ly)
    : Lx(Lx), Ly(Ly) {}

void SimulationBox::applyPBC(Particle &p) const {
    p.x = p.x - Lx * std::floor(p.x / Lx);
    p.y = p.y - Ly * std::floor(p.y / Ly);
}

double SimulationBox::minimumImageDistance(const Particle &p1, const Particle &p2) const {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;

    dx -= Lx * std::round(dx / Lx);
    dy -= Ly * std::round(dy / Ly);

    return std::sqrt(dx * dx + dy * dy);
}

double SimulationBox::getLx() const {
    return Lx;
}

double SimulationBox::getLy() const {
    return Ly;
}



void initializeParticles(std::vector<Particle> &particles, const SimulationBox &box, int N, bool random) {
    particles.clear();  // Clear the vector before initializing new particles
    particles.reserve(N);  // Reserve memory for N particles

    if (random) {
        srand(static_cast<unsigned>(time(0)));  // Seed the random number generator

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

        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                if (particles.size() < N) {
                    double x = i * spacingX;
                    double y = j * spacingY;
                    particles.emplace_back(x, y);
                }
            }
        }
    }
}
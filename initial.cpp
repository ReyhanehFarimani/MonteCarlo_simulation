#include "initial.h"

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

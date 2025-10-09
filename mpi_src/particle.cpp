#include "particle.h"

Particle::Particle(double x_init, double y_init, int id_init)
    : id(id_init), x(x_init), y(y_init) {}

void Particle::updatePosition(double dx, double dy) {
    x += dx; y += dy;
}

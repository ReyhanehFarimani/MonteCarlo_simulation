#ifndef PARTICLE_H
#define PARTICLE_H

class Particle {
public:
    int    id{0};
    double x{0.0}, y{0.0};
    double vx{0.0}, vy{0.0}; // reserved for future use

    Particle() = default;
    Particle(double x_init, double y_init, int id_init = 0);

    void updatePosition(double dx, double dy);
};

#endif // PARTICLE_H
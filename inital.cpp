#include <vector>
#include <cmath>

class Particle {
public:
    double x, y;  // Position of the particle

    // Constructor
    Particle(double x_init = 0.0, double y_init = 0.0)
        : x(x_init), y(y_init) {}

    // Method to update the position of the particle
    void updatePosition(double dx, double dy) {
        x += dx;
        y += dy;
    }
};


class SimulationBox {
private:
    double Lx, Ly;  // Dimensions of the simulation box

public:
    // Constructor
    SimulationBox(double Lx, double Ly)
        : Lx(Lx), Ly(Ly) {}

    // Apply periodic boundary conditions
    void applyPBC(Particle &p) {
        p.x = fmod(p.x + Lx, Lx);
        p.y = fmod(p.y + Ly, Ly);
    }   

   // Calculate minimum image distance considering periodic boundary conditions
    double minimum_image_distance(Particle &p1, Particle &p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;

    dx -= Lx * round(dx / Lx);
    dy -= Ly * round(dy / Ly);

    return sqrt(dx * dx + dy * dy);
    }

    // Getters for the dimensions (useful for other classes)
    double getLx() const { return Lx; }
    double getLy() const { return Ly; }
};


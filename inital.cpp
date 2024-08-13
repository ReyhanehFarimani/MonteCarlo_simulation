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

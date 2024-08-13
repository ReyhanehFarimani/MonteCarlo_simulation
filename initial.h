#ifndef INITIAL_H
#define INITIAL_H

#include <vector>
#include <cmath>

/**
 * @brief Represents a particle in the 2D simulation box.
 */
class Particle {
public:
    double x, y;  ///< Position of the particle

    /**
     * @brief Constructs a Particle with the given initial position.
     * @param x_init Initial x-coordinate of the particle.
     * @param y_init Initial y-coordinate of the particle.
     */
    Particle(double x_init = 0.0, double y_init = 0.0);

    /**
     * @brief Updates the position of the particle.
     * @param dx Displacement in the x-direction.
     * @param dy Displacement in the y-direction.
     */
    void updatePosition(double dx, double dy);
};

/**
 * @brief Manages the 2D simulation box, including periodic boundary conditions.
 */
class SimulationBox {
private:
    double Lx;  ///< Length of the simulation box in the x-direction.
    double Ly;  ///< Length of the simulation box in the y-direction.

public:
    /**
     * @brief Constructs a SimulationBox with the given dimensions.
     * @param Lx Length of the simulation box in the x-direction.
     * @param Ly Length of the simulation box in the y-direction.
     */
    SimulationBox(double Lx, double Ly);

    /**
     * @brief Applies periodic boundary conditions to a particle.
     * @param p The particle to which PBC will be applied.
     */
    void applyPBC(Particle &p) const;

    /**
     * @brief Calculates the minimum image distance between two particles considering PBC.
     * @param p1 The first particle.
     * @param p2 The second particle.
     * @return The shortest distance between p1 and p2 under periodic boundary conditions.
     */
    double minimumImageDistance(const Particle &p1, const Particle &p2) const;

    /**
     * @brief Gets the length of the simulation box in the x-direction.
     * @return The length in the x-direction.
     */
    double getLx() const;

    /**
     * @brief Gets the length of the simulation box in the y-direction.
     * @return The length in the y-direction.
     */
    double getLy() const;
};

#endif // INITIAL_H

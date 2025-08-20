#ifndef INITIAL_H
#define INITIAL_H
#include <string>
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
    double invLx;  ///< Precomputed 1/Lx for performance.
    double invLy;  ///< Precomputed 1/Ly for performance.
    double volume; ///< Precomputed volume of the simulation box.

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
     * @brief Calculates the minimum image distance squared between two particles considering PBC.
     * @param p1 The first particle.
     * @param p2 The second particle.
     * @return The shortest distance squared between p1 and p2 under periodic boundary conditions.
     */
    double minimumImageDistanceSquared(const Particle &p1, const Particle &p2) const;

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

    /**
     * @brief Gets the simulation box volume.
     * @return The simulation box volume.
     */
    double getV() const;

    /**
     * @brief set the simulation box in x_diection.
     * @param lx 
     */
    void setLx(double lx); 
    /**
     * @brief set the simulation box in y_diection.
     * @param ly
     */
    void setLy(double ly); 
    /**
     * @brief set the simulation box volume, keeping the ratio of the box.
     * @param v
     */
    void setV(double v); 
    /**
     * @brief recentering the center of mass based on the old box length:
     * @param p the particle
     * @param lx old lenth of the box
     * @param ly new length og the box
     */
    void recenter(Particle &p, double lx, double ly);
};

/**
 * @brief Initializes N particles in the simulation box.
 * @param particles A vector to store the particles.
 * @param box The simulation box where particles will be placed.
 * @param N The number of particles to initialize.
 * @param random If true, particles are placed randomly; if false, they are placed in a grid.
 */
void initializeParticles(std::vector<Particle> &particles, const SimulationBox &box, int N, bool random = true, unsigned int seed = 0);

/**
 * @brief Initializes N particles in the simulation box, from an xyz file
 * @param particles A vector to store the particles.
 * @param box The simulation box where particles will be placed.
 * @param N The number of particles to initialize.
 * @param filename_data the file from which we are reading simulation data.
 */
void initializeParticles_from_file(std::vector<Particle> &particles, const SimulationBox &box, int N, const std::string &filename_data);

#endif // INITIAL_H

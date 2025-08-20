// serial_src/cell_list.h
#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <vector>
#include "initial.h"
#include "input.h"

/**
 * @brief Spatial cell list for efficient neighbor lookup under PBC in 2D.
 */
class CellList {
public:
    /**
     * @param box Simulation box (defines size and PBC behavior).
     * @param rcut Interaction cutoff distance.
     */
    CellList(SimulationBox& box, double rcut);

    /**
     * @brief change the box of the simulation
     * @param box SimulationBox
     */
    void adjust_box(SimulationBox& box);

    /**
     * @brief Build the cell list for the current particle positions.
     * @param particles Vector of Particle objects.
     */
    void build(const std::vector<Particle>& particles);

    /**
     * @brief Get neighbor indices and squared distances of a given particle (within cutoff).
     * @param idx Index of the query particle.
     * @param particles Vector of Particle objects (for position lookup).
     * @return Vector of pairs (neighbor index, r_sq) excluding self.
     */
    std::vector<std::pair<int,double>> getNeighbors(int idx, const std::vector<Particle>& particles) const;
    /**
     * @brief Get neighbor indices and squared distances of a given particle (within cutoff).
     * @param particle the particle itself.
     * @param particles Vector of Particle objects (for position lookup).
     * @return Vector of pairs (neighbor index, r_sq) excluding self.
     */
    std::vector<std::pair<int,double>> getNeighbors2(const Particle& p, const std::vector<Particle>& particles) const;

    /**
     * @brief Adds a new particle to the appropriate cell.
     * @param p The particle to be added.
     * @param index Index of the particle in the main particle vector.
     */
    void addParticle(const Particle& p, int index);

    /**
     * @brief Removes a particle from its current cell.
     * @param index Index of the particle in the main particle vector.
     */
    void removeParticle(int index);

    private:
    SimulationBox& box_;
    double rcut_;
    double rcutsq_;
    int nx_, ny_;                    // number of cells in x and y
    double cell_dx_, cell_dy_;       // cell dimensions

    // cells_[i] holds indices of particles in cell i = cx + cy*nx_
    std::vector<std::vector<int>> cells_;
    std::vector<int> particle_cell_; // mapping from particle index to cell index

    /**
     * @brief Determine cell index for given (x,y) coordinate.
     */
    int cellIndex(double x, double y) const;
};

#endif // CELL_LIST_H

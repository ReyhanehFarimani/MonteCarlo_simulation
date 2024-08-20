#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include "initial.h"
#include "potential.h"
#include "logging.h"
struct CellListNode {
    int particleIndex;      ///< Index of the particle in the particles vector
    CellListNode* next;     ///< Pointer to the next node in the linked list

    CellListNode(int idx) : particleIndex(idx), next(nullptr) {}
};

enum class SimulationType {
    MonteCarloNVT,
    GCMC,
    // Other simulation types can be added here
    };

/**
 * @brief Selects the simulation type based on a string input.
 * 
 * @param potentialName The name of the potential type as a string.
 * @return The corresponding SimulationType enum value.
 */
SimulationType selectSimulationType(const std::string &simulationName);

/**
 * @brief Manages the entire simulation process, including initialization,
 * running the simulation, and logging results.
 */
class Simulation {
private:
    SimulationBox box;                   ///< The simulation box
    std::vector<Particle> particles;     ///< Particles in the simulation
    PotentialType potentialType;         ///< The type of potential used in the simulation
    double temperature;                  ///< Temperature of the simulation
    int numParticles;                    ///< Number of particles in the simulation
    double maxDisplacement;              ///< Max Displacement for the simulation
    double r2cut;                        ///< Squared distance cutoff for potential calculations
    unsigned int seed;                   ///< Seed for random number generation
    bool useCellList;                   ///< bool use cell list option for optimized energy compution
    int cellListUpdateFrequency;        ///< frquency of cell list updating.
   std::vector<CellListNode*> cellList;  ///< A vector of pointers to linked lists for each cell

public:
    double energy;

   

    /**
     * @brief Constructs a Simulation object with the specified parameters.
     * 
     * @param box The simulation box.
     * @param potentialType The type of potential used in the simulation.
     * @param temperature The temperature of the simulation.
     * @param numParticles The number of particles in the simulation.
     * @param maxDisplacement 
     * @param r2cut The squared distance cutoff for potential calculations.
     */
    Simulation(const SimulationBox &box, PotentialType potentialType, double temperature, int numParticles, double maxDisplacement, double r2cut,unsigned int seed, bool useCellList, int cellListUpdateFrequency);

    /**
     * @brief Initializes particles in the simulation box.
     * 
     * @param randomPlacement If true, particles are placed randomly; otherwise, in a grid.
     */
    void initializeParticles(bool randomPlacement);

    /**
     * @brief Sets the position of a specific particle.
     * @param index The index of the particle to modify.
     * @param x The new x-coordinate of the particle.
     * @param y The new y-coordinate of the particle.
     */
    void setParticlePosition(size_t index, double x, double y);

     /**
     * @brief Calculates and updates the total energy of the system.
     */
    double computeEnergy();



    /**
     * @brief Gets the current energy of the system.
     * @return The total energy of the system.
     */
    double getEnergy() const;

    /**
     * @brief Gets the number of particles in the simulation.
     * @return The number of particles.
     */
    int getNumParticles() const;

    /**
     * @brief Gets the temperature in the simulation.
     * @return Temperature.
     */
    double getTemperature() const;

    /**
     * @brief Perform a single Monte Carlo move.
     * @return True if the move is accepted, false otherwise.
     */
    bool monteCarloMove();
    
    /**
     * @brief Run the simulation.
     * @param numSteps The number of steps to run the simulation.
     * @param equilibrationTime The number of steps for equilibration before logging.
     * @param outputFrequency How often to log the results.
     * @param logger The logging object for output.
     * @param simType The type of simulation to run (e.g., Monte Carlo NVT).
     */
    void run(int numSteps, int equilibrationTime, int outputFrequency, Logging &logger, SimulationType simType);

    void buildCellList();

    /**
     * @brief Clears the current cell list.
     */
    void clearCellList();

    /**
     * @brief Computes the interaction energy of a particle with particles in its cell and neighboring cells.
     * @param particleIndex The index of the particle.
     * @return The computed interaction energy.
     */
    double computeLocalEnergy(int particleIndex) const;

    /**
     * @brief Computes the interaction forces of all particles.
     * @return The computed interaction froces.
     */
    double computeTotalForce() const;


};

#endif // SIMULATION_H

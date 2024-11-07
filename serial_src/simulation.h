#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <string>
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
    SimulationType simtype;             ///< The selected ensemble in the simulation
    PotentialType potentialType;         ///< The type of potential used in the simulation
    double temperature;                  ///< Temperature of the simulation
    int numParticles;                    ///< Number of particles in the simulation
    double maxDisplacement;              ///< Max Displacement for the simulation
    double r2cut;                        ///< Squared distance cutoff for potential calculations
    double rcut;                         ///<  cutoff for potential calculations
    unsigned int seed;                   ///< Seed for random number generation
    bool useCellList;                   ///< bool use cell list option for optimized energy compution
    int cellListUpdateFrequency;        ///< frquency of cell list updating.
    std::vector<CellListNode*> cellList;  ///< A vector of pointers to linked lists for each cell
    float f;                        ///< if one selectys the logarithmic potentials
    float f_d_prime;                    ///< if one selectys the thermal star
    float kappa;                        ///< if one selectys the thermal star
    double mu;                          ///< if one selects the GCMC!
    int numCellsX;                      ///< number of cells in the box in x direction.
    int numCellsY;                      ///< number of cells in the box in y direction.
    int maxNumParticle;                 ///<maximum number of particles in a box.

public:
    double energy;

   

/**
 * @brief Constructs a Simulation object with the specified parameters.
 * 
 * This constructor initializes the simulation environment, setting up the simulation box, interaction potential, temperature, particle count, and other key parameters. It also handles the seeding of the random number generator and prepares the system for simulation.
 * 
 * @param box The simulation box that defines the boundaries of the system.
 * @param potentialType The type of interaction potential used for particle interactions.
 * @param simtype The type of simulation being performed (e.g., Monte Carlo, GCMC).
 * @param temperature The temperature of the system, affecting the acceptance of moves.
 * @param numParticles The initial number of particles in the simulation.
 * @param maxDisplacement The maximum displacement allowed for particles in Monte Carlo steps.
 * @param r2cut The square of the cutoff distance beyond which interactions are neglected.
 * @param f_prime A parameter specific to certain potential types, such as the athermal/thermal star potential.
 * @param f_d_prime A parameter specific to certain potential types, such as the thermal star potential.
 * @param kappa A parameter specific to certain potential types, such as the thermal star potential.
 * @param mu The chemical potential, relevant for grand canonical simulations.
 * @param seed The seed for the random number generator. A value of zero will result in seeding with the current time.
 * @param useCellList Boolean indicating whether to use a cell list for efficient neighbor searching.
 * @param cellListUpdateFrequency The frequency (in simulation steps) at which the cell list should be updated.
 */
Simulation(const SimulationBox &box, PotentialType potentialType, SimulationType simtype, double temperature, int numParticles, double maxDisplacement, double r2cut, float f_prime, float f_d_prime, float kappa, double mu, unsigned int seed, bool useCellList, int cellListUpdateFrequency);

    /**
     * @brief Initializes particles in the simulation box.
     * 
     * @param randomPlacement If true, particles are placed randomly; otherwise, in a grid.
     * @param filename if not empty we'll read the data from the .xyz file.
     */
    void initializeParticles(bool randomPlacement, std::string const &filename);



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
     * @brief remove a particle from cell list.
     * @param particleIndex The index of the particle to be removed
     * @param cellIndex The index of particle to be removed.
     */
    void removeParticleFromCellList(int particleIndex, int cellIndex);

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
     * @brief Attempt to exchange a particle with a resorvoir.
     * @return True if the move is accepted, false otherwise.
     */
    bool monteCarloAddRemove();

    /**
     * @brief Run the simulation.
     * @param numSteps The number of steps to run the simulation.
     * @param equilibrationTime The number of steps for equilibration before logging.
     * @param outputFrequency How often to log the results.
     * @param logger The logging object for output.
     * @param simType The type of simulation to run (e.g., Monte Carlo NVT).
     */
    void run(int numSteps, int equilibrationTime, int outputFrequency, Logging &logger);

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

    /**
     * @brief Calculates the pressure of the system using the virial theorem.
     * 
     * The pressure is computed based on the kinetic energy and the virial of the system.
     * 
     * @return The calculated pressure.
     */
    double calculatePressure() const;

    /**
     * @brief Gets the pressure in the simulation.
     * @return Pressure.
     */
    double getPressure() const;
    /**
    * @brief Calculate tail correction for internal energy in 2D
    * @return Energy Correction.
    */
    double tail_correction_energy_2d() const;
    /**
    * @brief Calculate tail correction for pressure in 2D
    * @return Pressure Correction.
    */
    double tail_correction_pressure_2d() const;
};

#endif // SIMULATION_H

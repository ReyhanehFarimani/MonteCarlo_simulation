#include "simulation.h"
#include "iostream"
#include <cstdlib> // For srand() and rand()

/**
 * @brief Constructs a Simulation object with the specified parameters.
 * 
 * @param box The simulation box.
 * @param potentialType The type of potential used in the simulation.
 * @param temperature The temperature of the simulation.
 * @param numParticles The number of particles in the simulation.
 * @param maxDisplacement The time step for the simulation.
 * @param r2cut The squared distance cutoff for potential calculations.
 */
Simulation::Simulation(const SimulationBox &box, PotentialType potentialType, double temperature, int numParticles, double maxDisplacement, double r2cut, unsigned int seed, bool useCellList, int cellListUpdateFrequency)
    : box(box), potentialType(potentialType), temperature(temperature), numParticles(numParticles), maxDisplacement(maxDisplacement), r2cut(r2cut), energy(0.0), seed(seed), useCellList(useCellList), cellListUpdateFrequency(cellListUpdateFrequency){
    particles.resize(numParticles);
    if (seed != 0) {
        srand(seed); // Seed the random number generator
        
    } else {
        srand(static_cast<unsigned int>(time(nullptr))); // Use current time as seed
    }
    updateEnergy(); // Initialize the energy
}

/**
 * @brief Initializes particles in the simulation box.
 * 
 * @param randomPlacement If true, particles are placed randomly; otherwise, in a grid.
 */
void Simulation::initializeParticles(bool randomPlacement) {
    ::initializeParticles(particles, box, numParticles, randomPlacement, seed);
    if (useCellList) {
        buildCellList();
    }
    updateEnergy();
}

/**
 * @brief Sets the position of a specific particle.
 * @param index The index of the particle to modify.
 * @param x The new x-coordinate of the particle.
 * @param y The new y-coordinate of the particle.
 */
void Simulation::setParticlePosition(size_t index, double x, double y) {
    if (index < particles.size() + 1) {
        particles[index].x = x;
        particles[index].y = y;
    } else {
        std::cerr << "Error: Particle index out of bounds." << std::endl;
    }
    updateEnergy();
}


void Simulation::buildCellList() {
    // Calculate the number of cells along each dimension
    int numCellsX = static_cast<int>(box.getLx() / r2cut);
    int numCellsY = static_cast<int>(box.getLy() / r2cut);

    // Resize the cell list
    cellList.resize(numCellsX * numCellsY);

    // Clear previous cell list data
    for (auto &cell : cellList) {
        cell.clear();
    }

    // Assign particles to cells
    for (size_t i = 0; i < particles.size(); ++i) {
        int cellX = static_cast<int>(particles[i].x / r2cut);
        int cellY = static_cast<int>(particles[i].y / r2cut);
        int cellIndex = cellY * numCellsX + cellX;
        cellList[cellIndex].push_back(i);
    }
}

void Simulation::updateCellListIfNeeded(int step) {
    if (useCellList && step % cellListUpdateFrequency == 0) {
        buildCellList();
    }
}

bool Simulation::monteCarloMove() {

    double oldEnergy = energy; // Energy before the move
    // Randomly select a particle
    size_t particleIndex = rand() % (particles.size() + 1);
    Particle &p = particles[particleIndex];

    // Store the old position
    double oldX = p.x;
    double oldY = p.y;
    // Randomly displace the particle
    double dx = (rand() / double(RAND_MAX) - 0.5) * maxDisplacement;
    double dy = (rand() / double(RAND_MAX) - 0.5) * maxDisplacement;
    p.updatePosition(dx, dy);
    box.applyPBC(p);

    // Calculate the energy difference
    updateEnergy(); // Energy after the move
    
    double deltaE = energy - oldEnergy;
    
    // Decide whether to accept or reject the move
    if (deltaE < 0 || exp(-deltaE / temperature) > (rand() / double(RAND_MAX))) {
        // Accept the move
        return true;
    } else {
        // Reject the move, restore the old position
        p.x = oldX;
        p.y = oldY;
        energy = oldEnergy;
        return false;
    }
}

/**
 * @brief Calculates the total energy of the system.
 * 
 * @return The total energy of the system.
 */

void Simulation::updateEnergy() {
    energy = 0.0;
    
    for (size_t i = 0; i < particles.size() + 1; ++i) {
        for (size_t j = i + 1; j < particles.size() + 1; ++j) {
            double r2 = box.minimumImageDistanceSquared(particles[i], particles[j]);
            
            if (r2 < r2cut) {
                double potential = 0.0;
                switch (potentialType) {
                    case PotentialType::LennardJones:
                        potential = lennardJonesPotential(r2);
                        
                        break;
                    case PotentialType::WCA:
                        potential = wcaPotential(r2);
                        break;
                    case PotentialType::Yukawa:
                        potential = yukawaPotential(r2);
                        break;
                }
                energy += potential;
                
            }
        }
        
    } // Update the energy member variable
    
}

/**
 * @brief Runs the simulation for a specified number of steps.
 * 
 * @param numSteps The number of steps to run the simulation.
 * @param equilibrationTime The number of steps for equilibration before logging.
 * @param outputFrequency How often to log the results.
 * @param logger The logging object for output.
 */
void Simulation::run(int numSteps, int equilibrationTime, int outputFrequency, Logging &logger, SimulationType simType) {
    int acceptedMoves = 0;
    updateEnergy();
    for (int step = 0; step < numSteps; ++step) {
        if (simType == SimulationType::MonteCarloNVT) {
            if (monteCarloMove()) {
                acceptedMoves++;
            }
        }
        updateCellListIfNeeded(step);
        // Other simulation types can be added here as additional conditions
        
        // Optionally log data
        if (step >= equilibrationTime && step % outputFrequency == 0) {
            logger.logPositions_xyz(particles);
            updateEnergy();
            logger.logSimulationData(*this, step);
        }
    }

    if (simType == SimulationType::MonteCarloNVT) {
        std::cout << "Monte Carlo NVT simulation completed with " 
                  << acceptedMoves << " accepted moves out of " << numSteps << " steps." << std::endl;
    }
}

double Simulation::getEnergy() const {
    return energy;
}

int Simulation::getNumParticles() const {
    return numParticles;
}

double Simulation::getTemperature() const {
    return temperature;
}



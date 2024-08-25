#include "simulation.h"
#include "iostream"
#include <cstdlib> // For srand() and rand()
#include <cassert>
/**
 * @brief Initializes a Simulation object with the given parameters.
 * 
 * This constructor sets up the simulation environment, including the simulation box, interaction potential, temperature, and other relevant parameters. It also seeds the random number generator and computes the initial energy of the system.
 * 
 * @param box The simulation box that defines the boundaries of the system.
 * @param simtype The type of simulation being performed (e.g., Monte Carlo, GCMC).
 * @param potentialType The type of interaction potential used for particle interactions.
 * @param temperature The temperature of the system, which influences the acceptance of moves.
 * @param numParticles The initial number of particles in the simulation.
 * @param maxDisplacement The maximum displacement allowed for particle moves in Monte Carlo steps.
 * @param r2cut The square of the cutoff distance beyond which interactions are ignored.
 * @param f_prime An additional parameter related to the specific potential type (e.g., athermal star potential).
 * @param mu The chemical potential, relevant for grand canonical simulations.
 * @param seed The seed for the random number generator. If zero, the current time is used as the seed.
 * @param useCellList A boolean flag indicating whether to use a cell list for efficient neighbor searching.
 * @param cellListUpdateFrequency The frequency (in simulation steps) at which the cell list is updated.
 */
Simulation::Simulation(const SimulationBox &box, PotentialType potentialType, SimulationType simtype, double temperature, int numParticles, double maxDisplacement, double r2cut, float f_prime, double mu, unsigned int seed, bool useCellList, int cellListUpdateFrequency)
    : box(box), simtype(simtype), potentialType(potentialType), temperature(temperature), numParticles(numParticles), maxDisplacement(maxDisplacement), r2cut(r2cut * r2cut), f_prime(f_prime), mu(mu), energy(0.0), seed(seed), useCellList(useCellList), cellListUpdateFrequency(cellListUpdateFrequency) {
    particles.resize(numParticles);
    if (seed != 0) {
        srand(seed); // Seed the random number generator
    } else {
        srand(static_cast<unsigned int>(time(nullptr))); // Use the current time as the seed
    }
    energy = computeEnergy(); // Compute and store the initial energy of the system
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
    energy = computeEnergy();
}

/**
 * @brief Sets the position of a specific particle.
 * @param index The index of the particle to modify.
 * @param x The new x-coordinate of the particle.
 * @param y The new y-coordinate of the particle.
 */
void Simulation::setParticlePosition(size_t index, double x, double y) {
    if (index < particles.size()) {
        particles[index].x = x;
        particles[index].y = y;
    } else {
        std::cerr << "Error: Particle index out of bounds." << std::endl;
    }
    if (useCellList) {
        buildCellList();
    }
    energy = computeEnergy();
}

/**
 * @brief Perform a single Monte Carlo move.
 * @return True if the move is accepted, false otherwise.
 */
bool Simulation::monteCarloMove() {
    // Randomly select a particle
    size_t particleIndex = rand() % (particles.size());
    Particle &p = particles[particleIndex];

    // Store the old position
    double oldX = p.x;
    double oldY = p.y;
    // Calculate the initial local energy
    double initialEnergy;
    if (useCellList)
        initialEnergy = computeLocalEnergy(particleIndex);
    else
        initialEnergy = energy;
    // double test_e_new = computeEnergy();
    // Randomly displace the particle
    double dx = (rand() / double(RAND_MAX) - 0.5) * maxDisplacement;
    double dy = (rand() / double(RAND_MAX) - 0.5) * maxDisplacement;
    p.updatePosition(dx, dy);
    box.applyPBC(p);

    // Calculate the new local energy
    double newEnergy;
    if (useCellList)
        newEnergy = computeLocalEnergy(particleIndex);
    else
        newEnergy = computeEnergy();
    // double test_e_old = computeEnergy();
    // Calculate the energy difference
    double deltaE = newEnergy - initialEnergy;
    // std::cout<<deltaE<<"   "<<test_e_old - test_e_new<<std::endl;
    // Decide whether to accept or reject the move
    if (deltaE < 0 || exp(-deltaE / temperature) > (rand() / double(RAND_MAX))) {
        // Accept the move
        energy += deltaE;
        return true;
    } else {
        // Reject the move, restore the old position
        p.x = oldX;
        p.y = oldY;
        return false;
    }
}

/**
 * @brief Attempt to exchange a particle with a resorvoir.
 * @return True if the move is accepted, false otherwise.
 */
bool Simulation::monteCarloAddRemove() {
    double r_decide = rand() / double(RAND_MAX);
    if (r_decide < 0.5) { 
        // random position for the partile:
        double x = static_cast<double>(rand()) / RAND_MAX * box.getLx();
        double y = static_cast<double>(rand()) / RAND_MAX * box.getLy();
        Particle newParticle(x, y);

        double dE = 0;

        if (useCellList){
            //computing particle cell index:
            double rcut = sqrt(r2cut);
            int numCellsX = static_cast<int>(box.getLx() / rcut);
            int numCellsY = static_cast<int>(box.getLy() / rcut);
            // Ensure that there are enough cells to cover the box
            if (box.getLx() > numCellsX * rcut) numCellsX++;
            if (box.getLy() > numCellsY * rcut) numCellsY++;
            int cellX = static_cast<int>(newParticle.x / rcut);
            int cellY = static_cast<int>(newParticle.y / rcut);
            int cellIndex = cellY * numCellsX + cellX;
            assert(cellIndex >= 0 && cellIndex < cellList.size()); // Check bounds
            for (int offsetY = -2; offsetY <= 2; ++offsetY) {
                int neighborCellY = (cellY + offsetY + numCellsY) % numCellsY;
                for (int offsetX = -2; offsetX <= 2; ++offsetX) {
                    int neighborCellX = (cellX + offsetX + numCellsX) % numCellsX;
                    int neighborCellIndex = neighborCellY * numCellsX + neighborCellX;
                    for (CellListNode* node = cellList[neighborCellIndex]; node != nullptr; node = node->next) {
                        double r2 = box.minimumImageDistanceSquared(newParticle, particles[node->particleIndex]);
                        if (r2 < r2cut) {
                            switch (potentialType) {
                                case PotentialType::LennardJones:
                                    dE += lennardJonesPotential(r2);
                                    break;
                                case PotentialType::WCA:
                                    dE += wcaPotential(r2);
                                    break;
                                case PotentialType::Yukawa:
                                    dE += yukawaPotential(r2);
                                    break;
                                case PotentialType::AthermalStar:
                                    dE += athermalStarPotential(r2, f_prime);
                                    break;
                            }
                        }
                    }
                }
            }
        }//adding particle to the new cell and compute the energy (using cell list)
        else {
            for (size_t i = 0; i < particles.size(); ++i) {
                double r2 = box.minimumImageDistanceSquared(newParticle, particles[i]);
                if (r2 < r2cut) {
                    switch (potentialType) {
                        case PotentialType::LennardJones:
                            dE += lennardJonesPotential(r2);
                            break;
                        case PotentialType::WCA:
                            dE += wcaPotential(r2);
                            break;
                        case PotentialType::Yukawa:
                            dE += yukawaPotential(r2);
                            break;
                        case PotentialType::AthermalStar:
                            dE += athermalStarPotential(r2, f_prime);
                            break;
                    }
                }
            }
        }//adding particle to the new cell and compute the energy (no cell list)
        double acc = exp(- dE/ temperature) * mu * box.getLx() * box.getLy() / (particles.size());


    } //if adding a particle
    else {
        return 0;
    } //if for removing a particle
    return 0;
}

/**
 * @brief Calculates the total energy of the system.
 * 
 * @return The total energy of the system.
 */

double Simulation::computeEnergy() {
    double tmp_energy = 0.0;
    
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
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
                    case PotentialType::AthermalStar:
                        potential = athermalStarPotential(r2, f_prime);
                        break;
                }
                tmp_energy += potential;
                
            }
        }
        
    } 
    return tmp_energy;
    
}


/**
 * @brief Runs the simulation for a specified number of steps.
 * 
 * @param numSteps The number of steps to run the simulation.
 * @param equilibrationTime The number of steps for equilibration before logging.
 * @param outputFrequency How often to log the results.
 * @param logger The logging object for output.
 */
void Simulation::run(int numSteps, int equilibrationTime, int outputFrequency, Logging &logger, SimulationType simType = SimulationType::MonteCarloNVT) {
    if(simType == SimulationType::MonteCarloNVT){
        int acceptedMoves = 0;
        energy = computeEnergy();
        
        // Calculate energy at each step using the appropriate method
        for (int step = 0; step < numSteps + equilibrationTime; ++step) {
            if (useCellList){
            if (step%cellListUpdateFrequency == 0)
            {
                buildCellList();
                energy = computeEnergy();
            }
            }
            if (simType == SimulationType::MonteCarloNVT) {
                if (monteCarloMove()) {
                    acceptedMoves++;
                }
            }

            // Other simulation types can be added here as additional conditions
            
            // Optionally log data
            if (step >= equilibrationTime && step % outputFrequency == 0) {
                if (fabs(energy - computeEnergy())> 1){
                std::cerr<<energy<<"   local energy computation is not working.   "<<computeEnergy()<<std::endl;
                }   
                logger.logPositions_xyz(particles, box, r2cut);
                logger.logSimulationData(*this, step);
            }
        }
        std::cout << "Monte Carlo NVT simulation completed with " 
                    << acceptedMoves << " accepted moves out of " << (numSteps + equilibrationTime) << " steps." << std::endl;
        
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


void Simulation::clearCellList() {
    for (auto& cell : cellList) {
        CellListNode* current = cell;
        while (current != nullptr) {
            CellListNode* toDelete = current;
            current = current->next;
            delete toDelete;
        }
        cell = nullptr;
    }
}



void Simulation::buildCellList() {
    clearCellList();  // Clear the previous cell list
    double rcut = sqrt(r2cut);
    int numCellsX = static_cast<int>(box.getLx() / rcut);
    int numCellsY = static_cast<int>(box.getLy() / rcut);
    // Ensure that there are enough cells to cover the box
    if (box.getLx() > numCellsX * rcut) numCellsX++;
    if (box.getLy() > numCellsY * rcut) numCellsY++;
    // Resize the cells vector based on numCellsX and numCellsY
    cellList.resize(numCellsX * numCellsY, nullptr);

    for (size_t i = 0; i < particles.size(); ++i) {
        int cellX = static_cast<int>(particles[i].x / rcut);
        int cellY = static_cast<int>(particles[i].y / rcut);
        int cellIndex = cellY * numCellsX + cellX;
        assert(cellIndex >= 0 && cellIndex < cellList.size()); // Check bounds

        CellListNode* newNode = new CellListNode(i);
        newNode->next = cellList[cellIndex];
        cellList[cellIndex] = newNode;
    }
}

double Simulation::computeLocalEnergy(int particleIndex) const {
    double localEnergy = 0.0;
    const Particle& p = particles[particleIndex];

    double rcut = sqrt(r2cut);
    int numCellsX = static_cast<int>(box.getLx() / rcut);
    int numCellsY = static_cast<int>(box.getLy() / rcut);
    
    // Ensure that there are enough cells to cover the box
    if (box.getLx() > numCellsX * rcut) numCellsX++;
    if (box.getLy() > numCellsY * rcut) numCellsY++;


    // Get the current cell of the particle
    int cellX = static_cast<int>(p.x / rcut);
    int cellY = static_cast<int>(p.y / rcut);


    // Iterate over neighboring cells
    for (int offsetY = -2; offsetY <= 2; ++offsetY) {
        int neighborCellY = (cellY + offsetY + numCellsY) % numCellsY;
        for (int offsetX = -2; offsetX <= 2; ++offsetX) {
            int neighborCellX = (cellX + offsetX + numCellsX) % numCellsX;
            int neighborCellIndex = neighborCellY * numCellsX + neighborCellX;

            for (CellListNode* node = cellList[neighborCellIndex]; node != nullptr; node = node->next) {
                if (node->particleIndex != particleIndex) {
                    // std::cout<<node->particleIndex<<std::endl;
                    double r2 = box.minimumImageDistanceSquared(p, particles[node->particleIndex]);
                    if (r2 < r2cut) {
                        switch (potentialType) {
                            case PotentialType::LennardJones:
                                localEnergy += lennardJonesPotential(r2);
                                break;
                            case PotentialType::WCA:
                                localEnergy += wcaPotential(r2);
                                break;
                            case PotentialType::Yukawa:
                                localEnergy += yukawaPotential(r2);
                                break;
                            case PotentialType::AthermalStar:
                                localEnergy += athermalStarPotential(r2, f_prime);
                                break;
                        }
                    }
                }
            }
        }
    }
    // std::cout<<localEnergy<<std::endl;
    return localEnergy;
}


SimulationType selectSimulationType(const std::string &simulationName) {
    if (simulationName == "NVT") return SimulationType::MonteCarloNVT;
    if (simulationName == "GCMC") return SimulationType::GCMC;
    throw std::invalid_argument("Unknown potential type: " + simulationName);
}


double Simulation::computeTotalForce() const{
    double forceSum = 0.0;
    for (size_t p1 = 0; p1 < particles.size(); ++p1){
        for (size_t p2 = p1 + 1; p2 < particles.size(); ++p2){
            double r2 = box.minimumImageDistanceSquared(particles[p1], particles[p2]);
            if (r2<r2cut){
                switch (potentialType){
                    case PotentialType::LennardJones:
                        forceSum += lennardJonesForceDotR(r2);
                        break;
                    case PotentialType::WCA:
                        forceSum += wcaForceDotR(r2);
                        break;
                    case PotentialType::Yukawa:
                        forceSum += yukawaForceDotR(r2);
                        break;
                    case PotentialType::AthermalStar:
                        forceSum += athermalStarForceDotR(r2, f_prime);
                        break;
                }
            }
        }
    }
    return forceSum;
}

/**
 * @brief Calculates the pressure of the system using the virial theorem.
 * 
 * The pressure is calculated as P = (N * k_B * T + virial) / V,
 * where N is the number of particles, k_B is the Boltzmann constant,
 * T is the temperature, the virial is the sum of r_ij * f_ij over all pairs,
 * and V is the volume of the simulation box.
 * 
 * @return The calculated pressure.
 */
double Simulation::calculatePressure() const {
    double V = box.getLx() * box.getLy();
    double virialPressure = computeTotalForce() / (2 * V);
    double rho = numParticles / V;
    double idealPressure = rho * temperature;
    return idealPressure + virialPressure;
}


    
/**
* @brief Gets the pressure in the simulation.
* @return Pressure.
*/
double Simulation::getPressure() const {
    double P = calculatePressure();
    return P;
}

/**
* @brief Calculate tail correction for internal energy in 2D
* 
* in two dimension considering we have same densities, we are going to have,  
* the enrgy correction as follows:
* u_{tail} = \frac{1}{2} 2 \pi int_{r_c}^{\inf} dr r \rho u(r) = \pi * rho * int_{r_c}^{inf} dr r u(r) 
* 
* @return Energy Correction.
*/
double Simulation::tail_correction_energy_2d() const{
    double answer = M_PI * numParticles/(box.getLx() * box.getLy());
    switch (potentialType){
    case PotentialType::LennardJones:{
        double sr2 = 1/r2cut;
        double sr4 = sr2 * sr2;
        double sr10 = sr4 * sr4 * sr2;
        answer *= 4.0 * (sr4/4.0 - sr10/10.0);
        break;
    }
    case PotentialType::WCA:
        answer = 0;
        break;
    case PotentialType::AthermalStar:{
        double I  =  -exp(1 - r2cut)/4.0 * f_prime;
        answer *= I;
        break;
    }
    default:
        break;
    }
    return answer;
}

/**
* @brief Calculate tail correction for pressure in 2D.
* 
* g(r) = 1 -> P_cut = \frac{\rho^2 pi}{2} \int_{r_{cut}}^{\inf} f(r) r^2 dr
* 
* @return Pressure Correction.
*/
double Simulation::tail_correction_pressure_2d() const{
    double answer = M_PI * (numParticles/(box.getLx() * box.getLy())) * (numParticles/(box.getLx() * box.getLy())) /2;
    switch (potentialType){
    case PotentialType::LennardJones:{
        double sr2 = 1/r2cut;
        double sr4 = sr2 * sr2;
        double sr10 = sr4 * sr4 * sr2;
        answer *= 48.0 * (sr4/8.0 - sr10/10.0);
        break;
    }
    case PotentialType::WCA:
        answer = 0;
        break;
    case PotentialType::AthermalStar:{
        double I  =  -exp(1.0 - r2cut) / 2.0 * (1.0 + r2cut) * f_prime;
        answer *= I;
        break;
    }
    default:
        break;
    }
    return answer;
}

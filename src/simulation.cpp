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
    : box(box), simtype(simtype), potentialType(potentialType), temperature(temperature), numParticles(numParticles), maxDisplacement(maxDisplacement), rcut(r2cut), r2cut(r2cut * r2cut), f_prime(f_prime), mu(exp(mu/temperature)), energy(0.0), seed(seed), useCellList(useCellList), cellListUpdateFrequency(cellListUpdateFrequency) , numCellsX(1), numCellsY(1){
    particles.resize(numParticles);
    if (seed != 0) {
        srand(seed); // Seed the random number generator
    } else {
        srand(static_cast<unsigned int>(time(nullptr))); // Use the current time as the seed
    }
    energy = computeEnergy(); // Compute and store the initial energy of the system
    numCellsX = static_cast<int>(box.getLx() / rcut);
    numCellsY = static_cast<int>(box.getLy() / rcut);
    // Ensure that there are enough cells to cover the box
    if (box.getLx() > numCellsX * rcut) numCellsX++;
    if (box.getLy() > numCellsY * rcut) numCellsY++;
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
    double r_decide = 0;
    if (simtype == SimulationType::GCMC)
        r_decide = rand() / double(RAND_MAX);
    if (r_decide < 0.98){
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
    else {
        bool acceptance = monteCarloAddRemove();
        return acceptance;
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
        //computing particle cell index:
        int cellX = static_cast<int>(newParticle.x / rcut);
        int cellY = static_cast<int>(newParticle.y / rcut);
        int cellIndex = cellY * numCellsX + cellX;

        if (useCellList){
            
            for (int offsetY = -2; offsetY <= 2; ++offsetY) {
                int neighborCellY = (cellY + offsetY + numCellsY) % numCellsY;
                for (int offsetX = -2; offsetX <= 2; ++offsetX) {
                    int neighborCellX = (cellX + offsetX + numCellsX) % numCellsX;
                    int neighborCellIndex = neighborCellY * numCellsX + neighborCellX;
                    for (CellListNode* node = cellList[neighborCellIndex]; node != nullptr; node = node->next) {
                        double r2 = box.minimumImageDistanceSquared(newParticle, particles[node->particleIndex]);
                        if (r2 < r2cut) {
                            dE += computePairPotential(r2, potentialType, f_prime); 
                        }
                    }
                }
            }
        }//adding particle to the new cell and compute the energy (using cell list)
        else {
            for (size_t i = 0; i < particles.size(); ++i) {
                double r2 = box.minimumImageDistanceSquared(newParticle, particles[i]);
                if (r2 < r2cut) {
                    dE += computePairPotential(r2, potentialType, f_prime); 
                }
            }
        }//adding particle to the new cell and compute the energy (no cell list)
        double acc = exp( -dE/ temperature) * box.getV() / (particles.size() + 1) * mu;
        // std::cout<<acc<<std::endl;
        if (acc > (rand() / double(RAND_MAX))){
            particles.emplace_back(newParticle);
            numParticles ++;
            energy += dE;
            if (useCellList) {
                CellListNode* newNode = new CellListNode(numParticles - 1);
                newNode->next = cellList[cellIndex];
                cellList[cellIndex] = newNode;
            }
            return true;     
        } //adding a particle
        else 
            return false;

    } //if adding a particle
    else {
        if (particles.size()>5){
            //choosing a particle
            size_t particleIndex = rand() % (particles.size());
            Particle &p = particles[particleIndex];
            double dE = 0;
            //computing particle cell index:
            int cellX = static_cast<int>(p.x / rcut);
            int cellY = static_cast<int>(p.y / rcut);
            int cellIndex = cellY * numCellsX + cellX;

            if (useCellList){
                dE = -computeLocalEnergy(particleIndex);
            }//removing particle to the new cell and compute the energy (using cell list)
            else {
                for (size_t i = 0; i < particles.size(); ++i) {
                    if (i != particleIndex){
                        double r2 = box.minimumImageDistanceSquared(p, particles[i]);
                        if (r2 < r2cut) {
                            dE -= computePairPotential(r2, potentialType, f_prime); 
                        }
                    }
                }
            }//removing particle to the new cell and compute the energy (no cell list)
            double acc = exp( -dE/ temperature) / box.getV() * (particles.size()) * mu;
            
            if (acc > (rand() / double(RAND_MAX))){
                if (useCellList) {
                    removeParticleFromCellList(particleIndex, cellIndex);
                }
                if (particleIndex < particles.size()) 
                    particles.erase(particles.begin() + particleIndex);
                else
                    std::cerr<<"particle index out of bound!, unable to remove the particle!"<<std::endl;
                numParticles --;
                energy += dE;
                return true;     
            } //removing a particle
            else 
                return false;
            }
        return false;
    } //if removing a plarticle a particle
    
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
                double potential = computePairPotential(r2, potentialType, f_prime);
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
            if (monteCarloMove()) {
                acceptedMoves++;
                
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

double Simulation::getEnergy() const {
    return energy;
}

int Simulation::getNumParticles() const {
    assert(particles.size() == numParticles);
    return particles.size();
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
                        localEnergy += computePairPotential(r2, potentialType, f_prime);
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
                forceSum += computePairForce(r2, potentialType, f_prime);
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
    double V = box.getV();
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
    double answer = M_PI * numParticles/(box.getV());
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
    double answer = M_PI * (numParticles/(box.getV())) * (numParticles/(box.getV())) /2;
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


void Simulation::removeParticleFromCellList(int particleIndex, int cellIndex) {
    CellListNode* current = cellList[cellIndex];
    CellListNode* prev = nullptr;

    while (current != nullptr) {
        if (current->particleIndex == particleIndex) {
            if (prev == nullptr) {
                cellList[cellIndex] = current->next; // Removing head of the list
            } else {
                prev->next = current->next; // Bypassing the current node
            }
            delete current; // Free the node's memory
            break;
        }
        prev = current;
        current = current->next;
    }

    // Adjust the indices in the cell list to reflect the removed particle
    for (CellListNode* node : cellList) {
        if (node && node->particleIndex > particleIndex) {
            node->particleIndex--;
        }
    }
}



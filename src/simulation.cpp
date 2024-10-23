#include "simulation.h"
#include "iostream"
#include "mpi.h"
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
Simulation::Simulation(const SimulationBox &box, PotentialType potentialType, SimulationType simtype, double temperature, int numParticles, double maxDisplacement, double r2cut, float f_prime, float f_d_prime, float kappa, double mu, unsigned int seed, bool useCellList, int cellListUpdateFrequency)
    : box(box), simtype(simtype), potentialType(potentialType), temperature(temperature), numParticles(numParticles), maxDisplacement(maxDisplacement), rcut(r2cut), r2cut(r2cut * r2cut), f_prime(f_prime), f_d_prime(f_d_prime), kappa(kappa), mu(exp(mu/temperature)), energy(0.0), seed(seed), useCellList(useCellList), cellListUpdateFrequency(cellListUpdateFrequency) , numCellsX(1), numCellsY(1){
    // Initialize MPI variables
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    numCellsX = static_cast<int>(box.getLx() / rcut);
    numCellsY = static_cast<int>(box.getLy() / rcut);
    // Ensure that there are enough cells to cover the box
    if (box.getLx() > numCellsX * rcut) numCellsX++;
    if (box.getLy() > numCellsY * rcut) numCellsY++;
    maxNumParticle = int(box.getV() * 3);


    // Calculate local number of particles each rank will handle
    int local_numParticles = numParticles / world_size;
    int remainder = numParticles % world_size;
    
    // Rank 0 handles remainder if the number of particles isn't evenly divisible
    if (world_rank < remainder) {
        local_numParticles++;
    }

    particles.resize(local_numParticles);
    if (seed != 0) {
        srand(seed + world_rank * 13 + 71); // Seed the random number generator + word rank
    } else {
        srand(static_cast<unsigned int>(time(nullptr)) + world_rank * 13 + 71); // Use the current time as the seed and the world rank to make sure each rank is independent
        
    }
    // std::cout<<seed + world_rank * 13 + 71<<std::endl;
     // Compute and store the initial energy of the system
    
    
}


// In Simulation::initializeParticles
void Simulation::initializeParticles(bool randomPlacement, const std::string &filename) {
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (filename.empty()) {
        ::initializeParticles(particles, box, numParticles, randomPlacement, seed, world_rank, world_size);
    } else {
        ::initializeParticles_from_file(particles, box, numParticles, filename, world_rank, world_size);
    }

    // if (useCellList) {
    //     buildCellList();
    // }
    // Energy computation can be deferred or handled locally
    energy = computeEnergy();
    // std::cout<<"totalenergy"<<energy<<std::endl;
}

void Simulation::gatherAllParticles(std::vector<Particle> &allParticles) {
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // First, gather the number of particles on each rank
    int local_numParticles = particles.size();
    std::vector<int> allNumParticles(world_size);
    MPI_Gather(&local_numParticles, 1, MPI_INT, allNumParticles.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Now, gather the particles
    // Serialize the particles into a contiguous buffer
    std::vector<double> localParticleData(local_numParticles * 2);  // Each particle has x and y
    for (int i = 0; i < local_numParticles; ++i) {
        localParticleData[2 * i] = particles[i].x;
        localParticleData[2 * i + 1] = particles[i].y;
    }

    // Prepare to receive data on rank 0
    std::vector<int> displs;
    std::vector<int> recvCounts;
    std::vector<double> allParticleData;  // Will hold all particle data on rank 0

    if (world_rank == 0) {
        recvCounts.resize(world_size);
        displs.resize(world_size);
        int totalParticles = 0;
        for (int i = 0; i < world_size; ++i) {
            recvCounts[i] = allNumParticles[i] * 2;  // Each particle has x and y
            displs[i] = totalParticles;
            totalParticles += recvCounts[i];
        }
        allParticleData.resize(totalParticles);
    }

    // Gather the particle data
    MPI_Gatherv(localParticleData.data(), localParticleData.size(), MPI_DOUBLE,
                allParticleData.data(), recvCounts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // On rank 0, reconstruct the Particle objects
    if (world_rank == 0) {
        allParticles.clear();
        int totalParticles = allParticleData.size() / 2;
        allParticles.reserve(totalParticles);
        for (int i = 0; i < totalParticles; ++i) {
            double x = allParticleData[2 * i];
            double y = allParticleData[2 * i + 1];
            allParticles.emplace_back(x, y);
            // std::cout<<i<<" "<<x<<" "<<y<<std::endl;
            // std::cout<<allParticles.size()<<std::endl;
        }
    }
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
        if (numParticles < maxNumParticle){
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
                                dE += computePairPotential(r2, potentialType, f_prime, f_d_prime, kappa); 
                            }
                        }
                    }
                }
            }//adding particle to the new cell and compute the energy (using cell list)
            else {
                for (size_t i = 0; i < particles.size(); ++i) {
                    double r2 = box.minimumImageDistanceSquared(newParticle, particles[i]);
                    if (r2 < r2cut) {
                        dE += computePairPotential(r2, potentialType, f_prime, f_d_prime, kappa); 
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
        }// max number of particles

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
                            dE -= computePairPotential(r2, potentialType, f_prime, f_d_prime, kappa); 
                        }
                    }
                }
            }//removing particle to the new cell and compute the energy (no cell list)
            double acc = exp( -dE/ temperature) / box.getV() * (particles.size())/ mu;
            
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


double Simulation::computeLocalEnergy_parallel() {
    //helper function for computing the energy in parallel first we compute the energy in each subdomain!
    double localEnergy = 0.0;
    
    // Loop over all local particles in this rank
    for (size_t i = 0; i < particles.size(); ++i) {
        // std::cout<<particles[i].x<<std::endl;
        for (size_t j = i + 1; j < particles.size(); ++j) {
            double r2 = box.minimumImageDistanceSquared(particles[i], particles[j]);
            
            if (r2 < r2cut) {
                double potential = computePairPotential(r2, potentialType, f_prime, f_d_prime, kappa);
                localEnergy += potential;
            }
        }
    }
    // std::cout<<localEnergy<<std::endl;
    return localEnergy;
}

void Simulation::identifyBoundaryParticles(std::vector<Particle>& boundaryParticles) {
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    boundaryParticles.clear();

    // Compute the subdomain bounds for each rank
    double subdomain_x_min = world_rank * box.getLx() / world_size;
    double subdomain_x_max = (world_rank + 1) * box.getLx() / world_size;
    // std::cout<<subdomain_x_min<<std::endl;
    // Loop through all local particles and identify those near the subdomain boundaries
    for (const auto& particle : particles) {
        // Check if the particle is near the left or right boundary of the subdomain
        if ((particle.x < subdomain_x_min + rcut) || 
            (particle.x > subdomain_x_max - rcut)) {
            // If the particle is near the left or right boundary, add it to the boundaryParticles list
            boundaryParticles.push_back(particle);
            // std::cout<<particle.x<<" "<<particle.y<<std::endl;
        }
    }
    
}

void exchangeBoundaryParticlesWithNextRank(std::vector<Particle> &boundaryParticles, std::vector<Particle> &receivedParticles) {
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // The next rank (with periodic boundary)
    int nextRank = (world_rank + 1) % world_size;
    int prevRank = (world_rank - 1 + world_size) % world_size;

    // Step 1: Exchange the number of boundary particles
    int sendCount = boundaryParticles.size();
    int recvCount = 0;

    MPI_Sendrecv(&sendCount, 1, MPI_INT, nextRank, 0,
                 &recvCount, 1, MPI_INT, prevRank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Step 2: Resize the receivedParticles vector based on the number of particles to be received
    receivedParticles.resize(recvCount);

    // Step 3: Send boundary particles and receive boundary particles from neighboring ranks
    MPI_Sendrecv(boundaryParticles.data(), sendCount * sizeof(Particle), MPI_BYTE, 
                 nextRank, 1,
                 receivedParticles.data(), recvCount * sizeof(Particle), MPI_BYTE, 
                 prevRank, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}


double Simulation::computeBoundaryEnergy(const std::vector<Particle> &receivedParticles) {
    double boundaryEnergy = 0.0;

    // Loop over all local particles and received boundary particles
    for (const auto &localParticle : particles) {
        for (const auto &neighborParticle : receivedParticles) {
            double r2 = box.minimumImageDistanceSquared(localParticle, neighborParticle);
            // std::cout<<"r2"<<r2<<std::endl;
            if (r2 < r2cut) {
                double potential = computePairPotential(r2, potentialType, f_prime, f_d_prime, kappa);
                boundaryEnergy += potential;
            }
        }
    }

    return boundaryEnergy;
}


double Simulation::computeEnergy() {
    // Step 1: Compute local energy for particles in this rank
    double localEnergy = computeLocalEnergy_parallel();
    
    // Step 2: Exchange boundary particles with neighboring rank
    std::vector<Particle> boundaryParticles;  // Identify your boundary particles (e.g., those near boundary)
    std::vector<Particle> receivedParticles;
    // Identify boundary particles
    identifyBoundaryParticles(boundaryParticles);

    exchangeBoundaryParticlesWithNextRank(boundaryParticles, receivedParticles);

    // Step 3: Compute energy between local particles and received boundary particles
    double boundaryEnergy = computeBoundaryEnergy(receivedParticles);
    // std::cout<<boundaryEnergy<<std::endl;
    // Step 4: Sum the energies across all ranks using MPI_Allreduce
    double totalEnergy = 0.0;
    double rankEnergy = localEnergy + boundaryEnergy;
    MPI_Allreduce(&rankEnergy, &totalEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return totalEnergy;
}



/**
 * @brief Runs the simulation for a specified number of steps.
 * 
 * @param numSteps The number of steps to run the simulation.
 * @param equilibrationTime The number of steps for equilibration before logging.
 * @param outputFrequency How often to log the results.
 * @param logger The logging object for output.
 */
void Simulation::run(int numSteps, int equilibrationTime, int outputFrequency, Logging *logger) {
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        int acceptedMoves = 0;
        // energy = computeEnergy();
        
        // // Calculate energy at each step using the appropriate method
        for (int step = 0; step < numSteps + equilibrationTime; ++step) {
        //     if (useCellList){
        //     if (step%cellListUpdateFrequency == 0)
        //     {
        //         buildCellList();
        //         energy = computeEnergy();
        //     }
        //     }
        //     int n_run = particles.size();
        //     for(int i = 0; i<n_run; ++i){
        //         if (monteCarloMove()) {
        //             acceptedMoves++;
        //         }    
                
        //     }

            // Other simulation types can be added here as additional conditions
            
            // Optionally log data
            if (step >= equilibrationTime && step % outputFrequency == 0) {
                // if (fabs(energy - computeEnergy())> 1){
                // std::cerr<<energy<<"   local energy computation is not working.   "<<computeEnergy()<<std::endl;
                // }   
                MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes
                std::vector<Particle> allParticles;  // Will hold all particles on rank 0
                gatherAllParticles(allParticles);


                if (world_rank == 0){
                    logger->logPositions_xyz(allParticles, box, r2cut);
                }
                // logger.logSimulationData(*this, step);

                MPI_Barrier(MPI_COMM_WORLD);  // Ensure all processes finish writing
            }
        }
        std::cout << "Monte Carlo simulation completed with " 
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
                        localEnergy += computePairPotential(r2, potentialType, f_prime, f_d_prime, kappa);
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
                forceSum += computePairForce(r2, potentialType, f_prime, f_d_prime, kappa);
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
    case PotentialType::ThermalStar:{
        double r = sqrt(r2cut);
        double I  =  -exp(1 - r2cut)/4.0 * f_prime - f_d_prime/kappa/kappa * exp(-kappa * r) * (kappa * r + 1);
        answer *= I;
        break;
    }
    case PotentialType::Ideal:
        answer = 0;
        break;
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
    case PotentialType::ThermalStar:{
        double r = sqrt(r2cut);
        double I  =  -exp(1.0 - r2cut) / 2.0 * (1.0 + r2cut) * f_prime + f_d_prime / kappa/kappa * exp(-kappa * r) * (kappa * kappa* r2cut + 2 * kappa * r + 2);
        answer *= I;
        break;
    }
    case PotentialType::Ideal:
        answer = 0;
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



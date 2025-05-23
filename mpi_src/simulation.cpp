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
    : box(box), simtype(simtype), potentialType(potentialType), temperature(temperature), numParticles(numParticles), maxDisplacement(maxDisplacement), rcut(r2cut), r2cut(r2cut * r2cut), f_prime(f_prime), f_d_prime(f_d_prime), kappa(kappa), mu(exp(mu/temperature)), energy(0.0), seed(seed), useCellList(useCellList), cellListUpdateFrequency(cellListUpdateFrequency) , numCellsX(1), numCellsY(1), adjusted_rcut_x(rcut){
    // Initialize MPI variables
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Compute the subdomain size in the x-direction
    double Lx_total = box.getLx();
    double subdomain_x_length = Lx_total / world_size;

    // Adjust the number of cells in the x-direction to fit exactly into the subdomain
    numCellsX = static_cast<int>(subdomain_x_length / rcut);

    numCellsY = static_cast<int>(Lx_total / rcut);
    
    adjusted_rcut_x = subdomain_x_length / numCellsX;  // Adjusted cell length in x-direction
    

    // Ensure that there are enough cells to cover the box
    if (subdomain_x_length > numCellsX * adjusted_rcut_x) numCellsX++;
    if (box.getLy() > numCellsY * rcut) numCellsY++;

    numCellsX += 4; //to fit the boundary particles!

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

    MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes
    buildCellList_parallel();
    
    // Energy computation can be deferred or handled locally

    energy = computeEnergy();
    std::cout<<"totalenergy"<<energy<<std::endl;
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

void Simulation::reassignParticles(){
    // std::cout<<"reassigning"<<std::endl;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    double Lx = box.getLx();
    double invLx = 1.0 / Lx;

    // Compute the subdomain bounds for each rank
    double subdomain_x_min = world_rank * box.getLx() / world_size;
    double subdomain_x_max = (world_rank + 1) * box.getLx() / world_size;


    double domain_length = subdomain_x_max - subdomain_x_min;
    //assigning particles moved to left, assigning particles moved to right:
    std::vector<Particle> particlesToSendLeft, particlesToSendRight;
    std::vector<Particle> particlesToSendLeftBoundary, particlesToSendRightBoundary;
    //Identification of the particles:
    for (auto it = particles.begin(); it != particles.end();){
        //check for moving to the left from the min of the boundary
        double left = it->x - subdomain_x_min;
        left -= Lx * std::round(left * invLx);

        double right = it->x - subdomain_x_max;
        right -= Lx * std::round(right * invLx);

        if (left < 0){
            //checking if it has moved to much 
            if(left< -domain_length)
                std::cerr<<"The Domain-Update time should be increased!"<<std::endl;
            // Particle has crossed to the left subdomain
            // std::cout<<it->x<<" left "<<world_rank<<std::endl;
            particlesToSendLeft.push_back(*it);
            it = particles.erase(it);  // Remove particle from local list
        }

        else if(right > 0){
            if (right > domain_length)
                std::cerr<<"The Domain-Update time should be increased!"<<std::endl;
            // std::cout<<it->x<<" left "<<world_rank<<std::endl;
            particlesToSendRight.push_back(*it);
            it = particles.erase(it);
        }

        else {
            ++it;
        }
        
    }



    int leftNeighbor = (world_rank - 1 + world_size) % world_size;
    int rightNeighbor = (world_rank + 1) % world_size;

    // Send the number of particles to be sent to left and right neighbors
    int sendCountLeft = particlesToSendLeft.size();
    int sendCountRight = particlesToSendRight.size();


    // Receive counts from neighbors
    int recvCountLeft = 0, recvCountRight = 0;
    MPI_Sendrecv(&sendCountLeft, 1, MPI_INT, leftNeighbor, 0,
                 &recvCountRight, 1, MPI_INT, rightNeighbor, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Sendrecv(&sendCountRight, 1, MPI_INT, rightNeighbor, 1,
                 &recvCountLeft, 1, MPI_INT, leftNeighbor, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);


    // Prepare buffers to receive particles
    std::vector<Particle> receivedParticlesLeft(recvCountLeft);
    std::vector<Particle> receivedParticlesRight(recvCountRight);

    // Send particles and receive particles from neighbors
    MPI_Sendrecv(particlesToSendLeft.data(), sendCountLeft * sizeof(Particle), MPI_BYTE,
                 leftNeighbor, 2,
                 receivedParticlesRight.data(), recvCountRight * sizeof(Particle), MPI_BYTE,
                 rightNeighbor, 2,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Sendrecv(particlesToSendRight.data(), sendCountRight * sizeof(Particle), MPI_BYTE,
                 rightNeighbor, 3,
                 receivedParticlesLeft.data(), recvCountLeft * sizeof(Particle), MPI_BYTE,
                 leftNeighbor, 3,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Add received particles to the local particle list
    particles.insert(particles.end(), receivedParticlesLeft.begin(), receivedParticlesLeft.end());
    particles.insert(particles.end(), receivedParticlesRight.begin(), receivedParticlesRight.end());

    MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes
    buildCellList_parallel();
}


bool Simulation::monteCarloMove_parallel(int id) {
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    if (world_rank % 2 == id){
        if (particles.size() != 0){
            // Randomly select a particle
            size_t particleIndex = rand() % particles.size();
            Particle &p = particles[particleIndex];

            // Store the old position
            double oldX = p.x;
            double oldY = p.y;
            box.applyPBC(p);
            // std::cout<<"Montecarlo Move1:"<<computeLocalCellsEnergy_parallel(particleIndex)<<std::endl;
            double pre_update_energy = computeLocalCellsEnergy_parallel(particleIndex);

            // Propose a random displacement
            double dx = (rand() / double(RAND_MAX) - 0.5) * maxDisplacement;
            double dy = (rand() / double(RAND_MAX) - 0.5) * maxDisplacement;

            // Update the particle's position
            p.x += dx;
            p.y += dy;
            box.applyPBC(p);

            double post_update_energy = computeLocalCellsEnergy_parallel(particleIndex);

            double delta_E = post_update_energy - pre_update_energy;

            if (delta_E<0)
            {
                //aceepting:
                // energy += delta_E;
                return true;

            }
            else if(exp(-delta_E / temperature) > (rand() / double(RAND_MAX)))
            {
                // energy += delta_E;
                return true;
            }
            else
                {
                    p.x -= dx;
                    p.y -= dy;
                    box.applyPBC(p);

                }

                return false;
        }
    }
            
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
        if (1 == 1) {
            // If the particle is near the left or right boundary, add it to the boundaryParticles list
            boundaryParticles.push_back(particle);

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

    // std::cout<<sendCount<<std::endl;
    MPI_Sendrecv(&sendCount, 1, MPI_INT, nextRank, 0,
                 &recvCount, 1, MPI_INT, prevRank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Step 2: Resize the receivedParticles vector based on the number of particles to be received
    receivedParticles.resize(recvCount);
    // std::cout<<recvCount<<std::endl;
    // Step 3: Send boundary particles and receive boundary particles from neighboring ranks
    MPI_Sendrecv(boundaryParticles.data(), sendCount * sizeof(Particle), MPI_BYTE, 
                 nextRank, 1,
                 receivedParticles.data(), recvCount * sizeof(Particle), MPI_BYTE, 
                 prevRank, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}


double Simulation::computeBoundaryEnergy(const std::vector<Particle> &receivedParticles) {
    double boundaryEnergy = 0.0;
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Loop over all local particles and received boundary particles

    for (const auto &localParticle : particles) {
        
        for (const auto &neighborParticle : receivedParticles) {
            // std::cout<<world_rank<<" "<<localParticle.x<<" "<<neighborParticle.x<<std::endl;
            double r2 = box.minimumImageDistanceSquared(localParticle, neighborParticle);
            if (r2 < r2cut) {
                double potential = computePairPotential(r2, potentialType, f_prime, f_d_prime, kappa);
                boundaryEnergy += potential;
            }
            // std::cout<<boundaryEnergy<<std::endl;
        }
    }

    return boundaryEnergy;
}


double Simulation::computeEnergy() {
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Step 1: Compute local energy for particles in this rank
    double localEnergy = computeLocalEnergy_parallel();
    // std::cout<<world_rank<<" "<<localEnergy<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);  // Ensure all processes are synched
    // std::cout<<localEnergy<<std::endl;
    // Step 2: Exchange boundary particles with neighboring rank
    std::vector<Particle> boundaryParticles;  // Identify your boundary particles (e.g., those near boundary)
    std::vector<Particle> receivedParticles;
    // Identify boundary particles
    identifyBoundaryParticles(boundaryParticles);
    
    MPI_Barrier(MPI_COMM_WORLD);  // Ensure all processes are synched
    exchangeBoundaryParticlesWithNextRank(boundaryParticles, receivedParticles);
    MPI_Barrier(MPI_COMM_WORLD);  // Ensure all processes are synched
    // Step 3: Compute energy between local particles and received boundary particles
    double boundaryEnergy = computeBoundaryEnergy(receivedParticles);

    MPI_Barrier(MPI_COMM_WORLD);  // Ensure all processes are synched
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
        float acceptedMoves = 0;
        // energy = computeEnergy();
        
        // // Calculate energy at each step using the appropriate method
        for (int step = 0; step < numSteps + equilibrationTime; ++step) {
            for (int tmp = 0 ; tmp<numParticles; tmp++){
                acceptedMoves += monteCarloMove_parallel(0);
                MPI_Barrier(MPI_COMM_WORLD); 
                acceptedMoves += monteCarloMove_parallel(1);
                MPI_Barrier(MPI_COMM_WORLD); 
            }
            if (step%cellListUpdateFrequency == 0)
            {
                reassignParticles();
                MPI_Barrier(MPI_COMM_WORLD);
            }


            // Other simulation types can be added here as additional conditions
            
            // Optionally log data
            if (step % outputFrequency == 0) {
                // if (fabs(energy - computeEnergy())> 1){
                // std::cerr<<energy<<"   local energy computation is not working.   "<<computeEnergy()<<std::endl;
                // }   
                reassignParticles();
                MPI_Barrier(MPI_COMM_WORLD);
                energy = computeEnergy();
                std::cout<<energy<<std::endl;
                MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes
                std::vector<Particle> allParticles;  // Will hold all particles on rank 0
                gatherAllParticles(allParticles);


                if (world_rank == 0){
                    logger->logPositions_xyz(allParticles, box, r2cut);
                }
                // logger.logSimulationData(*this, step);

                MPI_Barrier(MPI_COMM_WORLD);  // Ensure all processes are synched
            }
        }
        std::cout << "Monte Carlo simulation completed with " 
                    <<  " accepted moves out of " << (acceptedMoves/(numSteps + equilibrationTime)) << " steps." << std::endl;
        
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




 void Simulation::buildCellList_parallel() {
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int leftNeighbor = (world_rank - 1 + world_size) % world_size;
    int rightNeighbor = (world_rank + 1) % world_size;

    double Lx_total = box.getLx();
    double invLx = 1/Lx_total;
    double subdomain_x_length = Lx_total / world_size;

    std::vector<Particle> sentParticlesLeft, sentParticlesRight;
    sentParticlesLeft.clear();
    sentParticlesRight.clear();

    clearCellList();
    cellList.resize((numCellsX) * numCellsY, nullptr);

    for (size_t i = 0; i < particles.size(); ++i) {
        double local_x = particles[i].x - world_rank * subdomain_x_length;
        int cellX = static_cast<int>(local_x / adjusted_rcut_x);
        int cellY = static_cast<int>(particles[i].y / rcut);

        int cellIndex = cellY * (numCellsX) + (cellX + 2);
        assert(cellIndex >= 0 && cellIndex < cellList.size());

        CellListNode* newNode = new CellListNode(i);
        newNode->next = cellList[cellIndex];
        cellList[cellIndex] = newNode;

        if (cellX < 2) {
            // std::cout<<"sending a particle to the left boundary:"<<local_x<<std::endl;
            sentParticlesLeft.push_back(particles[i]);  // Near left boundary
        } else if (cellX > numCellsX - 7) {
            //  std::cout<<"sending a particle to the right boundary:"<<local_x<<std::endl;
            sentParticlesRight.push_back(particles[i]);  // Near right boundary
        }
    }

    boundaryParticlesLeft.clear();
    boundaryParticlesRight.clear();

    int sendCountLeft = sentParticlesLeft.size();
    int sendCountRight = sentParticlesRight.size();
    int recvCountLeft = 0, recvCountRight = 0;

    MPI_Sendrecv(&sendCountLeft, 1, MPI_INT, leftNeighbor, 0,
                 &recvCountRight, 1, MPI_INT, rightNeighbor, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Sendrecv(&sendCountRight, 1, MPI_INT, rightNeighbor, 1,
                 &recvCountLeft, 1, MPI_INT, leftNeighbor, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    boundaryParticlesLeft.resize(recvCountLeft);
    boundaryParticlesRight.resize(recvCountRight);

    MPI_Sendrecv(sentParticlesLeft.data(), sendCountLeft * sizeof(Particle), MPI_BYTE,
                 leftNeighbor, 2,
                 boundaryParticlesRight.data(), recvCountRight * sizeof(Particle), MPI_BYTE,
                 rightNeighbor, 2,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Sendrecv(sentParticlesRight.data(), sendCountRight * sizeof(Particle), MPI_BYTE,
                 rightNeighbor, 3,
                 boundaryParticlesLeft.data(), recvCountLeft * sizeof(Particle), MPI_BYTE,
                 leftNeighbor, 3,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);


    for (size_t i = 0; i < boundaryParticlesLeft.size(); ++i) {
        double local_x = boundaryParticlesLeft[i].x - world_rank * subdomain_x_length;
        //applying rhe local X boudary condtion is missing!
        local_x -= Lx_total * std::round(local_x * invLx);
        int cellX = static_cast<int>(local_x / adjusted_rcut_x) - 1;
        // std::cout<<"cellX in Left boundary particles: "<<cellX<<std::endl;
        int cellY = static_cast<int>(boundaryParticlesLeft[i].y / rcut);

        int cellIndex = cellY * (numCellsX) + (cellX + 2);
        assert(cellIndex >= 0 && cellIndex < cellList.size());
        
        CellListNode* newNode = new CellListNode((-2 * i - 2));
        // leftBoundaryCounter++;
        newNode->next = cellList[cellIndex];
        cellList[cellIndex] = newNode;

        // std::cout << "Left boundaryRank " << world_rank << " - Particle " << i
        //   << " at (" << boundaryParticlesLeft[i].x << ", " << boundaryParticlesLeft[i].y
        //   << ") assigned to Cell (" << cellX << ", " << cellY << ") with CellIndex " << cellIndex << "\n";


    }

    for (size_t i = 0; i < boundaryParticlesRight.size(); ++i) {
        
        double local_x = boundaryParticlesRight[i].x - world_rank * subdomain_x_length;
        //applying rhe local X boudary condtion is missing!
        local_x -= Lx_total * std::round(local_x * invLx);
        int cellX = static_cast<int>(local_x / adjusted_rcut_x);
        int cellY = static_cast<int>(boundaryParticlesRight[i].y / rcut);
        // std::cout<<"cellX in right boundary particles: "<<cellX<<std::endl;
        int cellIndex = cellY * (numCellsX) + (cellX + 2);
        assert(cellIndex >= 0 && cellIndex < cellList.size());

        CellListNode* newNode = new CellListNode((-2 * i - 1));
        // rightBoundaryCounter++;
        newNode->next = cellList[cellIndex];
        cellList[cellIndex] = newNode;




    }
}


double Simulation::computeLocalCellsEnergy_parallel(int particleIndex) const {
    double localEnergy = 0;
    const Particle& p = particles[particleIndex];

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Compute the subdomain size in the x-direction
    double Lx_total = box.getLx();
    double invLx = 1/Lx_total;
    double subdomain_x_length = Lx_total / world_size;

    double local_x = p.x - world_rank * subdomain_x_length;
    local_x -= Lx_total * std::round(local_x * invLx);
    int cellX = static_cast<int>(local_x / adjusted_rcut_x) + 2;
    int cellY = static_cast<int>(p.y / rcut);

    // for the case of the cells does not need to make contact with cells of other domains:
    if ((cellX>3) && (cellX < numCellsX - 4)){
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
    }
    
    else if (cellX == 3){
        for (int offsetY = -2; offsetY <= 2; ++offsetY) {
            int neighborCellY = (cellY + offsetY + numCellsY) % numCellsY;
            for (int offsetX = -1; offsetX <= 2; ++offsetX) {
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

            //doing offset -2 by hand:
            int offsetX = -2;
            int neighborCellX = (cellX + offsetX + numCellsX) % numCellsX;
            int neighborCellIndex = neighborCellY * numCellsX + neighborCellX;


            
            
            for (CellListNode* node = cellList[neighborCellIndex]; node != nullptr; node = node->next) {

                int computed_node = static_cast<int>(-((node->particleIndex + 2) / 2));
                
                double r2 = box.minimumImageDistanceSquared(p, boundaryParticlesLeft[computed_node]);
                if (r2 < r2cut) {
                    localEnergy += computePairPotential(r2, potentialType, f_prime, f_d_prime, kappa);
                    }

                }

        }


    }

        else if (cellX == 2){
        for (int offsetY = -2; offsetY <= 2; ++offsetY) {
            int neighborCellY = (cellY + offsetY + numCellsY) % numCellsY;
            for (int offsetX = 0; offsetX <= 2; ++offsetX) {
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

            //doing offset -2, -1 by hand:
            for (int offsetX = -2; offsetX <= -1; ++offsetX) {
                int neighborCellX = (cellX + offsetX + numCellsX) % numCellsX;
                int neighborCellIndex = neighborCellY * numCellsX + neighborCellX;

                for (CellListNode* node = cellList[neighborCellIndex]; node != nullptr; node = node->next) {
                    // std::cout<<"mpi node of neighbors"<<node->particleIndex<<std::endl;
                    int computed_node = static_cast<int>(-((node->particleIndex + 2) / 2));
                    // std::cout<<"computed node of neighbors"<<computed_node<<std::endl;
                    double r2 = box.minimumImageDistanceSquared(p, boundaryParticlesLeft[computed_node]);
                    if (r2 < r2cut) {
                        localEnergy += computePairPotential(r2, potentialType, f_prime, f_d_prime, kappa);
                        }

            //     std::cout << "Rank " << world_rank << " - Computing local energy for Particle " << particleIndex 
            //   << " at (" << p.x << ", " << p.y << ") assigned to Cell (" << cellX << ", " << cellY << ")\n";
            //     std::cout << "Checking boundary right, Particle index " << computed_node 
            //     << " position: (" << boundaryParticlesLeft[computed_node].x << ", "
            //     << boundaryParticlesLeft[computed_node].y << "), r^2: " << r2 <<" neighbor cell x and y"<<neighborCellX<<" "<<neighborCellY<<"\n";
                
                        
                    }
                    
            }

        }


    }


        else if (cellX == 1){
        for (int offsetY = -2; offsetY <= 2; ++offsetY) {
            int neighborCellY = (cellY + offsetY + numCellsY) % numCellsY;
            for (int offsetX = 1; offsetX <= 2; ++offsetX) {
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

            //doing offset -1, 0 by hand:
            for (int offsetX = -1; offsetX <= 0; ++offsetX) {
                int neighborCellX = (cellX + offsetX + numCellsX) % numCellsX;
                int neighborCellIndex = neighborCellY * numCellsX + neighborCellX;

                for (CellListNode* node = cellList[neighborCellIndex]; node != nullptr; node = node->next) {
                    // std::cout<<"mpi node of neighbors"<<node->particleIndex<<std::endl;
                    int computed_node = static_cast<int>(-((node->particleIndex + 2) / 2));
                    // std::cout<<"computed node of neighbors"<<computed_node<<std::endl;
                    double r2 = box.minimumImageDistanceSquared(p, boundaryParticlesLeft[computed_node]);
                    if (r2 < r2cut) {
                        localEnergy += computePairPotential(r2, potentialType, f_prime, f_d_prime, kappa);
                        }

                    }
            }

        }


    }
    
    else if (cellX == numCellsX - 4){
        for (int offsetY = -2; offsetY <= 2; ++offsetY) {
            int neighborCellY = (cellY + offsetY + numCellsY) % numCellsY;
            for (int offsetX = -2; offsetX <= 1; ++offsetX) {
                int neighborCellX = (cellX + offsetX + numCellsX) % numCellsX;
                int neighborCellIndex = neighborCellY * (numCellsX) + neighborCellX;

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

            //doing offset 2 by hand:
            int offsetX = 2;
            int neighborCellX = (cellX + offsetX + numCellsX) % numCellsX;
            int neighborCellIndex = neighborCellY * numCellsX + neighborCellX;

            for (CellListNode* node = cellList[neighborCellIndex]; node != nullptr; node = node->next) {
                // std::cout<<"mpi node of neighbors"<<node->particleIndex<<std::endl;
                int computed_node = static_cast<int>(-((node->particleIndex + 1) / 2));
                // std::cout<<"computed node of neighbors"<<computed_node<<std::endl;
                double r2 = box.minimumImageDistanceSquared(p, boundaryParticlesRight[computed_node]);
                if (r2 < r2cut) {
                    localEnergy += computePairPotential(r2, potentialType, f_prime, f_d_prime, kappa);
                    }
                }

        }


    }

        else if (cellX == numCellsX - 3){
        for (int offsetY = -2; offsetY <= 2; ++offsetY) {
            int neighborCellY = (cellY + offsetY + numCellsY) % numCellsY;
    // Inside this loop, check individual particle positions
            for (int offsetX = -2; offsetX <= 0; ++offsetX) {
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

            //doing offset 1, 2 by hand:
            for (int offsetX = 1; offsetX <= 2; ++offsetX) {
                int neighborCellX = (cellX + offsetX + numCellsX) % numCellsX;
                int neighborCellIndex = (neighborCellY) * numCellsX + neighborCellX;

                for (CellListNode* node = cellList[neighborCellIndex]; node != nullptr; node = node->next) {
                    // std::cout<<"mpi node of neighbors"<<node->particleIndex<<std::endl;
                    int computed_node = static_cast<int>(-((node->particleIndex + 1) / 2));
                    // std::cout<<"computed node of neighbors"<<computed_node<<std::endl;
                    double r2 = box.minimumImageDistanceSquared(p, boundaryParticlesRight[computed_node]);
                    if (r2 < r2cut) {
                        localEnergy += computePairPotential(r2, potentialType, f_prime, f_d_prime, kappa);
                        }

            //     std::cout << "Rank " << world_rank << " - Computing local energy for Particle " << particleIndex 
            //   << " at (" << p.x << ", " << p.y << ") assigned to Cell (" << cellX << ", " << cellY << ")\n";
            //     std::cout << "Checking boundary right, Particle index " << computed_node 
            //     << " position: (" << boundaryParticlesRight[computed_node].x << ", "
            //     << boundaryParticlesRight[computed_node].y << "), r^2: " << r2 <<" neighbor cell x and y"<<neighborCellX<<" "<<neighborCellY<<"\n";
                    }
            }

        }


    }


        else if (cellX == numCellsX - 2){
        for (int offsetY = -2; offsetY <= 2; ++offsetY) {
            int neighborCellY = (cellY + offsetY + numCellsY) % numCellsY;
            for (int offsetX = -2; offsetX <= -1; ++offsetX) {
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

            //doing offset 1, 0 by hand:
            for (int offsetX = 0; offsetX <= 1; ++offsetX) {
                int neighborCellX = (cellX + offsetX + numCellsX) % numCellsX;
                int neighborCellIndex = neighborCellY * numCellsX + neighborCellX;

                for (CellListNode* node = cellList[neighborCellIndex]; node != nullptr; node = node->next) {
                    // std::cout<<"mpi node of neighbors"<<node->particleIndex<<std::endl;
                    int computed_node = static_cast<int>(-((node->particleIndex + 1) / 2));
                    // std::cout<<"computed node of neighbors"<<computed_node<<std::endl;
                    double r2 = box.minimumImageDistanceSquared(p, boundaryParticlesRight[computed_node]);
                    if (r2 < r2cut) {
                        localEnergy += computePairPotential(r2, potentialType, f_prime, f_d_prime, kappa);
                        }
                    }
            }

        }


    }
    else{
        std::cout<<"cellx"<<cellX<<std::endl;
        std::cerr<<"Particles are outside the decomposed CPU area! Please Repeat the simulation using higher decomposition reassignment rate."<<std::endl;
        // MPI_Finalize();
        // exit(EXIT_SUCCESS);  
    }


    return localEnergy;
}

SimulationType selectSimulationType(const std::string &simulationName) {
    if (simulationName == "NVT") return SimulationType::MonteCarloNVT;
    if (simulationName == "GCMC") return SimulationType::GCMC;
    throw std::invalid_argument("Unknown simulation type: " + simulationName);
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



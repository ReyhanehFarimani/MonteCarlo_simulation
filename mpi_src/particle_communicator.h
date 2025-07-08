#ifndef PARTICLE_COMMUNICATOR_H
#define PARTICLE_COMMUNICATOR_H

#include "initial.h"
#include "domain_decomposition.h"
#include <vector>
#include <mpi.h>

/**
 * @brief Handles ghost particle exchange and particle migration between MPI ranks.
 */
class ParticleCommunicator {
public:
    /**
     * @brief Constructs the communicator with MPI context and cutoff radius.
     * @param comm The MPI communicator.
     * @param rcut The interaction cutoff radius.
     */
    ParticleCommunicator(MPI_Comm comm, double rcut);

    /**
     * @brief Exchange ghost particles across all 8 neighbors.
     * Called after a cell list rebuild.
     * @param dd Domain decomposition helper.
     * @param localParticles Particles owned by this rank.
     * @param ghostParticles Will be filled with ghost particles from neighbors.
     */
    void exchangeGhostParticles(
        const DomainDecomposition& dd,
        const std::vector<Particle>& localParticles,
        std::vector<Particle>& ghostParticles
    );

    /**
     * @brief Migrate particles that left the current rank's domain.
     * Called after MC moves.
     * @param dd Domain decomposition helper.
     * @param localParticles Input/output: particles will be updated after migration.
     */
    void migrateParticles(
        const DomainDecomposition& dd,
        std::vector<Particle>& localParticles
    );

private:
    MPI_Comm comm_;
    double rcut_;
    int rank_, size_;

    /**
     * @brief Pack edge-layer particles for ghost exchange per neighbor.
     */
    void packParticlesForGhostExchange(
        const DomainDecomposition& dd,
        const std::vector<Particle>& localParticles,
        std::vector<std::vector<Particle>>& sendBuffers
    );

    /**
     * @brief Unpack received particles into the ghost vector.
     */
    void unpackParticlesToGhostVector(
        const std::vector<std::vector<Particle>>& recvBuffers,
        std::vector<Particle>& ghostParticles
    );

    /**
     * @brief Perform all non-blocking MPI send/recv operations with neighbors.
     */
    void sendAndReceive(
        const std::vector<std::vector<Particle>>& sendBuffers,
        std::vector<std::vector<Particle>>& recvBuffers,
        const std::vector<int>& neighborRanks
    );
};

#endif


#include "particle_communicator.h"
#include <mpi.h>
#include <algorithm>

ParticleCommunicator::ParticleCommunicator(MPI_Comm comm, double rcut)
    : comm_(comm), rcut_(rcut)
{
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);
}

void ParticleCommunicator::exchangeGhostParticles(
    const DomainDecomposition& dd,
    const std::vector<Particle>& localParticles,
    std::vector<Particle>& ghostParticles)
{
    std::vector<int> neighborRanks = dd.getNeighborRanks();
    const int nNeighbors = static_cast<int>(neighborRanks.size());

    std::vector<std::vector<Particle>> sendBuffers(nNeighbors);
    std::vector<std::vector<Particle>> recvBuffers(nNeighbors);

    packParticlesForGhostExchange(dd, localParticles, sendBuffers);
    sendAndReceive(sendBuffers, recvBuffers, neighborRanks);
    unpackParticlesToGhostVector(recvBuffers, ghostParticles);
}

void ParticleCommunicator::packParticlesForGhostExchange(
    const DomainDecomposition& dd,
    const std::vector<Particle>& localParticles,
    std::vector<std::vector<Particle>>& sendBuffers)
{
    const int dims[2] = {dd.getPx(), dd.getPy()};
    const double lx = dd.getLocalBox().getLx();
    const double ly = dd.getLocalBox().getLy();
    const double ox = dd.getOriginX();
    const double oy = dd.getOriginY();

    const int dx[8] = {-1, 1,  0,  0, -1, 1, -1, 1};
    const int dy[8] = { 0, 0, -1,  1, -1, -1,  1, 1};

    for (size_t dir = 0; dir < 8; ++dir) {
        for (const Particle& p : localParticles) {
            bool send = false;
            if (dx[dir] == -1 && p.x < ox + rcut_) send = true;
            if (dx[dir] ==  1 && p.x > ox + lx - rcut_) send = true;
            if (dy[dir] == -1 && p.y < oy + rcut_) send = true;
            if (dy[dir] ==  1 && p.y > oy + ly - rcut_) send = true;

            if ((dx[dir] != 0 && dy[dir] != 0)) {
                bool inX = (dx[dir] == -1 && p.x < ox + rcut_) || (dx[dir] == 1 && p.x > ox + lx - rcut_);
                bool inY = (dy[dir] == -1 && p.y < oy + rcut_) || (dy[dir] == 1 && p.y > oy + ly - rcut_);
                send = inX && inY;
            }

            if (send) {
                sendBuffers[dir].push_back(p);
            }
        }
    }
}

void ParticleCommunicator::flatten(const std::vector<Particle>& in, std::vector<double>& out)
{
    out.resize(3 * in.size());
    for (size_t i = 0; i < in.size(); ++i) {
        out[3 * i + 0] = in[i].x;
        out[3 * i + 1] = in[i].y;
        out[3 * i + 2] = static_cast<double>(in[i].id);
    }
}

void ParticleCommunicator::inflate(const std::vector<double>& in, std::vector<Particle>& out)
{
    size_t count = in.size() / 3;
    out.resize(count);
    for (size_t i = 0; i < count; ++i) {
        out[i].x = in[3 * i + 0];
        out[i].y = in[3 * i + 1];
        out[i].id = static_cast<int>(in[3 * i + 2]);
    }
}

void ParticleCommunicator::sendAndReceive(
    const std::vector<std::vector<Particle>>& sendBuffers,
    std::vector<std::vector<Particle>>& recvBuffers,
    const std::vector<int>& neighborRanks)
{
    const int nNeighbors = static_cast<int>(neighborRanks.size());
    std::vector<MPI_Request> requests(2 * nNeighbors);

    std::vector<int> sendCounts(nNeighbors);
    std::vector<int> recvCounts(nNeighbors);
    std::vector<MPI_Request> countRequests(2 * nNeighbors);

    // Exchange counts first
    for (int i = 0; i < nNeighbors; ++i) {
        sendCounts[i] = static_cast<int>(sendBuffers[i].size());
        MPI_Irecv(&recvCounts[i], 1, MPI_INT, neighborRanks[i], 99, comm_, &countRequests[i]);
        MPI_Isend(&sendCounts[i], 1, MPI_INT, neighborRanks[i], 99, comm_, &countRequests[nNeighbors + i]);
    }
    MPI_Waitall(2 * nNeighbors, countRequests.data(), MPI_STATUSES_IGNORE);

    // Flatten data to double buffers
    std::vector<std::vector<double>> sendFlat(nNeighbors), recvFlat(nNeighbors);
    for (int i = 0; i < nNeighbors; ++i) {
        flatten(sendBuffers[i], sendFlat[i]);
    }

    // Post receives
    for (int i = 0; i < nNeighbors; ++i) {
        if (recvCounts[i] > 0) {
            recvFlat[i].resize(3 * recvCounts[i]);
            MPI_Irecv(recvFlat[i].data(), 3 * recvCounts[i], MPI_DOUBLE, neighborRanks[i], 0, comm_, &requests[i]);
        } else {
            requests[i] = MPI_REQUEST_NULL;
        }
    }

    // Post sends
    for (int i = 0; i < nNeighbors; ++i) {
        if (!sendFlat[i].empty()) {
            MPI_Isend(sendFlat[i].data(), static_cast<int>(sendFlat[i].size()), MPI_DOUBLE, neighborRanks[i], 0, comm_, &requests[nNeighbors + i]);
        } else {
            requests[nNeighbors + i] = MPI_REQUEST_NULL;
        }
    }

    MPI_Waitall(2 * nNeighbors, requests.data(), MPI_STATUSES_IGNORE);

    // Inflate back
    for (int i = 0; i < nNeighbors; ++i) {
        inflate(recvFlat[i], recvBuffers[i]);
    }
}

void ParticleCommunicator::unpackParticlesToGhostVector(
    const std::vector<std::vector<Particle>>& recvBuffers,
    std::vector<Particle>& ghostParticles)
{
    ghostParticles.clear();
    for (const auto& buf : recvBuffers) {
        ghostParticles.insert(ghostParticles.end(), buf.begin(), buf.end());
    }
}

void ParticleCommunicator::migrateParticles(
    const DomainDecomposition& dd,
    std::vector<Particle>& localParticles)
{
    // Placeholder: implement migration after displacement step
    // Check which particles moved outside the local domain and forward them
}
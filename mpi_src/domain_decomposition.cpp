#include "domain_decomposition.h"
#include <cmath>
#include <stdexcept>

DomainDecomposition::DomainDecomposition(const SimulationBox& globalBox, double rcut, MPI_Comm comm)
    : globalBox_(globalBox),
      rcut_(rcut),
      comm_(comm),
      localBox_(1.0, 1.0)  // Temporary, will be overwritten below
{
    // Get basic MPI info (from original communicator)
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);

    // Determine 2D grid layout (Px × Py) as square-ish as possible
    Px_ = static_cast<int>(std::sqrt(size_));
    while (size_ % Px_ != 0) --Px_;
    Py_ = size_ / Px_;

    if (Px_ * Py_ != size_) {
        throw std::runtime_error("Could not find valid Px × Py decomposition");
    }

    // Set up periodic Cartesian topology
    int dims[2] = {Px_, Py_};
    int periods[2] = {1, 1};  // PBC in both directions
    int reorder = 0;
    MPI_Cart_create(comm_, 2, dims, periods, reorder, &cart_comm_);

    // Get this rank's coordinates in the Cartesian grid
    int cart_rank;
    MPI_Comm_rank(cart_comm_, &cart_rank);  // Rank within Cartesian communicator

    int coords[2];
    MPI_Cart_coords(cart_comm_, cart_rank, 2, coords);
    coordX_ = coords[0];
    coordY_ = coords[1];

    // Assign subdomain dimensions
    double lx = globalBox_.getLx() / Px_;
    double ly = globalBox_.getLy() / Py_;
    localBox_ = SimulationBox(lx, ly);

    // Compute global offset of this domain
    originX_ = coordX_ * lx;
    originY_ = coordY_ * ly;
}

const SimulationBox& DomainDecomposition::getLocalBox() const {
    return localBox_;
}

int DomainDecomposition::getRank() const {
    return rank_;
}

int DomainDecomposition::getSize() const {
    return size_;
}

int DomainDecomposition::getPx() const {
    return Px_;
}

int DomainDecomposition::getPy() const {
    return Py_;
}

int DomainDecomposition::getCoordX() const {
    return coordX_;
}

int DomainDecomposition::getCoordY() const {
    return coordY_;
}

double DomainDecomposition::getOriginX() const {
    return originX_;
}

double DomainDecomposition::getOriginY() const {
    return originY_;
}

std::pair<double, double> DomainDecomposition::getGlobalPosition(const Particle& p) const {
    double gx = originX_ + p.x;
    double gy = originY_ + p.y;
    return {gx, gy};
}

MPI_Comm DomainDecomposition::getCartComm() const {
    return cart_comm_;
}

std::vector<int> DomainDecomposition::getNeighborRanks() const {
    std::vector<int> neighbors;

    const int directions[8][2] = {
        {-1,  0}, {1,  0},  // W, E
        { 0, -1}, {0,  1},  // S, N
        {-1, -1}, {1, -1},  // SW, SE
        {-1,  1}, {1,  1}   // NW, NE
    };

    int dims[2] = {Px_, Py_};

    for (const auto& dir : directions) {
        int coords[2] = {
            (coordX_ + dir[0] + dims[0]) % dims[0],
            (coordY_ + dir[1] + dims[1]) % dims[1]
        };

        int neighbor_rank;
        MPI_Cart_rank(cart_comm_, coords, &neighbor_rank);
        neighbors.push_back(neighbor_rank);
    }

    return neighbors;
}

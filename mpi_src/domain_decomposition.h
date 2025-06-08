#ifndef DOMAIN_DECOMPOSITION_H
#define DOMAIN_DECOMPOSITION_H

#include "initial.h"   // For SimulationBox
#include <vector>
#include <mpi.h>

/**
 * @brief Handles domain splitting and neighbor management in a 2D MPI grid.
 */
class DomainDecomposition {
public:
    DomainDecomposition(const SimulationBox& globalBox, double rcut, MPI_Comm comm);

    const SimulationBox& getLocalBox() const;

    std::vector<int> getNeighborRanks() const;

    int getRank() const;
    int getSize() const;

    int getPx() const;
    int getPy() const;
    int getCoordX() const;
    int getCoordY() const;

    MPI_Comm getCartComm() const;  // Getter for the Cartesian communicator

private:
    SimulationBox globalBox_;
    SimulationBox localBox_;

    double rcut_;
    int rank_, size_;
    int Px_, Py_;
    int coordX_, coordY_;
    MPI_Comm comm_;        // Original communicator (MPI_COMM_WORLD)
    MPI_Comm cart_comm_;   // New Cartesian communicator (returned by MPI_Cart_create)
};

#endif

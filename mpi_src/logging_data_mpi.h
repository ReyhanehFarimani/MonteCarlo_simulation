#ifndef LOGGING_DATA_MPI_H
#define LOGGING_DATA_MPI_H

#include <string>
#include <fstream>
#include <mpi.h>
#include "particle.h"
#include "simulation_box.h"
#include "cell_list_parallel.h"
#include "particle_exchange.h"
#include "thermodynamic_calculator_parallel.h"

/**
 * @brief Parallel logger for thermodynamic scalars (CSV), rank 0 writer.
 *
 * Columns: timestep, N_global, rho, T, U, W, P, V
 */
class LoggingDataMPI {
public:
    LoggingDataMPI(const std::string& basename,
                   MPI_Comm comm = MPI_COMM_WORLD,
                   bool append = false);

    ~LoggingDataMPI();

    /** Write one CSV row of global thermodynamic data (rank 0 only). */
    void log_step(const std::vector<Particle>& owned,
                  const SimulationBox& box,
                  const CellListParallel& clp,
                  const ParticleExchange& pex,
                  const ThermodynamicCalculatorParallel& thermo,
                  int timestep);

    void close();

private:
    void ensure_open_();

    std::string fname_;
    MPI_Comm    comm_;
    int         rank_{0};
    bool        append_{false};
    bool        header_written_{false};

    std::ofstream ofs_;
};

#endif // LOGGING_DATA_MPI_H

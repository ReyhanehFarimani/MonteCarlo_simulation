#include "logging_data_mpi.h"
#include <iomanip>
#include <stdexcept>

LoggingDataMPI::LoggingDataMPI(const std::string& basename, MPI_Comm comm, bool append)
    : fname_(basename + ".thermo.csv"), comm_(comm), append_(append)
{
    MPI_Comm_rank(comm_, &rank_);
}

LoggingDataMPI::~LoggingDataMPI() { close(); }

void LoggingDataMPI::ensure_open_() {
    if (rank_ != 0) return;
    if (ofs_.is_open()) return;

    std::ios_base::openmode mode = std::ios::out;
    mode |= (append_ ? std::ios::app : std::ios::trunc);
    ofs_.open(fname_, mode);
    if (!ofs_) throw std::runtime_error("LoggingDataMPI: cannot open file: " + fname_);

    if (!append_) header_written_ = false;
}

void LoggingDataMPI::close() {
    if (rank_ == 0 && ofs_.is_open()) {
        ofs_.flush();
        ofs_.close();
    }
}

void LoggingDataMPI::log_step(const std::vector<Particle>& owned,
                              const SimulationBox& box,
                              const CellListParallel& clp,
                              const ParticleExchange& pex,
                              const ThermodynamicCalculatorParallel& thermo,
                              int timestep)
{
    // Global particle count
    int nloc = static_cast<int>(owned.size());
    int nglob = 0;
    MPI_Allreduce(&nloc, &nglob, 1, MPI_INT, MPI_SUM, comm_);

    // Compute thermodynamic quantities (these already do the needed MPI reductions internally)
    const double U = thermo.totalEnergy(owned, box, clp, pex);
    const double W = thermo.totalVirial(owned, box, clp, pex);
    const double P = thermo.pressure(owned, box, clp, pex);
    const double T = thermo.getTemperature();
    const double V = box.getV();
    const double rho = (V > 0.0) ? (static_cast<double>(nglob) / V) : 0.0;

    // Rank 0 writes a CSV row
    if (rank_ == 0) {
        ensure_open_();
        if (!header_written_) {
            ofs_ << "timestep,N_global,rho,T,U,W,P,V\n";
            header_written_ = true;
        }
        ofs_ << std::setprecision(17)
             << timestep << ","
             << nglob    << ","
             << rho      << ","
             << T        << ","
             << U        << ","
             << W        << ","
             << P        << ","
             << V        << "\n";
        ofs_.flush();
    }
}

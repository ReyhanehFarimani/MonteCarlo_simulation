#ifndef LOGGING_TRAJ_MPI_H
#define LOGGING_TRAJ_MPI_H

#include <string>
#include <vector>
#include <fstream>
#include <mpi.h>
#include "particle.h"
#include "simulation_box.h"

class LoggingTrajMPI {
public:
    enum class Mode { PerRank, GatherRoot, MPIIO };

    LoggingTrajMPI(std::string basename,
                   Mode mode,
                   MPI_Comm comm = MPI_COMM_WORLD,
                   bool append = false,
                   bool rank_as_type = false);

    ~LoggingTrajMPI();

    // Append one XYZ frame (2D; z=0). Positions are GLOBAL.
    void log_xyz(const std::vector<Particle>& local, const SimulationBox& box);

    // Append one LAMMPS dump frame (2D; z=0). Positions are GLOBAL.
    void log_dump(const std::vector<Particle>& local, const SimulationBox& box, int timestep);

    void close();

private:
    // Common
    std::string basename_;
    Mode        mode_;
    MPI_Comm    comm_;
    int         rank_{0}, size_{1};
    bool        rank_as_type_{false};
    bool        append_{false};

    // ---------- PerRank ----------
    std::ofstream perrank_xyz_;
    std::ofstream perrank_dump_;
    std::string   perrank_xyz_name_;
    std::string   perrank_dump_name_;

    // ---------- GatherRoot (rank 0 only) ----------
    std::ofstream gather_xyz_;
    std::ofstream gather_dump_;
    std::string   gather_xyz_name_;
    std::string   gather_dump_name_;

    // ---------- MPI-IO (separate files/offsets) ----------
    MPI_File   mpiio_xyz_{MPI_FILE_NULL};
    MPI_File   mpiio_dump_{MPI_FILE_NULL};
    MPI_Offset mpiio_xyz_off_{0};
    MPI_Offset mpiio_dump_off_{0};
    bool       mpiio_xyz_open_{false};
    bool       mpiio_dump_open_{false};

    // Helpers
    static std::string rank_filename(const std::string& base, int rank, const char* ext);

    void ensure_perrank_xyz_open_();
    void ensure_perrank_dump_open_();

    void ensure_gather_xyz_open_();
    void ensure_gather_dump_open_();

    void ensure_mpiio_xyz_open_();
    void ensure_mpiio_dump_open_();

    void gather_counts_(int nloc, std::vector<int>& counts, std::vector<int>& displs) const;

    // Implementations per mode
    void write_xyz_perrank_(const std::vector<Particle>& local, const SimulationBox& box);
    void write_xyz_gather_ (const std::vector<Particle>& local, const SimulationBox& box);
    void write_xyz_mpiio_  (const std::vector<Particle>& local, const SimulationBox& box);

    void write_dump_perrank_(const std::vector<Particle>& local, const SimulationBox& box, int ts);
    void write_dump_gather_ (const std::vector<Particle>& local, const SimulationBox& box, int ts);
    void write_dump_mpiio_  (const std::vector<Particle>& local, const SimulationBox& box, int ts);

    // MPI-IO framed writes (text): header on rank 0, bodies on all ranks
    void mpiio_write_frame_(MPI_File fh, MPI_Offset& ofs,
                            const std::string& header,
                            const std::string& local_body);
};

#endif // LOGGING_TRAJ_MPI_H

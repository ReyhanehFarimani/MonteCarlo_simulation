#include "logging_traj_mpi.h"
#include <numeric>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <cstring>

static std::string with_ext(const std::string& base, const char* ext) {
    const size_t n = base.size(), e = std::strlen(ext);
    return (n >= e && base.substr(n - e) == ext) ? base : (base + ext);
}

LoggingTrajMPI::LoggingTrajMPI(std::string basename, Mode mode, MPI_Comm comm, bool append, bool rank_as_type)
    : basename_(std::move(basename)), mode_(mode), comm_(comm),
      rank_as_type_(rank_as_type), append_(append)
{
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);
}

LoggingTrajMPI::~LoggingTrajMPI() { close(); }

void LoggingTrajMPI::close() {
    // PerRank
    if (perrank_xyz_.is_open())  perrank_xyz_.close();
    if (perrank_dump_.is_open()) perrank_dump_.close();
    // Gather
    if (rank_ == 0) {
        if (gather_xyz_.is_open())  gather_xyz_.close();
        if (gather_dump_.is_open()) gather_dump_.close();
    }
    // MPI-IO
    if (mpiio_xyz_open_)  { MPI_File_close(&mpiio_xyz_);  mpiio_xyz_open_  = false; mpiio_xyz_off_  = 0; }
    if (mpiio_dump_open_) { MPI_File_close(&mpiio_dump_); mpiio_dump_open_ = false; mpiio_dump_off_ = 0; }
}

std::string LoggingTrajMPI::rank_filename(const std::string& base, int rank, const char* ext) {
    std::ostringstream oss;
    oss << base << ".rank" << std::setfill('0') << std::setw(4) << rank << ext;
    return oss.str();
}

/* ---------- open helpers ---------- */
void LoggingTrajMPI::ensure_perrank_xyz_open_() {
    if (perrank_xyz_.is_open()) return;
    perrank_xyz_name_ = rank_filename(basename_, rank_, ".xyz");
    perrank_xyz_.open(perrank_xyz_name_, append_ ? (std::ios::out | std::ios::app)
                                                 : (std::ios::out | std::ios::trunc));
    if (!perrank_xyz_) throw std::runtime_error("Cannot open per-rank XYZ: " + perrank_xyz_name_);
}
void LoggingTrajMPI::ensure_perrank_dump_open_() {
    if (perrank_dump_.is_open()) return;
    perrank_dump_name_ = rank_filename(basename_, rank_, ".dump");
    perrank_dump_.open(perrank_dump_name_, append_ ? (std::ios::out | std::ios::app)
                                                   : (std::ios::out | std::ios::trunc));
    if (!perrank_dump_) throw std::runtime_error("Cannot open per-rank DUMP: " + perrank_dump_name_);
}

void LoggingTrajMPI::ensure_gather_xyz_open_() {
    if (rank_ != 0) return;
    if (gather_xyz_.is_open()) return;
    gather_xyz_name_ = with_ext(basename_, ".xyz");
    gather_xyz_.open(gather_xyz_name_, append_ ? (std::ios::out | std::ios::app)
                                               : (std::ios::out | std::ios::trunc));
    if (!gather_xyz_) throw std::runtime_error("Cannot open gather XYZ: " + gather_xyz_name_);
}
void LoggingTrajMPI::ensure_gather_dump_open_() {
    if (rank_ != 0) return;
    if (gather_dump_.is_open()) return;
    gather_dump_name_ = with_ext(basename_, ".dump");
    gather_dump_.open(gather_dump_name_, append_ ? (std::ios::out | std::ios::app)
                                                 : (std::ios::out | std::ios::trunc));
    if (!gather_dump_) throw std::runtime_error("Cannot open gather DUMP: " + gather_dump_name_);
}

void LoggingTrajMPI::ensure_mpiio_xyz_open_() {
    if (mpiio_xyz_open_) return;
    const std::string name = with_ext(basename_, ".xyz");
    int amode = MPI_MODE_WRONLY | MPI_MODE_CREATE; // manual ofs, no APPEND
    if (MPI_File_open(comm_, name.c_str(), amode, MPI_INFO_NULL, &mpiio_xyz_) != MPI_SUCCESS)
        throw std::runtime_error("Cannot open MPI-IO XYZ: " + name);
    if (append_) { MPI_Offset sz = 0; MPI_File_get_size(mpiio_xyz_, &sz); mpiio_xyz_off_ = sz; }
    mpiio_xyz_open_ = true;
}
void LoggingTrajMPI::ensure_mpiio_dump_open_() {
    if (mpiio_dump_open_) return;
    const std::string name = with_ext(basename_, ".dump");
    int amode = MPI_MODE_WRONLY | MPI_MODE_CREATE; // manual ofs, no APPEND
    if (MPI_File_open(comm_, name.c_str(), amode, MPI_INFO_NULL, &mpiio_dump_) != MPI_SUCCESS)
        throw std::runtime_error("Cannot open MPI-IO DUMP: " + name);
    if (append_) { MPI_Offset sz = 0; MPI_File_get_size(mpiio_dump_, &sz); mpiio_dump_off_ = sz; }
    mpiio_dump_open_ = true;
}

/* ---------- shared helpers ---------- */
void LoggingTrajMPI::gather_counts_(int nloc, std::vector<int>& counts, std::vector<int>& displs) const {
    counts.assign(size_, 0);
    MPI_Gather(&nloc, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm_);
    if (rank_ == 0) {
        displs.resize(size_);
        displs[0] = 0;
        for (int i = 1; i < size_; ++i) displs[i] = displs[i-1] + counts[i-1];
    }
}

/* ===================== XYZ ===================== */
void LoggingTrajMPI::log_xyz(const std::vector<Particle>& local, const SimulationBox& box) {
    switch (mode_) {
        case Mode::PerRank:    write_xyz_perrank_(local, box); break;
        case Mode::GatherRoot: write_xyz_gather_(local, box);  break;
        case Mode::MPIIO:      write_xyz_mpiio_(local, box);   break;
    }
}

void LoggingTrajMPI::write_xyz_perrank_(const std::vector<Particle>& local, const SimulationBox& /*box*/) {
    ensure_perrank_xyz_open_();
    perrank_xyz_ << local.size() << "\n";
    perrank_xyz_ << "rank=" << rank_ << "\n";
    for (const auto& p : local) perrank_xyz_ << "1 " << p.x << " " << p.y << "\n";
    perrank_xyz_.flush();
}

void LoggingTrajMPI::write_xyz_gather_(const std::vector<Particle>& local, const SimulationBox& box) {
    if (rank_ == 0) ensure_gather_xyz_open_();

    const int nloc = (int)local.size();
    std::vector<int> counts, displs;
    gather_counts_(nloc, counts, displs);

    std::vector<double> send(2ull * nloc);
    for (int i = 0; i < nloc; ++i) { send[2*i] = local[i].x; send[2*i+1] = local[i].y; }

    std::vector<int> counts2, displs2;
    std::vector<double> recv;
    int Ntot = 0;
    if (rank_ == 0) {
        counts2.resize(size_);
        displs2.resize(size_);
        for (int i = 0; i < size_; ++i) counts2[i] = 2 * counts[i];
        displs2[0] = 0;
        for (int i = 1; i < size_; ++i) displs2[i] = displs2[i-1] + counts2[i-1];
        Ntot = std::accumulate(counts.begin(), counts.end(), 0);
        recv.resize(2ull * Ntot);
    }

    MPI_Gatherv(send.data(), 2*nloc, MPI_DOUBLE,
                rank_==0 ? recv.data() : nullptr,
                rank_==0 ? counts2.data() : nullptr,
                rank_==0 ? displs2.data() : nullptr,
                MPI_DOUBLE, 0, comm_);

    if (rank_ == 0) {
        gather_xyz_ << Ntot << "\n";
        gather_xyz_ << std::setprecision(10) << "Lx=" << box.getLx() << " Ly=" << box.getLy() << "\n";
        for (int i = 0; i < Ntot; ++i)
            gather_xyz_ << "1 " << recv[2*i] << " " << recv[2*i+1] << "\n";
        gather_xyz_.flush();
    }
}

void LoggingTrajMPI::write_xyz_mpiio_(const std::vector<Particle>& local, const SimulationBox& box) {
    ensure_mpiio_xyz_open_();

    int nloc = (int)local.size(), Ntot = 0;
    MPI_Allreduce(&nloc, &Ntot, 1, MPI_INT, MPI_SUM, comm_);

    std::string header;
    if (rank_ == 0) {
        std::ostringstream h;
        h << Ntot << "\n";
        h << std::setprecision(10) << "Lx=" << box.getLx() << " Ly=" << box.getLy() << "\n";
        header = h.str();
    }

    std::ostringstream oss;
    for (const auto& p : local) oss << "1 " << p.x << " " << p.y << "\n";

    mpiio_write_frame_(mpiio_xyz_, mpiio_xyz_off_, header, oss.str());
}

/* ===================== DUMP ===================== */
void LoggingTrajMPI::log_dump(const std::vector<Particle>& local, const SimulationBox& box, int timestep) {
    switch (mode_) {
        case Mode::PerRank:    write_dump_perrank_(local, box, timestep); break;
        case Mode::GatherRoot: write_dump_gather_(local, box, timestep);  break;
        case Mode::MPIIO:      write_dump_mpiio_ (local, box, timestep);  break;
    }
}

void LoggingTrajMPI::write_dump_perrank_(const std::vector<Particle>& local, const SimulationBox& box, int ts) {
    ensure_perrank_dump_open_();
    perrank_dump_ << "ITEM: TIMESTEP\n" << ts << "\n";
    perrank_dump_ << "ITEM: NUMBER OF ATOMS\n" << local.size() << "\n";
    perrank_dump_ << "ITEM: BOX BOUNDS pp pp ff\n";
    perrank_dump_ << 0.0 << " " << box.getLx() << "\n";
    perrank_dump_ << 0.0 << " " << box.getLy() << "\n";
    perrank_dump_ << 0.0 << " " << 0.0 << "\n";
    perrank_dump_ << "ITEM: ATOMS " << (rank_as_type_ ? "id type x y z\n" : "id x y z\n");
    int id = 1;
    for (const auto& p : local) {
        if (rank_as_type_) perrank_dump_ << id++ << " " << (rank_ + 1) << " " << p.x << " " << p.y << " " << 0.0 << "\n";
        else               perrank_dump_ << id++ << " " << p.x << " " << p.y << " " << 0.0 << "\n";
    }
    perrank_dump_.flush();
}

void LoggingTrajMPI::write_dump_gather_(const std::vector<Particle>& local, const SimulationBox& box, int ts) {
    if (rank_ == 0) ensure_gather_dump_open_();

    const int nloc = (int)local.size();
    std::vector<int> counts, displs;
    gather_counts_(nloc, counts, displs);

    std::vector<double> send(3ull * nloc);
    for (int i = 0; i < nloc; ++i) {
        send[3*i+0] = local[i].x;
        send[3*i+1] = local[i].y;
        send[3*i+2] = 0.0;
    }

    std::vector<int> counts3, displs3;
    std::vector<double> recv;
    int Ntot = 0;
    if (rank_ == 0) {
        counts3.resize(size_);
        displs3.resize(size_);
        for (int i = 0; i < size_; ++i) counts3[i] = 3 * counts[i];
        displs3[0] = 0;
        for (int i = 1; i < size_; ++i) displs3[i] = displs3[i-1] + counts3[i-1];
        Ntot = std::accumulate(counts.begin(), counts.end(), 0);
        recv.resize(3ull * Ntot);
    }

    MPI_Gatherv(send.data(), 3*nloc, MPI_DOUBLE,
                rank_==0 ? recv.data() : nullptr,
                rank_==0 ? counts3.data() : nullptr,
                rank_==0 ? displs3.data() : nullptr,
                MPI_DOUBLE, 0, comm_);

    if (rank_ == 0) {
        gather_dump_ << "ITEM: TIMESTEP\n" << ts << "\n";
        gather_dump_ << "ITEM: NUMBER OF ATOMS\n" << Ntot << "\n";
        gather_dump_ << "ITEM: BOX BOUNDS pp pp ff\n";
        gather_dump_ << 0.0 << " " << box.getLx() << "\n";
        gather_dump_ << 0.0 << " " << box.getLy() << "\n";
        gather_dump_ << 0.0 << " " << 0.0 << "\n";
        gather_dump_ << "ITEM: ATOMS " << (rank_as_type_ ? "id type x y z\n" : "id x y z\n");

        int id = 1;
        int r = 0;
        int upto = (size_ > 0 ? counts[0] : 0);
        for (int i = 0; i < Ntot; ++i) {
            while (rank_as_type_ && i >= upto && r+1 < size_) { ++r; upto += counts[r]; }
            if (rank_as_type_) {
                const int type = r + 1;
                gather_dump_ << id++ << " "
                             << type << " "
                             << recv[3*i+0] << " " << recv[3*i+1] << " " << recv[3*i+2] << "\n";
            } else {
                gather_dump_ << id++ << " "
                             << recv[3*i+0] << " " << recv[3*i+1] << " " << recv[3*i+2] << "\n";
            }
        }
        gather_dump_.flush();
    }
}

void LoggingTrajMPI::write_dump_mpiio_(const std::vector<Particle>& local, const SimulationBox& box, int ts) {
    ensure_mpiio_dump_open_();

    int nloc = (int)local.size(), Ntot = 0;
    MPI_Allreduce(&nloc, &Ntot, 1, MPI_INT, MPI_SUM, comm_);

    std::string header;
    if (rank_ == 0) {
        std::ostringstream h;
        h << "ITEM: TIMESTEP\n" << ts << "\n";
        h << "ITEM: NUMBER OF ATOMS\n" << Ntot << "\n";
        h << "ITEM: BOX BOUNDS pp pp ff\n";
        h << 0.0 << " " << box.getLx() << "\n";
        h << 0.0 << " " << box.getLy() << "\n";
        h << 0.0 << " " << 0.0 << "\n";
        h << "ITEM: ATOMS " << (rank_as_type_ ? "id type x y z\n" : "id x y z\n");
        header = h.str();
    }

    // global 1-based IDs to avoid duplicates across ranks
    int id_base = 0;
    MPI_Exscan(&nloc, &id_base, 1, MPI_INT, MPI_SUM, comm_);
    if (rank_ == 0) id_base = 0;

    std::ostringstream oss;
    oss << std::setprecision(17);
    const int type = rank_as_type_ ? (rank_ + 1) : 0;
    for (int i=0;i<nloc;++i) {
        const int gid = id_base + 1 + i;
        if (rank_as_type_) oss << gid << " " << type << " " << local[i].x << " " << local[i].y << " 0\n";
        else               oss << gid << " "          << local[i].x << " " << local[i].y << " 0\n";
    }

    mpiio_write_frame_(mpiio_dump_, mpiio_dump_off_, header, oss.str());
}

/* ---------- MPI-IO frame writer ---------- */
void LoggingTrajMPI::mpiio_write_frame_(MPI_File fh, MPI_Offset& ofs,
                                        const std::string& header,
                                        const std::string& local_body)
{
    // header size (rank 0 only)
    int hdr_sz = (rank_ == 0) ? static_cast<int>(header.size()) : 0;
    MPI_Bcast(&hdr_sz, 1, MPI_INT, 0, comm_);

    // body sizes and prefix offsets
    int my_sz = static_cast<int>(local_body.size());
    int my_off = 0;
    MPI_Exscan(&my_sz, &my_off, 1, MPI_INT, MPI_SUM, comm_);
    if (rank_ == 0) my_off = 0;

    int sum_sz = 0;
    MPI_Allreduce(&my_sz, &sum_sz, 1, MPI_INT, MPI_SUM, comm_);

    // header at ofs (collective; only rank 0 provides a buffer)
    if (hdr_sz > 0) {
        MPI_Status st;
        const char* hbuf = (rank_ == 0) ? header.data() : nullptr;
        const int    hsz = (rank_ == 0) ? hdr_sz       : 0;
        MPI_File_write_at_all(fh, ofs, hbuf, hsz, MPI_CHAR, &st);
    }

    // bodies at ofs + hdr_sz + prefix
    if (my_sz > 0) {
        MPI_Status st;
        MPI_File_write_at_all(fh, ofs + hdr_sz + my_off, local_body.data(), my_sz, MPI_CHAR, &st);
    }

    ofs += static_cast<MPI_Offset>(hdr_sz + sum_sz);
}

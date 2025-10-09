#ifndef PARTICLE_EXCHANGE_H
#define PARTICLE_EXCHANGE_H

#include <mpi.h>
#include <array>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <limits>
#include <cmath>
#include <type_traits>

#include "particle.h"
#include "simulation_box.h"
#include "cell_list_parallel.h"

// POD payload for MPI
struct PackedParticle {
    double x, y;
    std::int64_t id;
};

class ParticleExchange {
public:
    struct Params {
        bool enable_incremental{false}; // per-move ghost updates (optional)
    };

    // Assumptions (enforced):
    // - Halo width == one interior cell in each direction:
    //     halo_wx = cl.dxCell(), halo_wy = cl.dyCell()
    // - Ghosts occupy the ring cells only (i==0 or nx_in+1, j==0 or ny_in+1).
    ParticleExchange(MPI_Comm comm,
                     const SimulationBox& box,
                     const SimulationBox::Decomposition& decomp,
                     int rank,
                     const CellListParallel& cl,
                     const Params& params = {false});

    // Rebuild ghost ring from current owned border particles.
    void refreshGhosts(const std::vector<Particle>& owned, CellListParallel& cl);

    // Incremental (optional): queue a moved owned particle; later flush.
    void queueGhostUpdateCandidate(const Particle& moved);
    void flushGhostUpdates(CellListParallel& cl);

    // Migrate owners that left [x0,x1)×[y0,y1); returns # sent out.
    // NOTE: Does NOT update ghosts. Call refreshGhosts() afterwards.
    int migrate(std::vector<Particle>& owned, CellListParallel& cl);

    [[nodiscard]] const std::vector<PackedParticle>& getGhosts() const noexcept { return ghosts_; }

private:
    // ---- MPI / geometry (immutable after ctor) ----
    MPI_Comm comm_;
    int rank_{0}, size_{1};
    int Px_{1}, Py_{1}, rx_{0}, ry_{0};
    double Lx_{0.0}, Ly_{0.0};
    double x0_{0.0}, x1_{0.0}, y0_{0.0}, y1_{0.0}; // local owner bounds
    double halo_wx_{0.0}, halo_wy_{0.0};          // = cl.dxCell(), cl.dyCell()
    bool incremental_{false};

    // neighbor ranks: E,W,N,S, NE,NW,SE,SW (fixed order)
    std::array<int,8> nbr_{};

    // ---- Ghost store (authoritative mirror on this rank) ----
    std::vector<PackedParticle> ghosts_;                // dense
    std::unordered_map<std::int64_t,int> ghost_idx_;    // id -> index in ghosts_

    // ---- Exchange scratch ----
    std::array<std::vector<PackedParticle>,8> send_buf_{}, recv_buf_{};
    std::array<int,8> sc_{}, rc_{};

    // Migration scratch
    std::vector<std::vector<PackedParticle>> mig_send_; // size=size_
    std::vector<PackedParticle> mig_recv_flat_;
    std::vector<int> a2_sc_, a2_rc_, a2_sd_, a2_rd_;

private:
    // Wrap-safe distance helpers (noexcept)
    inline double d_to_left_(double x)  const noexcept {
        double d = std::fmod(x - x0_ + Lx_, Lx_); if (d < 0) d += Lx_; return d;
    }
    inline double d_to_right_(double x) const noexcept {
        double d = std::fmod(x1_ - x + Lx_, Lx_); if (d < 0) d += Lx_; return d;
    }
    inline double d_to_bottom_(double y) const noexcept {
        double d = std::fmod(y - y0_ + Ly_, Ly_); if (d < 0) d += Ly_; return d;
    }
    inline double d_to_top_(double y) const noexcept {
        double d = std::fmod(y1_ - y + Ly_, Ly_); if (d < 0) d += Ly_; return d;
    }

    // setup
    void readGeom_(const SimulationBox& box,
                   const SimulationBox::Decomposition& decomp,
                   int rank);
    void buildNeighbors_() noexcept;

    // halo tests (strict/loose pairing to avoid double-sends when strips meet)
    bool in_left_border_ (double x) const noexcept;  // <
    bool in_right_border_(double x) const noexcept;  // <=
    bool in_bottom_border_(double y) const noexcept; // <
    bool in_top_border_   (double y) const noexcept; // <=

    // batch refresh helpers
    void pack_border_owned_(const std::vector<Particle>& owned);
    void exchange_counts_payloads_();
    void replace_ghost_store_from_recv_();
    void push_store_into_celllist_(CellListParallel& cl);

    // incremental updates
    void merge_incremental_into_store_();

    // migration helpers
    inline bool inside_owner_(double x, double y) const noexcept {
        return (x >= x0_ && x < x1_ && y >= y0_ && y < y1_);
    }
    int  dest_rank_for_(double x, double y) const noexcept;
    void alltoallv_migration_();

    // utils
    static inline double wrap_(double a, double L) noexcept {
        a = std::fmod(a, L); if (a < 0) a += L; return a;
    }

    // MPI datatype
    static MPI_Datatype mpi_packed_particle_type(); // cached
};

#endif // PARTICLE_EXCHANGE_H

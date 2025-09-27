/*
================================================================================
 ParticleExchange — Ghost/Halo Exchange and Owner Migration for 2D MPI Domains
================================================================================

ASCII OVERVIEW (Periodic Px × Py process grid)
----------------------------------------------

        N
   NW   ↑   NE
     \  |  /
  W ←-- R --→ E
     /  |  \
   SW   ↓   SE
        S

- Each MPI rank "R" owns a rectangular subdomain [x0,x1) × [y0,y1).
- Global box is periodic with size Lx × Ly, decomposed into Px × Py tiles.
- We maintain:
  - OWNERS: particles whose home lies in R's subdomain.
  - GHOSTS: read-only mirrors of *border* owners from 8 adjacent neighbors.
- Ghosts occupy a *one-cell-thick ring* in the cell list, with halo widths:
    halo_wx = cl.dxCell(), halo_wy = cl.dyCell()

HALO REFRESH (batch)
--------------------
1) PACK: scan OWNERS near each face/corner (<= 1 interior cell) into 8 send buffers.
2) EXCHANGE: two-phase with 8 neighbors (counts, then payloads).
3) REPLACE: rebuild local ghost store (dedup by particle id).
4) PUSH: clear and refill ghost bins in CellListParallel, then rebuild ghost bins.

INCREMENTAL UPDATES (optional)
------------------------------
- Instead of full rebuild every step, you can queue moved owners and flush:
  queueGhostUpdateCandidate() + flushGhostUpdates()
- flush performs an exchange and *in-place* ghost updates (no duplication),
  then re-bins ghost cells.

MIGRATION (owners)
------------------
- Owners that leave [x0,x1) × [y0,y1) are sent to their destination rank
  (computed from global wrapped coordinates).
- Uses Alltoallv over PackedParticle arrays, then:
  - Remove out-migrants locally
  - Append in-migrants (wrapped to [0,L))
  - Rebuild interior bins in CellListParallel

CORRECTNESS PITFALLS ADDRESSED
------------------------------
A) MPI TAGS: Both peers must use matching tags. We use uniform tags for counts
   and payloads across all 8 neighbors to avoid mismatches.

B) Alltoallv ELEMENT SIZE: We define an explicit MPI_Datatype for PackedParticle
   so counts/displacements are in *elements*, not bytes.

C) POD STABILITY: PackedParticle is trivially copyable (static_assert). If you
   ever change layout, update the MPI datatype builder accordingly.

================================================================================
*/

#include "particle_exchange.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>

// =========================== Utility: MPI datatype ============================
// Build and cache an MPI_Datatype for PackedParticle so Alltoallv uses element
// counts/displacements (cleaner and safer than byte arithmetic).
namespace {
    inline MPI_Datatype build_mpi_packed_particle_type() {
        static_assert(std::is_trivially_copyable<PackedParticle>::value,
                      "PackedParticle must be trivially copyable.");

        MPI_Datatype dtype;
        // Describe the struct layout to MPI using relative displacements.
        int blocklen[3] = {1, 1, 1};
        MPI_Aint disp[3];
        MPI_Datatype types[3] = { MPI_DOUBLE, MPI_DOUBLE, MPI_LONG_LONG };

        PackedParticle probe{};
        MPI_Aint base;
        MPI_Get_address(&probe, &base);
        MPI_Get_address(&probe.x, &disp[0]);
        MPI_Get_address(&probe.y, &disp[1]);
        MPI_Get_address(&probe.id, &disp[2]);
        for (int k=0; k<3; ++k) disp[k] -= base;

        MPI_Type_create_struct(3, blocklen, disp, types, &dtype);
        MPI_Type_commit(&dtype);
        return dtype;
    }

    inline MPI_Datatype mpi_packed_particle_type() {
        static MPI_Datatype t = build_mpi_packed_particle_type();
        return t;
    }
}

// ============================ ctor / geometry =================================

ParticleExchange::ParticleExchange(MPI_Comm comm,
                                   const SimulationBox& box,
                                   const SimulationBox::Decomposition& decomp,
                                   int rank,
                                   const CellListParallel& cl,
                                   const Params& params)
: comm_(comm), incremental_(params.enable_incremental)
{
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);

    readGeom_(box, decomp, rank);
    buildNeighbors_();

    // Lock halo width to ONE interior cell.
    halo_wx_ = cl.dxCell();
    halo_wy_ = cl.dyCell();
    if (halo_wx_ <= 0.0 || halo_wy_ <= 0.0) {
        throw std::runtime_error("ParticleExchange: invalid halo widths (dx/dy).");
    }

    // Some light initial reservations to curb realloc churn.
    ghosts_.reserve(256);
    for (auto& v : send_buf_) v.reserve(128);
    for (auto& v : recv_buf_) v.reserve(128);
}

void ParticleExchange::readGeom_(const SimulationBox& box,
                                 const SimulationBox::Decomposition& decomp,
                                 int rank)
{
    Lx_ = box.getLx();
    Ly_ = box.getLy();
    Px_ = decomp.Px; Py_ = decomp.Py;

    auto rc = decomp.coordsOf(rank);
    rx_ = rc.first; ry_ = rc.second;

    double x0,x1,y0,y1;
    decomp.localBounds(rank, Lx_, Ly_, x0, x1, y0, y1);
    x0_ = x0; x1_ = x1; y0_ = y0; y1_ = y1;
}

void ParticleExchange::buildNeighbors_()
{
    auto id = [&](int ix, int iy){
        ix = (ix%Px_ + Px_)%Px_; // periodic wrap
        iy = (iy%Py_ + Py_)%Py_;
        return iy*Px_ + ix;
    };
    // Fixed order: E,W,N,S, NE,NW,SE,SW
    nbr_[0] = id(rx_+1, ry_);   // E
    nbr_[1] = id(rx_-1, ry_);   // W
    nbr_[2] = id(rx_,   ry_+1); // N
    nbr_[3] = id(rx_,   ry_-1); // S
    nbr_[4] = id(rx_+1, ry_+1); // NE
    nbr_[5] = id(rx_-1, ry_+1); // NW
    nbr_[6] = id(rx_+1, ry_-1); // SE
    nbr_[7] = id(rx_-1, ry_-1); // SW
}

// =================== border tests (one-cell ring around owners) =================
//
// Tests classify an *owned* particle as belonging to the one-cell-thick border
// along each face. We include periodic wrap checks to be robust near global
// boundaries; wrap_ normalizes when we finally store ghosts.
//
bool ParticleExchange::in_left_border_(double x) const {
    return (x - x0_) <= halo_wx_ || (x > x1_ && (x0_ + Lx_ - x) <= halo_wx_);
}
bool ParticleExchange::in_right_border_(double x) const {
    return (x1_ - x) <= halo_wx_ || (x < x0_ && (x + Lx_ - x1_) <= halo_wx_);
}
bool ParticleExchange::in_bottom_border_(double y) const {
    return (y - y0_) <= halo_wy_ || (y > y1_ && (y0_ + Ly_ - y) <= halo_wy_);
}
bool ParticleExchange::in_top_border_(double y) const {
    return (y1_ - y) <= halo_wy_ || (y < y0_ && (y + Ly_ - y1_) <= halo_wy_);
}

// ============================ batch HALO refresh ===============================
//
// 1) pack_border_owned_   — classify border owners into 8 directional buffers
// 2) exchange_counts_payloads_ — 8-neighbor exchange (counts, then payloads)
// 3) replace_ghost_store_from_recv_ — rebuild ghost store with dedup by id
// 4) push_store_into_celllist_      — feed to CellListParallel and rebuild ghost bins
//
void ParticleExchange::pack_border_owned_(const std::vector<Particle>& owned)
{
    for (int i=0;i<8;++i){ send_buf_[i].clear(); sc_[i]=0; }

    for (const auto& p : owned) {
        const bool L = in_left_border_ (p.x);
        const bool R = in_right_border_(p.x);
        const bool B = in_bottom_border_(p.y);
        const bool T = in_top_border_   (p.y);

        if (R) send_buf_[0].push_back({p.x,p.y,p.id}); // E
        if (L) send_buf_[1].push_back({p.x,p.y,p.id}); // W
        if (T) send_buf_[2].push_back({p.x,p.y,p.id}); // N
        if (B) send_buf_[3].push_back({p.x,p.y,p.id}); // S
        if (R && T) send_buf_[4].push_back({p.x,p.y,p.id}); // NE
        if (L && T) send_buf_[5].push_back({p.x,p.y,p.id}); // NW
        if (R && B) send_buf_[6].push_back({p.x,p.y,p.id}); // SE
        if (L && B) send_buf_[7].push_back({p.x,p.y,p.id}); // SW
    }
    for (int i=0;i<8;++i) sc_[i] = (int)send_buf_[i].size();
}

void ParticleExchange::exchange_counts_payloads_()
{
    // Use uniform tags for both directions; prevents mismatch (E vs W, etc.)
    constexpr int TAG_CNT = 100; // counts phase
    constexpr int TAG_PAY = 200; // payload phase

    // --- Phase 1: exchange counts with all 8 neighbors (nonblocking) ---
    std::array<MPI_Request,16> rq{};
    for (int i=0;i<8;++i) {
        MPI_Irecv(&rc_[i], 1, MPI_INT, nbr_[i], TAG_CNT, comm_, &rq[i]);
        MPI_Isend(&sc_[i], 1, MPI_INT, nbr_[i], TAG_CNT, comm_, &rq[8+i]);
    }
    MPI_Waitall(16, rq.data(), MPI_STATUSES_IGNORE);

    // Prepare receive buffers
    for (int i=0;i<8;++i) recv_buf_[i].resize(rc_[i]);

    // --- Phase 2: exchange payloads using the PackedParticle MPI datatype ---
    std::array<MPI_Request,16> rq2{};
    MPI_Datatype T = mpi_packed_particle_type();

    for (int i=0;i<8;++i) {
        MPI_Irecv(recv_buf_[i].data(), rc_[i], T, nbr_[i], TAG_PAY, comm_, &rq2[i]);
        MPI_Isend(send_buf_[i].data(), sc_[i], T, nbr_[i], TAG_PAY, comm_, &rq2[8+i]);
    }
    MPI_Waitall(16, rq2.data(), MPI_STATUSES_IGNORE);
}

void ParticleExchange::replace_ghost_store_from_recv_()
{
    ghosts_.clear();
    ghost_idx_.clear();

    // Deduplicate by id; last copy wins (e.g., face + corner overlap)
    for (int i=0;i<8;++i) {
        for (const auto& q : recv_buf_[i]) {
            auto it = ghost_idx_.find(q.id);
            if (it == ghost_idx_.end()) {
                int k = (int)ghosts_.size();
                ghosts_.push_back(q);
                ghost_idx_.emplace(q.id, k);
            } else {
                ghosts_[it->second] = q;
            }
        }
    }
}

void ParticleExchange::push_store_into_celllist_(CellListParallel& cl)
{
    // Normalize positions into [0,L) and push as ghosts into one-cell ring bins
    cl.clearGhosts();
    for (const auto& g : ghosts_) {
        cl.addGhost(Particle{ wrap_(g.x, Lx_), wrap_(g.y, Ly_), static_cast<int>(g.id) });
    }
    cl.rebuildGhostBins();
}

void ParticleExchange::refreshGhosts(const std::vector<Particle>& owned, CellListParallel& cl)
{
    pack_border_owned_(owned);
    exchange_counts_payloads_();
    replace_ghost_store_from_recv_();
    push_store_into_celllist_(cl);
}

// ======================== incremental (per-move) updates =======================
//
// queueGhostUpdateCandidate(): classify a moved owner into 8 directional buffers.
// flushGhostUpdates():
//   - exchange with neighbors based on queued candidates
//   - in-place merge into ghost store (no duplication)
//   - re-bin ghosts in CellListParallel
//
void ParticleExchange::queueGhostUpdateCandidate(const Particle& p)
{
    if (!incremental_) return;

    const bool L = in_left_border_ (p.x);
    const bool R = in_right_border_(p.x);
    const bool B = in_bottom_border_(p.y);
    const bool T = in_top_border_   (p.y);

    if (R) send_buf_[0].push_back({p.x,p.y,p.id}); // E
    if (L) send_buf_[1].push_back({p.x,p.y,p.id}); // W
    if (T) send_buf_[2].push_back({p.x,p.y,p.id}); // N
    if (B) send_buf_[3].push_back({p.x,p.y,p.id}); // S
    if (R && T) send_buf_[4].push_back({p.x,p.y,p.id}); // NE
    if (L && T) send_buf_[5].push_back({p.x,p.y,p.id}); // NW
    if (R && B) send_buf_[6].push_back({p.x,p.y,p.id}); // SE
    if (L && B) send_buf_[7].push_back({p.x,p.y,p.id}); // SW
}

void ParticleExchange::merge_incremental_into_store_()
{
    for (int i=0;i<8;++i) {
        for (const auto& q : recv_buf_[i]) {
            auto it = ghost_idx_.find(q.id);
            if (it == ghost_idx_.end()) {
                int k = (int)ghosts_.size();
                ghosts_.push_back(q);
                ghost_idx_.emplace(q.id, k);
            } else {
                ghosts_[it->second] = q; // modify in place (no duplication)
            }
        }
    }
}

void ParticleExchange::flushGhostUpdates(CellListParallel& cl)
{
    if (!incremental_) return;

    for (int i=0;i<8;++i) sc_[i] = (int)send_buf_[i].size();

    exchange_counts_payloads_();     // fills rc_ + recv_buf_
    merge_incremental_into_store_(); // modify existing ghosts in place

    for (int i=0;i<8;++i) send_buf_[i].clear();

    push_store_into_celllist_(cl);   // re-bin ghosts in ring cells
}

// ================================ migration ===================================
//
// Steps:
// 1) For each owner outside local bounds, compute destination rank via wrapped
//    coordinates → grid cell index.
// 2) Pack per-destination send buffers (PackedParticle).
// 3) Alltoallv over the PackedParticle MPI datatype.
// 4) Erase out-migrants locally; append in-migrants (wrapped to [0,L)).
// 5) Rebuild interior bins for owners.
//
int ParticleExchange::dest_rank_for_(double x, double y) const
{
    x = wrap_(x, Lx_); y = wrap_(y, Ly_);
    const double w = Lx_ / Px_, h = Ly_ / Py_;
    // Map to grid cell indices; clamp at upper edge to Px_-1 / Py_-1
    int ix = std::min((int)(x / w), Px_-1);
    int iy = std::min((int)(y / h), Py_-1);
    return iy*Px_ + ix;
}

void ParticleExchange::alltoallv_migration_()
{
    MPI_Datatype T = mpi_packed_particle_type();

    // Counts (elements)
    a2_sc_.assign(size_, 0);
    for (int r=0; r<size_; ++r) a2_sc_[r] = (int)mig_send_[r].size();

    a2_rc_.assign(size_, 0);
    MPI_Alltoall(a2_sc_.data(), 1, MPI_INT, a2_rc_.data(), 1, MPI_INT, comm_);

    // Displacements (elements)
    a2_sd_.assign(size_, 0);
    a2_rd_.assign(size_, 0);
    int stot = 0, rtot = 0;
    for (int r=0; r<size_; ++r) {
        a2_sd_[r] = stot; stot += a2_sc_[r];
        a2_rd_[r] = rtot; rtot += a2_rc_[r];
    }

    // Flatten sends; prepare receive flat buffer
    std::vector<PackedParticle> send_flat(stot);
    for (int r=0; r<size_; ++r) {
        std::copy(mig_send_[r].begin(), mig_send_[r].end(),
                  send_flat.begin() + a2_sd_[r]);
    }
    mig_recv_flat_.resize(rtot);

    // Exchange in *elements* using the explicit datatype
    MPI_Alltoallv(send_flat.data(), a2_sc_.data(), a2_sd_.data(), T,
                  mig_recv_flat_.data(), a2_rc_.data(), a2_rd_.data(), T,
                  comm_);
}

int ParticleExchange::migrate(std::vector<Particle>& owned, CellListParallel& cl)
{
    mig_send_.assign(size_, {});
    int out = 0;

    // Partition out-migrants
    for (const auto& p : owned) {
        const bool inside = (p.x>=x0_ && p.x<x1_ && p.y>=y0_ && p.y<y1_);
        if (!inside) {
            int dest = dest_rank_for_(p.x, p.y);
            mig_send_[dest].push_back({p.x, p.y, p.id});
            ++out;
        }
    }

    alltoallv_migration_();

    // Remove local out-migrants
    owned.erase(std::remove_if(owned.begin(), owned.end(), [&](const Particle& p){
        return !(p.x>=x0_ && p.x<x1_ && p.y>=y0_ && p.y<y1_);
    }), owned.end());

    // Append in-migrants (wrap to [0,L))
    owned.reserve(owned.size() + mig_recv_flat_.size());
    for (const auto& q : mig_recv_flat_) {
        owned.push_back(Particle{ wrap_(q.x, Lx_), wrap_(q.y, Ly_), static_cast<int>(q.id) });
    }

    // Rebuild owner bins
    cl.buildInterior(owned);

    return out;
}

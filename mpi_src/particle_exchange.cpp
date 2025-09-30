#include "particle_exchange.h"

#include <algorithm>
#include <stdexcept>

// =========================== Utility: MPI datatype ============================
// Build and cache an MPI_Datatype for PackedParticle. Use MPI_Type_match_size
// for portable 64-bit integer selection.
namespace {
    inline MPI_Datatype build_mpi_packed_particle_type_portable() {
        static_assert(std::is_trivially_copyable<PackedParticle>::value,
                      "PackedParticle must be trivially copyable.");

        MPI_Datatype int64_dt;
        MPI_Type_match_size(MPI_TYPECLASS_INTEGER, 8, &int64_dt);

        int blocklen[3] = {1, 1, 1};
        MPI_Aint disp[3], base;

        PackedParticle probe{};
        MPI_Get_address(&probe,    &base);
        MPI_Get_address(&probe.x,  &disp[0]);
        MPI_Get_address(&probe.y,  &disp[1]);
        MPI_Get_address(&probe.id, &disp[2]);
        for (int k=0; k<3; ++k) disp[k] -= base;

        MPI_Datatype dtype;
        MPI_Datatype types[3] = { MPI_DOUBLE, MPI_DOUBLE, int64_dt };
        MPI_Type_create_struct(3, blocklen, disp, types, &dtype);
        MPI_Type_commit(&dtype);
        return dtype;
    }
}
MPI_Datatype ParticleExchange::mpi_packed_particle_type() {
    static MPI_Datatype t = build_mpi_packed_particle_type_portable();
    return t;
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

    // Consistency checks: communicator vs decomposition
    if (size_ != Px_ * Py_) {
        throw std::runtime_error("ParticleExchange: MPI comm size != Px*Py");
    }
    const int expected_lin = ry_ * Px_ + rx_;
    if (expected_lin != rank_) {
        throw std::runtime_error("ParticleExchange: Decomposition rank != MPI rank");
    }

    // Lock halo width to ONE interior cell; enforce invariants.
    halo_wx_ = cl.dxCell();
    halo_wy_ = cl.dyCell();

    if (halo_wx_ < cl.rcut() || halo_wy_ < cl.rcut() || halo_wx_ <= 0.0 || halo_wy_ <= 0.0) {
        throw std::runtime_error("ParticleExchange: halo width must be >= rcut and > 0");
    }
    if ((cl.nxInterior() % 2) || (cl.nyInterior() % 2)) {
        throw std::runtime_error("ParticleExchange: CellList interior counts must be even");
    }

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

void ParticleExchange::buildNeighbors_() noexcept
{
    auto id = [&](int ix, int iy){
        ix = (ix%Px_ + Px_)%Px_;
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
// Use strict/loose pairing to avoid duplicate sends when bands meet.
bool ParticleExchange::in_left_border_ (double x) const noexcept { return d_to_left_ (x)  <  halo_wx_; }
bool ParticleExchange::in_right_border_(double x) const noexcept { return d_to_right_(x) <= halo_wx_; }
bool ParticleExchange::in_bottom_border_(double y) const noexcept{ return d_to_bottom_(y) <  halo_wy_; }
bool ParticleExchange::in_top_border_   (double y) const noexcept{ return d_to_top_   (y) <= halo_wy_; }

// ============================ batch HALO refresh ===============================

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
    constexpr int TAG_CNT = 100; // counts
    constexpr int TAG_PAY = 200; // payloads

    // --- Phase 1: counts ---
    std::array<MPI_Request,16> rq{};
    for (int i=0;i<8;++i) {
        MPI_Irecv(&rc_[i], 1, MPI_INT, nbr_[i], TAG_CNT, comm_, &rq[i]);
        MPI_Isend(&sc_[i], 1, MPI_INT, nbr_[i], TAG_CNT, comm_, &rq[8+i]);
    }
    MPI_Waitall(16, rq.data(), MPI_STATUSES_IGNORE);

    for (int i=0;i<8;++i) recv_buf_[i].resize(rc_[i]);

    // --- Phase 2: payloads ---
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
    cl.clearGhosts();
    for (const auto& g : ghosts_) {
        if (g.id < std::numeric_limits<int>::min() || g.id > std::numeric_limits<int>::max()) {
            throw std::runtime_error("ParticleExchange: ghost id out of 32-bit range for CellListParallel::Particle.id");
        }
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
                ghosts_[it->second] = q; // modify in place
            }
        }
    }
}

void ParticleExchange::flushGhostUpdates(CellListParallel& cl)
{
    if (!incremental_) return;

    for (int i=0;i<8;++i) sc_[i] = (int)send_buf_[i].size();

    exchange_counts_payloads_();     // fills rc_ + recv_buf_
    merge_incremental_into_store_(); // modify store
    for (int i=0;i<8;++i) send_buf_[i].clear();

    push_store_into_celllist_(cl);   // re-bin ghosts
}

// ================================ migration ===================================

int ParticleExchange::dest_rank_for_(double x, double y) const noexcept
{
    x = wrap_(x, Lx_); y = wrap_(y, Ly_);
    const double w = Lx_ / Px_, h = Ly_ / Py_;
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

    // Exchange using the explicit datatype (counts are in elements)
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
        if (!inside_owner_(p.x, p.y)) {
            int dest = dest_rank_for_(p.x, p.y);
            mig_send_[dest].push_back({p.x, p.y, p.id});
            ++out;
        }
    }

    alltoallv_migration_();

    // Remove local out-migrants
    owned.erase(std::remove_if(owned.begin(), owned.end(),
                   [&](const Particle& p){ return !inside_owner_(p.x, p.y); }),
                owned.end());

    // Append in-migrants (wrap to [0,L))
    owned.reserve(owned.size() + mig_recv_flat_.size());
    for (const auto& q : mig_recv_flat_) {
        owned.emplace_back(Particle{ wrap_(q.x, Lx_), wrap_(q.y, Ly_), static_cast<int>(q.id) });
    }

    // Rebuild owner bins
    cl.buildInterior(owned);

    // Caller should call refreshGhosts(owned, cl) next.
    return out;
}

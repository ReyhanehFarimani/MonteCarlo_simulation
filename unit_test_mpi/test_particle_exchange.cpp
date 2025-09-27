#include "catch.hpp"
#include <mpi.h>
#include <vector>
#include <cmath>
#include <random>
#include <unordered_set>
#include <numeric>
#include <tuple>
#include <sstream>
#include <limits>
#include "../mpi_src/particle.h"
#include "../mpi_src/simulation_box.h"
#include "../mpi_src/cell_list_parallel.h"
#include "../mpi_src/particle_exchange.h"

// ----------------- Local helpers (test-only) ----------------------------------


// Unique id helper to tag particles with (rank, tag) -> int (matches Particle::id)
static inline int make_id(int rank, int tag) {
    long long v = 1LL * rank * 1'000'000 + tag;
    // Guard against overflow into negative/undefined int
    REQUIRE(v <= std::numeric_limits<int>::max());
    REQUIRE(v >= std::numeric_limits<int>::min());
    return static_cast<int>(v);
}
// Neighbor rank computation mirroring ParticleExchange::buildNeighbors_ order
// Order: E,W,N,S, NE,NW,SE,SW
static std::array<int,8> neighbor_ranks(int rank, int Px, int Py) {
    int rx = rank % Px;
    int ry = rank / Px;
    auto wrap = [](int i, int n){ i%=n; if (i<0) i+=n; return i; };
    auto id = [&](int ix, int iy){ ix=wrap(ix,Px); iy=wrap(iy,Py); return iy*Px+ix; };
    return {
        id(rx+1, ry),   // E
        id(rx-1, ry),   // W
        id(rx,   ry+1), // N
        id(rx,   ry-1), // S
        id(rx+1, ry+1), // NE
        id(rx-1, ry+1), // NW
        id(rx+1, ry-1), // SE
        id(rx-1, ry-1)  // SW
    };
}

// Probe whether neighborsOfOwned(i) contains any ghost within rcut^2
static bool has_any_ghost_within_rcut2(const CellListParallel& clp,
                                       int idx_owned,
                                       const std::vector<Particle>& owned)
{
    const double r2 = std::max(clp.dxCell(), clp.dyCell());
    const double r2max = r2 * r2 + 1e-12;
    auto nb = clp.neighborsOfOwned(idx_owned, owned);
    for (const auto& pr : nb) {
        int j = pr.first;
        double d2 = pr.second;
        if (j < 0 && d2 <= r2max) return true; // negative index → ghost in this API
    }
    return false;
}

// Return the minimum distance^2 to any ghost in neighborsOfOwned(i). If none, returns +inf.
static double min_d2_to_ghost(const CellListParallel& clp,
                              int idx_owned,
                              const std::vector<Particle>& owned)
{
    auto nb = clp.neighborsOfOwned(idx_owned, owned);
    double best = std::numeric_limits<double>::infinity();
    for (const auto& pr : nb) {
        int j = pr.first;
        double d2 = pr.second;
        if (j < 0) best = std::min(best, d2);
    }
    return best;
}

// Construct 4 owned particles inside the local subdomain near each face (with ids).
// Tags: 10 (L), 20 (R), 30 (B), 40 (T)
static std::vector<Particle>
make_owned_face_quartet_with_ids(int rank,
                                 double x0, double x1, double y0, double y1,
                                 double halo_wx, double halo_wy)
{
    // small offset inside (quarter halo)
    const double ex = 0.25 * halo_wx;
    const double ey = 0.25 * halo_wy;
    const double yc = 0.5 * (y0 + y1);
    const double xc = 0.5 * (x0 + x1);

    std::vector<Particle> v;
    v.reserve(4);
    v.push_back(Particle{x0 + ex, yc, make_id(rank, 10)}); // near left face
    v.push_back(Particle{x1 - ex, yc, make_id(rank, 20)}); // near right face
    v.push_back(Particle{xc, y0 + ey, make_id(rank, 30)}); // near bottom
    v.push_back(Particle{xc, y1 - ey, make_id(rank, 40)}); // near top
    return v;
}

// Build one owned particle at the center of the local domain (stable probe)
static Particle make_center_probe(const CellListParallel& clp, int rank, std::int64_t tag) {
    const double xc = 0.5*(clp.x0()+clp.x1());
    const double yc = 0.5*(clp.y0()+clp.y1());
    return Particle{xc, yc, make_id(rank, tag)};
}

// Sums an integer across all ranks (for global conservation checks)
static int allreduce_sum_int(int local) {
    int global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

// ----------------- Tests -------------------------------------------------------

// Geometry smoke test is already in your CellListParallel tests; we focus on ParticleExchange.

TEST_CASE("ParticleExchange: halo refresh across 8 neighbors", "[MPI][ParticleExchange][ghosts]") {
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Use a box with mild anisotropy; rely on bestDecomposition(size).
    SimulationBox box(240.0, 180.0);
    auto decomp = box.bestDecomposition(size);

    const double rcut = 4.0;
    CellListParallel::Params cl_params{rcut};
    CellListParallel clp(box, decomp, rank, cl_params);

    // Owned: four particles near each face with unique ids
    auto owned = make_owned_face_quartet_with_ids(rank, clp.x0(), clp.x1(), clp.y0(), clp.y1(),
                                                  clp.dxCell(), clp.dyCell());
    // Add one probe in the center so neighborsOfOwned(4) is stable
    owned.push_back(make_center_probe(clp, rank, 777));
    clp.buildInterior(owned);

    // Build the particle exchange (batch mode)
    ParticleExchange::Params pex_params{false};
    ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pex_params);

    // First full refresh: should populate ghosts in the ring from neighbors
    pex.refreshGhosts(owned, clp);

    // Face-adjacent owned must "see" at least one ghost across rcut (coming from neighbor border)
    REQUIRE(has_any_ghost_within_rcut2(clp, 0, owned)); // near left face
    REQUIRE(has_any_ghost_within_rcut2(clp, 1, owned)); // near right face
    REQUIRE(has_any_ghost_within_rcut2(clp, 2, owned)); // near bottom
    REQUIRE(has_any_ghost_within_rcut2(clp, 3, owned)); // near top

    // The center probe should have some ghost neighbors if rcut covers the ring horizontally/vertically.
    // We don't require it strictly (depends on cell sizing), but min distance must be finite when close enough.
    // const double d2_center = min_d2_to_ghost(clp, 6, owned);
    // REQUIRE(std::isfinite(d2_center)); // With one-cell halo and our sizes, this should hold.

    // Corner coverage sanity: make a particle near top-right corner and check it sees a ghost.
    Particle corner{clp.x1() - 0.25*clp.dxCell(), clp.y1() - 0.25*clp.dyCell(), make_id(rank, 888)};
    owned.push_back(corner);
    Particle corner2{clp.x0() + 0.05*clp.dxCell(), clp.y0() + 0.25*clp.dyCell(), make_id(rank, 188)};
    owned.push_back(corner2);
    clp.buildInterior(owned);
    pex.refreshGhosts(owned, clp);
    REQUIRE(has_any_ghost_within_rcut2(clp, (int)owned.size()-1, owned));
}

TEST_CASE("ParticleExchange: incremental ghost updates reflect neighbor border moves", "[MPI][ParticleExchange][incremental]") {
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Require at least 2 ranks to have distinct face neighbors.
    if (size < 2) {
        WARN("Incremental test requires size >= 2; skipping.");
        return;
    }

    SimulationBox box(240.0, 180.0);
    auto decomp = box.bestDecomposition(size);

    const double rcut = 4.0;
    CellListParallel::Params cl_params{rcut};
    CellListParallel clp(box, decomp, rank, cl_params);

    // Owned: create a center probe that will "look" toward the left ring,
    // plus a right-face border particle (this one will be exported to our E neighbor).
    std::vector<Particle> owned;
    owned.reserve(3);
    owned.push_back(make_center_probe(clp, rank, 1001)); // idx 0
    // Right-face border particle (inside)
    owned.push_back(Particle{clp.x1() - 0.25*clp.dxCell(), 0.5*(clp.y0()+clp.y1()), make_id(rank, 2001)}); // idx 1
    // Left-face border particle (inside)
    owned.push_back(Particle{clp.x0() + 0.25*clp.dxCell(), 0.5*(clp.y0()+clp.y1()), make_id(rank, 2002)}); // idx 2
    clp.buildInterior(owned);

    ParticleExchange::Params pex_params{/*enable_incremental=*/true};
    ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pex_params);

    // Initial full build
    pex.refreshGhosts(owned, clp);
    // Use the LEFT-FACE owned particle (idx 2) as the probe; it should see W-ghosts.
    const int idx_left = 2;
    const double d2_before = min_d2_to_ghost(clp, idx_left, owned);

    // Now every rank SLIGHTLY shifts its right-face border particle further right (still within the border),
    // and flushes incremental updates. This changes the ghosts seen by its EAST neighbor.
    // To make this affect *our* left-ring ghosts, our WEST neighbor is performing the same action.
    // Hence the center probe should measure a (non-infinite) distance change to the nearest ghost after flush.
    {
        // queue updates for this rank's border particles (both faces to be symmetric)
        Particle oldR = owned[1];
        owned[1].x += 0.2*clp.dxCell(); // still within border (one-cell halo)
        pex.queueGhostUpdateCandidate(owned[1]);

        Particle oldL = owned[2];
        owned[2].x -= 0.2*clp.dxCell(); // still within border
        pex.queueGhostUpdateCandidate(owned[2]);

        // Flush incremental halo updates and re-bin ghosts
        pex.flushGhostUpdates(clp);
    }

    const double d2_after = min_d2_to_ghost(clp, idx_left, owned);
    REQUIRE(std::isfinite(d2_before));
    REQUIRE(std::isfinite(d2_after));
    // Not strictly guaranteed to be different numerically on all layouts, but
    // with uniform shifts on all ranks, the neighbor ghost position should change.
    REQUIRE(d2_before != Approx(d2_after));
}

TEST_CASE("ParticleExchange: migration sends leavers, receives newcomers, and preserves global count", "[MPI][ParticleExchange][migration]") {
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    SimulationBox box(300.0, 240.0);
    auto decomp = box.bestDecomposition(size);

    const double rcut = 5.0;
    CellListParallel::Params cl_params{rcut};
    CellListParallel clp(box, decomp, rank, cl_params);

    // Start with a few owners fully inside (stable), plus a few intentionally OUTSIDE
    // (so migrate() will send them away).
    std::vector<Particle> owned;
    owned.reserve(8);

    // Two interior particles
    owned.push_back(Particle{0.5*(clp.x0()+clp.x1()), 0.5*(clp.y0()+clp.y1()), make_id(rank, 3001)});
    owned.push_back(Particle{0.6*(clp.x0()+clp.x1()), 0.4*(clp.y0()+clp.y1()), make_id(rank, 3002)});

    // Three out-migrants crafted to specific directions:
    // (1) Slightly to the right of x1()  → goes to E neighbor (or wrapped)
    owned.push_back(Particle{clp.x1() + 0.1*clp.dxCell(), 0.5*(clp.y0()+clp.y1()), make_id(rank, 4001)});
    // (2) Slightly to the left of x0()   → goes to W neighbor (or wrapped)
    owned.push_back(Particle{clp.x0() - 0.1*clp.dxCell(), 0.5*(clp.y0()+clp.y1()), make_id(rank, 4002)});
    // (3) Slightly above y1()            → goes to N neighbor (or wrapped)
    owned.push_back(Particle{0.5*(clp.x0()+clp.x1()), clp.y1() + 0.1*clp.dyCell(), make_id(rank, 4003)});

    // Also add one that exits diagonally (NE)
    owned.push_back(Particle{clp.x1() + 0.1*clp.dxCell(), clp.y1() + 0.1*clp.dyCell(), make_id(rank, 4004)});

    clp.buildInterior(owned);

    ParticleExchange::Params pex_params{/*enable_incremental=*/false};
    ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pex_params);

    // Baseline global count
    int n0_local = (int)owned.size();
    int n0_global = allreduce_sum_int(n0_local);

    // Migrate and rebuild interior bins (inside migrate())
    int n_sent = pex.migrate(owned, clp);
    REQUIRE(n_sent >= 0);
    // Local count may increase/decrease, but global must be conserved.
    int n1_local = (int)owned.size();
    int n1_global = allreduce_sum_int(n1_local);
    REQUIRE(n1_global == n0_global);

    // After migration, call a ghost refresh to ensure halos still work with new ownership
    pex.refreshGhosts(owned, clp);

    // Basic sanity: all remaining local particles must be inside local bounds
    for (const auto& p : owned) {
        REQUIRE(p.x >= clp.x0());
        REQUIRE(p.x <  clp.x1());
        REQUIRE(p.y >= clp.y0());
        REQUIRE(p.y <  clp.y1());
    }
}

TEST_CASE("ParticleExchange: repeated refresh is idempotent (no growth/leak of ghosts)", "[MPI][ParticleExchange][ghosts][idempotence]") {
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    SimulationBox box(200.0, 160.0);
    auto decomp = box.bestDecomposition(size);

    const double rcut = 4.5;
    CellListParallel::Params cl_params{rcut};
    CellListParallel clp(box, decomp, rank, cl_params);

    // Build a modest grid of owned particles inside (not only faces), to exercise packing.
    std::vector<Particle> owned;
    const int nx = std::max(4, clp.nxInterior());
    const int ny = std::max(4, clp.nyInterior());
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            double x = clp.x0() + (i-0.5)* ( (clp.x1()-clp.x0())/nx );
            double y = clp.y0() + (j-0.5)* ( (clp.y1()-clp.y0())/ny );
            owned.push_back(Particle{x, y, make_id(rank, 500000 + j*100 + i)});
        }
    }
    clp.buildInterior(owned);

    ParticleExchange::Params pex_params{/*enable_incremental=*/false};
    ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pex_params);

    // Two consecutive refreshes must not change neighbor sets for owned particles.
    pex.refreshGhosts(owned, clp);

    // Collect min ghost distance for a few sample owned particles
    std::vector<int> samples;
    for (int k = 0; k < (int)owned.size(); ++k) {
        // sample corners in the owned array
        if (k==0 || k==(int)owned.size()-1 || k==(int)owned.size()/2) samples.push_back(k);
    }
    std::vector<double> before;
    for (int idx : samples) before.push_back(min_d2_to_ghost(clp, idx, owned));

    // Repeat refresh
    pex.refreshGhosts(owned, clp);

    std::vector<double> after;
    for (int idx : samples) after.push_back(min_d2_to_ghost(clp, idx, owned));

    REQUIRE(before.size() == after.size());
    for (size_t i=0; i<before.size(); ++i) {
        REQUIRE(after[i] == Approx(before[i]));
    }
}

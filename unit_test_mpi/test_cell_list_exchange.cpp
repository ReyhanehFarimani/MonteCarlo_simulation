
#include "catch.hpp"
#include <mpi.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

#include "../mpi_src/particle.h"
#include "../mpi_src/simulation_box.h"
#include "../mpi_src/cell_list_parallel.h"
#include "../mpi_src/particle_exchange.h"

// --------- helpers (minimal, deterministic) -----------------------------------

static inline int make_id(int rank, int tag) {
    long long v = 1LL * rank * 1'000'000 + tag;
    REQUIRE(v <= std::numeric_limits<int>::max());
    return (int)v;
}

static bool has_ghost_neighbor_within(const CellListParallel& cl,
                                      int owned_idx,
                                      const std::vector<Particle>& owned,
                                      double r2max)
{
    auto nb = cl.neighborsOfOwned(owned_idx, owned);
    for (auto& pr : nb) {
        if (pr.first < 0 && pr.second <= r2max + 1e-12) return true;
    }
    return false;
}

// construct 4 owned particles inside the local subdomain near each face (+ a center probe)
static std::vector<Particle>
owned_face_quartet_plus_center(const CellListParallel& cl, int rank)
{
    const double x0 = cl.x0(), x1 = cl.x1();
    const double y0 = cl.y0(), y1 = cl.y1();
    const double ex = 0.25 * cl.dxCell();
    const double ey = 0.25 * cl.dyCell();
    const double yc = 0.5 * (y0 + y1);
    const double xc = 0.5 * (x0 + x1);

    std::vector<Particle> v;
    v.reserve(5);
    v.push_back(Particle{x0 + ex, yc, make_id(rank, 10)}); // near left
    v.push_back(Particle{x1 - ex, yc, make_id(rank, 20)}); // near right
    v.push_back(Particle{xc, y0 + ey, make_id(rank, 30)}); // near bottom
    v.push_back(Particle{xc, y1 - ey, make_id(rank, 40)}); // near top
    v.push_back(Particle{xc, yc,      make_id(rank, 50)}); // center probe
    return v;
}

// ------------------------------- test -----------------------------------------

TEST_CASE("CellListParallel + ParticleExchange: basic co-operation", "[MPI][integration]") {
    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // geometry
    const double Lx = 200.0, Ly = 160.0;
    const double rcut = 4.0;

    SimulationBox box(Lx, Ly);
    auto decomp = box.bestDecomposition(size);

    // cell-list (enforces dx,dy >= rcut)
    CellListParallel::Params clp_params{rcut};
    CellListParallel cl(box, decomp, rank, clp_params);

    // owned: 4 face-adjacent + 1 center
    auto owned = owned_face_quartet_plus_center(cl, rank);
    cl.buildInterior(owned);

    // particle-exchange (one-cell halo, >= rcut)
    ParticleExchange::Params pex_params{/*enable_incremental=*/true};
    ParticleExchange px(MPI_COMM_WORLD, box, decomp, rank, cl, pex_params);

    // 1) full halo refresh — ghosts should appear across each face
    px.refreshGhosts(owned, cl);
    const double r2 = rcut * rcut;

    REQUIRE(has_ghost_neighbor_within(cl, 0, owned, r2)); // left-owned sees W ghosts
    REQUIRE(has_ghost_neighbor_within(cl, 1, owned, r2)); // right-owned sees E ghosts
    REQUIRE(has_ghost_neighbor_within(cl, 2, owned, r2)); // bottom-owned sees S ghosts
    REQUIRE(has_ghost_neighbor_within(cl, 3, owned, r2)); // top-owned sees N ghosts

    // snapshot min ghost distance for the left-owned particle
    auto min_d2_to_ghost = [&](int idx)->double{
        double best = std::numeric_limits<double>::infinity();
        for (auto& pr : cl.neighborsOfOwned(idx, owned))
            if (pr.first < 0) best = std::min(best, pr.second);
        return best;
    };
    const int idx_left = 0;
    const double d2_before = min_d2_to_ghost(idx_left);
    REQUIRE(std::isfinite(d2_before));

    // 2) incremental updates — shift both left and right face particles *inside* border, flush
    //    This should move ghost mirrors at neighbors and alter our min ghost distance too (by symmetry).
    {
        // left particle move slightly more left (still inside owner tile)
        Particle oldL = owned[0];
        owned[0].x -= 0.2 * cl.dxCell(); box.applyPBC(owned[0]); // inside border
        px.queueGhostUpdateCandidate(owned[0]);

        // right particle move slightly more right (still inside)
        Particle oldR = owned[1];
        owned[1].x += 0.2 * cl.dxCell(); box.applyPBC(owned[1]);
        px.queueGhostUpdateCandidate(owned[1]);

        // update bins for owners and flush halo
        cl.onAcceptedMove(0, oldL, owned[0]);
        cl.onAcceptedMove(1, oldR, owned[1]);
        px.flushGhostUpdates(cl);
    }

    const double d2_after = min_d2_to_ghost(idx_left);
    REQUIRE(std::isfinite(d2_after));
    // The nearest ghost position should have changed due to incremental propagation
    REQUIRE(d2_after != Approx(d2_before));

    // 3) migration — send one intentionally outside, then rebuild halo and check invariants
    //    Put particle slightly to the right of x1() to force migration to E (or wrap if single rank).
    owned.push_back(Particle{cl.x1() + 0.1 * cl.dxCell(), 0.5*(cl.y0()+cl.y1()), make_id(rank, 99)});
    cl.buildInterior(owned);

    int sent = px.migrate(owned, cl); // rebuilds owners' interior bins
    REQUIRE(sent >= 0);

    // After migration, refresh ghosts (API contract)
    px.refreshGhosts(owned, cl);

    // All remaining local particles must be inside local bounds
    for (const auto& p : owned) {
        REQUIRE(p.x >= cl.x0());
        REQUIRE(p.x <  cl.x1());
        REQUIRE(p.y >= cl.y0());
        REQUIRE(p.y <  cl.y1());
    }

    // Face-adjacent particles still see some ghost within rcut after migration+refresh
    REQUIRE(has_ghost_neighbor_within(cl, 0, owned, r2));
    REQUIRE(has_ghost_neighbor_within(cl, 1, owned, r2));
    REQUIRE(has_ghost_neighbor_within(cl, 2, owned, r2));
    REQUIRE(has_ghost_neighbor_within(cl, 3, owned, r2));

    // 4) idempotence — a second refresh should not change min ghost distances for stable configuration
    const double d2_before2 = min_d2_to_ghost(idx_left);
    px.refreshGhosts(owned, cl);
    const double d2_after2  = min_d2_to_ghost(idx_left);
    REQUIRE(d2_after2 == Approx(d2_before2));
}

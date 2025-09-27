#include "catch.hpp"
#include <mpi.h>
#include <vector>
#include <cmath>
#include <random>
#include <unordered_set>

#include "../mpi_src/particle.h"
#include "../mpi_src/simulation_box.h"
#include "../mpi_src/cell_list_parallel.h"

// ----------------- Helpers -----------------------------------------------------

static std::vector<Particle>
make_face_owned(double x0, double x1, double y0, double y1, double rcut) {
    const double eps = 0.25 * rcut;
    const double yc  = 0.5 * (y0 + y1);
    const double xc  = 0.5 * (x0 + x1);

    std::vector<Particle> v;
    v.reserve(4);
    v.push_back(Particle{x0 + eps, yc});      // left face (inside)
    v.push_back(Particle{x1 - eps, yc});      // right face (inside)
    v.push_back(Particle{xc, y0 + eps});      // bottom face (inside)
    v.push_back(Particle{xc, y1 - eps});      // top face (inside)
    return v;
}

static void add_face_ghosts(CellListParallel& clp, double Lx, double Ly,
                            double x0, double x1, double y0, double y1, double rcut)
{
    const double eps = 0.25 * rcut;
    const double yc  = 0.5 * (y0 + y1);
    const double xc  = 0.5 * (x0 + x1);

    auto wrap = [](double x, double L) {
        double k = std::floor(x / L);
        x -= k * L;
        while (x < 0) x += L;
        while (x >= L) x -= L;
        return x;
    };

    clp.addGhost(Particle{wrap(x0 - eps, Lx), yc}); // left ghost (outside)
    clp.addGhost(Particle{wrap(x1 + eps, Lx), yc}); // right ghost
    clp.addGhost(Particle{xc, wrap(y0 - eps, Ly)}); // bottom ghost
    clp.addGhost(Particle{xc, wrap(y1 + eps, Ly)}); // top ghost
    clp.rebuildGhostBins();
}

// Build one owned particle at the center of every interior cell.
static std::vector<Particle>
make_owned_grid_centers(const CellListParallel& clp)
{
    std::vector<Particle> v;
    v.reserve(clp.nxInterior()*clp.nyInterior());
    const double x0 = clp.x0(), y0 = clp.y0();
    const double dx = clp.dxCell(), dy = clp.dyCell();
    for (int j = 1; j <= clp.nyInterior(); ++j) {
        for (int i = 1; i <= clp.nxInterior(); ++i) {
            const double cx = x0 + (i - 0.5)*dx;
            const double cy = y0 + (j - 0.5)*dy;
            v.push_back(Particle{cx, cy});
        }
    }
    return v;
}

// Count neighbors within rcut^2 using the class query
static int count_neighbors_owned_i(const CellListParallel& clp,
                                   int idx,
                                   const std::vector<Particle>& owned)
{
    auto nb = clp.neighborsOfOwned(idx, owned);
    return (int)nb.size();
}

// Compare per-particle neighbor counts before/after rebuild
static void expect_neighbors_unchanged_after_rebuild(CellListParallel& clp,
                                                     std::vector<Particle>& owned)
{
    std::vector<int> before(owned.size(), 0);
    for (int i=0;i<(int)owned.size();++i) before[i] = count_neighbors_owned_i(clp, i, owned);

    // trigger rebuild
    clp.buildInterior(owned);

    std::vector<int> after(owned.size(), 0);
    for (int i=0;i<(int)owned.size();++i) after[i] = count_neighbors_owned_i(clp, i, owned);

    REQUIRE(before.size()==after.size());
    for (size_t i=0;i<before.size();++i) REQUIRE(before[i]==after[i]);
}

// Check that two interior cells are not edge-adjacent (Manhattan distance >= 2) for same parity.
static bool non_adjacent_same_parity(std::pair<int,int> a, std::pair<int,int> b) {
    const int dx = std::abs(a.first  - b.first);
    const int dy = std::abs(a.second - b.second);
    return (dx + dy) >= 2;
}

// ----------------- Tests -------------------------------------------------------

TEST_CASE("CellListParallel geometry and invariants", "[MPI][CellListParallel][geom]") {
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    SimulationBox box(123.0, 87.0);
    auto decomp = box.bestDecomposition(size);

    const double rcut = 3.5;
    CellListParallel::Params params{rcut};
    CellListParallel clp(box, decomp, rank, params);

    REQUIRE(clp.nxInterior() >= 2);
    REQUIRE(clp.nyInterior() >= 2);
    REQUIRE((clp.nxInterior() % 2) == 0);
    REQUIRE((clp.nyInterior() % 2) == 0);
    REQUIRE(clp.nxTotal() == clp.nxInterior() + 2);
    REQUIRE(clp.nyTotal() == clp.nyInterior() + 2);

    const double Llocx = clp.x1() - clp.x0();
    const double Llocy = clp.y1() - clp.y0();
    const int nmax_x = (int)std::floor(Llocx / rcut);
    const int nmax_y = (int)std::floor(Llocy / rcut);

    if (nmax_x >= 2) REQUIRE(clp.dxCell() >= rcut - 1e-12);
    else             REQUIRE(clp.nxInterior() == 2);
    if (nmax_y >= 2) REQUIRE(clp.dyCell() >= rcut - 1e-12);
    else             REQUIRE(clp.nyInterior() == 2);

    // interior mapping
    int ix=-9, iy=-9;
    const double xc = 0.5*(clp.x0()+clp.x1());
    const double yc = 0.5*(clp.y0()+clp.y1());
    REQUIRE(clp.mapToLocalCell(xc, yc, ix, iy));
    REQUIRE(ix >= 1);
    REQUIRE(ix <= clp.nxInterior());
    REQUIRE(iy >= 1);
    REQUIRE(iy <= clp.nyInterior());
}

TEST_CASE("Ghost ring mapping and neighbor hits across faces", "[MPI][CellListParallel][ghosts][neighbors]") {
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    SimulationBox box(200.0, 160.0);
    auto decomp = box.bestDecomposition(size);

    const double rcut = 4.0;
    CellListParallel::Params params{rcut};
    CellListParallel clp(box, decomp, rank, params);

    // Owned: 4 particles near each face (inside)
    auto owned = make_face_owned(clp.x0(), clp.x1(), clp.y0(), clp.y1(), rcut);
    clp.buildInterior(owned);

    // Ghosts: one per face just outside the domain
    clp.clearGhosts();
    add_face_ghosts(clp, box.getLx(), box.getLy(), clp.x0(), clp.x1(), clp.y0(), clp.y1(), rcut);

    // Mapping tests: faces land on ring cells
    {
        int ix, iy;
        REQUIRE(clp.mapToLocalCell(clp.x0() - 0.25*rcut, 0.5*(clp.y0()+clp.y1()), ix, iy));
        REQUIRE(ix == 0);
        REQUIRE(clp.mapToLocalCell(clp.x1() + 0.25*rcut, 0.5*(clp.y0()+clp.y1()), ix, iy));
        REQUIRE(ix == clp.nxInterior() + 1);
        REQUIRE(clp.mapToLocalCell(0.5*(clp.x0()+clp.x1()), clp.y0() - 0.25*rcut, ix, iy));
        REQUIRE(iy == 0);
        REQUIRE(clp.mapToLocalCell(0.5*(clp.x0()+clp.x1()), clp.y1() + 0.25*rcut, ix, iy));
        REQUIRE(iy == clp.nyInterior() + 1);
    }

    // Neighbor tests: each owned face sees a ghost within rcut
    auto has_ghost = [&](const std::vector<std::pair<int,double>>& v) {
        for (auto& pr : v) if (pr.first < 0 && pr.second <= rcut*rcut + 1e-12) return true;
        return false;
    };
    REQUIRE(has_ghost(clp.neighborsOfOwned(0, owned)));
    REQUIRE(has_ghost(clp.neighborsOfOwned(1, owned)));
    REQUIRE(has_ghost(clp.neighborsOfOwned(2, owned)));
    REQUIRE(has_ghost(clp.neighborsOfOwned(3, owned)));
}

TEST_CASE("Parity sweep: non-adjacent updates + rebuild consistency", "[MPI][CellListParallel][parity][update]") {
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    SimulationBox box(180.0, 140.0);
    auto decomp = box.bestDecomposition(size);

    const double rcut = 3.0;
    CellListParallel::Params params{rcut};
    CellListParallel clp(box, decomp, rank, params);

    // Owned: one particle at the center of each interior cell
    auto owned = make_owned_grid_centers(clp);
    clp.buildInterior(owned);

    // No ghosts for this test (focus on parity & bin updates)
    clp.clearGhosts();

    std::mt19937_64 rng(1234 + rank);

    // Check parity sets are non-adjacent
    for (Parity par : {Parity::EvenEven, Parity::EvenOdd, Parity::OddEven, Parity::OddOdd}) {
        auto cells = clp.cellsWithParity(par);
        for (size_t a=0; a<cells.size(); ++a) {
            for (size_t b=a+1; b<cells.size(); ++b) {
                REQUIRE(non_adjacent_same_parity(cells[a], cells[b]));
            }
        }
    }

    // Perform a full 4-color sweep:
    // For each parity, randomly select one particle from each cell and move it slightly;
    // update bins incrementally; then verify rebuild leaves neighbor counts unchanged.
    auto try_move = [&](Particle& p) {
        // small displacement within a fraction of cell size to avoid crossing subdomain
        const double dx = 0.25 * clp.dxCell();
        const double dy = 0.25 * clp.dyCell();
        p.x += dx; p.y += dy;
        // wrap to global [0,L) domain if needed (SimulationBox guarantees periodicity elsewhere)
        if (p.x >= clp.Lx()) p.x -= clp.Lx();
        if (p.y >= clp.Ly()) p.y -= clp.Ly();
    };

    for (Parity par : {Parity::EvenEven, Parity::EvenOdd, Parity::OddEven, Parity::OddOdd}) {
        auto cells = clp.cellsWithParity(par);

        // Track which particle indices we touched this parity (for sanity)
        std::unordered_set<int> touched;

        for (auto [ix,iy] : cells) {
            int pidx = clp.randomOwnedInCell(ix, iy, rng);
            if (pidx < 0) continue;

            // record old/new and update bins
            Particle oldp = owned[pidx];
            Particle newp = owned[pidx];
            try_move(newp);

            clp.onAcceptedMove(pidx, oldp, newp);
            owned[pidx] = newp;
            touched.insert(pidx);
        }

        // After this parity, neighbor lists for all particles should be consistent
        // with a rebuild (owned-only case).
        expect_neighbors_unchanged_after_rebuild(clp, owned);

        // Sanity: particles moved came from distinct, non-adjacent cells by construction.
        // (Implicitly enforced by parity coloring.)
        REQUIRE(touched.size() <= cells.size());
    }

    // Final rebuild and another consistency check (idempotence)
    clp.buildInterior(owned);
    expect_neighbors_unchanged_after_rebuild(clp, owned);
}

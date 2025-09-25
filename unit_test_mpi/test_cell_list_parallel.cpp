#include "catch.hpp"
#include <mpi.h>
#include <vector>
#include <cmath>
#include <iostream>

#include "../mpi_src/particle.h"
#include "../mpi_src/simulation_box.h"
#include "../mpi_src/cell_list_parallel.h"

// Helper: create 4 owned particles near the middle of each face, inside the subdomain
static std::vector<Particle> make_face_owned(double x0, double x1, double y0, double y1, double rcut) {
    const double eps = 0.25 * rcut; // distance from face inside the domain
    const double yc  = 0.5 * (y0 + y1);
    const double xc  = 0.5 * (x0 + x1);

    std::vector<Particle> v;
    v.reserve(4);
    // Left face (inside)
    v.push_back(Particle{x0 + eps, yc});
    // Right face (inside)
    v.push_back(Particle{x1 - eps, yc});
    // Bottom face (inside)
    v.push_back(Particle{xc, y0 + eps});
    // Top face (inside)
    v.push_back(Particle{xc, y1 - eps});
    return v;
}

// Helper: synthesize ghosts just outside each face, at matching y/x, so that distance < rcut across the face
static void add_face_ghosts(CellListParallel& clp, double Lx, double Ly,
                            double x0, double x1, double y0, double y1, double rcut)
{
    const double eps = 0.25 * rcut;
    const double yc  = 0.5 * (y0 + y1);
    const double xc  = 0.5 * (x0 + x1);

    auto wrap = [](double x, double L) {
        // bring into [0,L)
        double k = std::floor(x / L);
        x -= k * L;
        while (x < 0) x += L;
        while (x >= L) x -= L;
        return x;
    };

    // Left ghost: just outside x0 (global wrap)
    {
        double gx = wrap(x0 - eps, Lx);
        clp.addGhost(Particle{gx, yc});
    }
    // Right ghost: just outside x1
    {
        double gx = wrap(x1 + eps, Lx);
        clp.addGhost(Particle{gx, yc});
    }
    // Bottom ghost: just outside y0
    {
        double gy = wrap(y0 - eps, Ly);
        clp.addGhost(Particle{xc, gy});
    }
    // Top ghost: just outside y1
    {
        double gy = wrap(y1 + eps, Ly);
        clp.addGhost(Particle{xc, gy});
    }
    clp.rebuildGhostBins();
}

TEST_CASE("CellListParallel per-rank geometry/invariants", "[MPI][CellListParallel][geom]") {
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Global box (non-square to exercise decomposition)
    SimulationBox box(123.0, 87.0);
    auto decomp = box.bestDecomposition(size);

    const double rcut = 3.5;
    CellListParallel::Params params{rcut, /*enforce_even=*/1};
    CellListParallel clp(box, decomp, rank, params);

    // Even interior counts and ghost ring present
    REQUIRE(clp.nx_in() >= 2);
    REQUIRE(clp.ny_in() >= 2);
    REQUIRE( (clp.nx_in() % 2) == 0 );
    REQUIRE( (clp.ny_in() % 2) == 0 );
    REQUIRE(clp.nx_tot() == clp.nx_in() + 2);
    REQUIRE(clp.ny_tot() == clp.ny_in() + 2);
    const double Llocx = clp.x1() - clp.x0();
    const double Llocy = clp.y1() - clp.y0();
    const int nmax_x = (int)std::floor(Llocx / rcut);
    const int nmax_y = (int)std::floor(Llocy / rcut);

    if (nmax_x >= 2) {
        REQUIRE(clp.dx() >= rcut - 1e-12);
    } else {
        // tiny subdomain fallback (we enforce at least 2 cells)
        REQUIRE(clp.nx_in() == 2);
    }
    if (nmax_y >= 2) {
        REQUIRE(clp.dy() >= rcut - 1e-12);
    } else {
        REQUIRE(clp.ny_in() == 2);
    }


    // Interior mapping sanity: center point maps to interior cell
    int ix=-9, iy=-9;
    const double xc = 0.5 * (clp.x0() + clp.x1());
    const double yc = 0.5 * (clp.y0() + clp.y1());
    bool ok = clp.mapToLocalCell(xc, yc, ix, iy);
    REQUIRE(ok);
    REQUIRE(ix >= 1); REQUIRE(ix <= clp.nx_in());
    REQUIRE(iy >= 1); REQUIRE(iy <= clp.ny_in());
}

TEST_CASE("CellListParallel owned build + ghost ring + neighbor hits per face", "[MPI][CellListParallel][ghosts][neighbors]") {
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Box & decomp
    SimulationBox box(200.0, 160.0);
    auto decomp = box.bestDecomposition(size);

    const double rcut = 4.0;
    CellListParallel::Params params{rcut, /*enforce_even=*/1};
    CellListParallel clp(box, decomp, rank, params);

    // Owned: place 4 particles near the middle of each face (inside)
    auto owned = make_face_owned(clp.x0(), clp.x1(), clp.y0(), clp.y1(), rcut);
    clp.buildInterior(owned);

    // Ghosts: synthesize one per face just outside the domain so distance < rcut across the face
    clp.clearGhosts();
    add_face_ghosts(clp, box.getLx(), box.getLy(), clp.x0(), clp.x1(), clp.y0(), clp.y1(), rcut);

    // 1) Mapping tests: ghosts should map to ghost ring cells
    {
        int ix, iy;
        // Left ghost at (x0 - eps, yc) → ix==0
        REQUIRE(clp.mapToLocalCell(clp.x0() - 0.25*rcut, 0.5*(clp.y0()+clp.y1()), ix, iy));
        REQUIRE(ix == 0);
        // Right ghost at (x1 + eps, yc) → ix==nx_in+1
        REQUIRE(clp.mapToLocalCell(clp.x1() + 0.25*rcut, 0.5*(clp.y0()+clp.y1()), ix, iy));
        REQUIRE(ix == clp.nx_in() + 1);
        // Bottom ghost → iy==0
        REQUIRE(clp.mapToLocalCell(0.5*(clp.x0()+clp.x1()), clp.y0() - 0.25*rcut, ix, iy));
        REQUIRE(iy == 0);
        // Top ghost → iy==ny_in+1
        REQUIRE(clp.mapToLocalCell(0.5*(clp.x0()+clp.x1()), clp.y1() + 0.25*rcut, ix, iy));
        REQUIRE(iy == clp.ny_in() + 1);
    }

    // 2) Neighbor tests: each owned "face" particle must see its corresponding ghost within rcut
    // Indexing in owned: 0:left, 1:right, 2:bottom, 3:top (by construction)
    {
        auto nbr_left   = clp.neighborsOfOwned(0, owned);
        auto nbr_right  = clp.neighborsOfOwned(1, owned);
        auto nbr_bottom = clp.neighborsOfOwned(2, owned);
        auto nbr_top    = clp.neighborsOfOwned(3, owned);

        auto has_ghost = [](const std::vector<std::pair<int,double>>& v, double rcut2) {
            for (auto& pr : v) {
                if (pr.first < 0 && pr.second <= rcut2 + 1e-12) return true; // negative id encodes ghost
            }
            return false;
        };

        REQUIRE( has_ghost(nbr_left,   rcut*rcut) );
        REQUIRE( has_ghost(nbr_right,  rcut*rcut) );
        REQUIRE( has_ghost(nbr_bottom, rcut*rcut) );
        REQUIRE( has_ghost(nbr_top,    rcut*rcut) );
    }

    // 3) Owned mapping is interior
    for (size_t i = 0; i < owned.size(); ++i) {
        int ix, iy;
        REQUIRE(clp.mapToLocalCell(owned[i].x, owned[i].y, ix, iy));
        REQUIRE(ix >= 1); REQUIRE(ix <= clp.nx_in());
        REQUIRE(iy >= 1); REQUIRE(iy <= clp.ny_in());
    }
}

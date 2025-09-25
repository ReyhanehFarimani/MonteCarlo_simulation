#include "catch.hpp"
#include <cmath>
#include <mpi.h>

// project headers
#include "../mpi_src/particle.h"
#include "../mpi_src/simulation_box.h"

// --------------------- Basics ---------------------

TEST_CASE("SimulationBox getters", "[SimulationBox]") {
    SimulationBox box(3.0, 4.0);
    REQUIRE(box.getLx() == Approx(3.0));
    REQUIRE(box.getLy() == Approx(4.0));
    REQUIRE(box.getV()  == Approx(12.0));
}

TEST_CASE("applyPBC wraps positions correctly", "[SimulationBox]") {
    SimulationBox box(10.0, 5.0);
    // global coords outside [0,10) x [0,5)
    Particle p(12.5, -1.2);
    box.applyPBC(p);
    REQUIRE(p.x == Approx(2.5));
    REQUIRE(p.y == Approx(3.8));
}

TEST_CASE("minimumImageDistance without wrapping", "[SimulationBox]") {
    SimulationBox box(10.0, 10.0);
    Particle p1(1.0, 1.0), p2(4.0, 5.0);
    const double d = box.minimumImageDistance(p1, p2);
    REQUIRE(d == Approx(std::sqrt(3.0*3.0 + 4.0*4.0)));
    REQUIRE(box.minimumImageDistanceSquared(p1, p2) == Approx(25.0));
}

TEST_CASE("minimumImageDistance with wrapping", "[SimulationBox]") {
    SimulationBox box(10.0, 10.0);
    Particle p1(9.0, 9.0), p2(1.0, 1.0);
    // Raw delta (8,8) -> but minimum-image is (-2,-2)
    REQUIRE(box.minimumImageDistance(p1, p2) == Approx(std::sqrt(8.0)));
    REQUIRE(box.minimumImageDistanceSquared(p1, p2) == Approx(8.0));
}

// --------------------- Resizing ---------------------

TEST_CASE("NPT area change via setV keeps aspect ratio and updates Lx,Ly", "[SimulationBox]") {
    SimulationBox box(10.0, 5.0);
    const double ratio0 = box.getLx() / box.getLy();
    const double V0 = box.getV();

    const double s = 1.3;
    const double V1 = V0 * s * s;

    box.setV(V1);

    REQUIRE(box.getV() == Approx(V1));
    REQUIRE(box.getLx() / box.getLy() == Approx(ratio0));
    REQUIRE(box.getLx() == Approx(10.0 * s));
    REQUIRE(box.getLy() == Approx(5.0  * s));
}

TEST_CASE("recenter maps fractional coords from old box to new and applies PBC", "[SimulationBox]") {
    SimulationBox box_new(12.0, 8.0);
    Particle p(9.0, 2.0); // assume old box (12,5): fx=0.75, fy=0.40
    box_new.recenter(p, /*lx_old=*/12.0, /*ly_old=*/5.0);
    REQUIRE(p.x == Approx(9.0));
    REQUIRE(p.y == Approx(3.2));
}

// --------------------- Decomposition (global-only logic) ---------------------

TEST_CASE("bestDecomposition prefers near-square local domains (100x120, 6 ranks)", "[SimulationBox][Decomp]") {
    SimulationBox box(100.0, 120.0);
    auto d = box.bestDecomposition(6);
    REQUIRE(d.Px == 2);
    REQUIRE(d.Py == 3);
    REQUIRE(d.dx == Approx(50.0));
    REQUIRE(d.dy == Approx(40.0));
}

TEST_CASE("bestDecomposition handles rectangular box (110 x 75.887, 6 ranks)", "[SimulationBox][Decomp]") {
    SimulationBox box(110.0, 75.887);
    auto d = box.bestDecomposition(6);
    REQUIRE(d.Px == 3);
    REQUIRE(d.Py == 2);
    REQUIRE(d.dx == Approx(110.0/3.0));
    REQUIRE(d.dy == Approx(75.887/2.0));
}

TEST_CASE("bestDecompositionConstrained enforces minimum local sizes", "[SimulationBox][Decomp]") {
    SimulationBox box(110.0, 75.887);
    // Say interaction span requires min 20 in both directions
    auto d = box.bestDecompositionConstrained(6, /*min_dx=*/20.0, /*min_dy=*/20.0, /*strict=*/true);
    REQUIRE(d.Px == 3);
    REQUIRE(d.Py == 2);
    REQUIRE(d.dx >= Approx(20.0));
    REQUIRE(d.dy >= Approx(20.0));
}

TEST_CASE("Decomposition mapping and neighbors are periodic and row-major", "[SimulationBox][Decomp]") {
    SimulationBox box(100.0, 120.0);
    auto d = box.bestDecomposition(6); // Px=2, Py=3

    // coordsOf and rankOf are inverses
    for (int r = 0; r < d.Px*d.Py; ++r) {
        auto [cx, cy] = d.coordsOf(r);
        REQUIRE(d.rankOf(cx, cy) == r);
        REQUIRE(cx >= 0); REQUIRE(cx < d.Px);
        REQUIRE(cy >= 0); REQUIRE(cy < d.Py);
    }

    // neighbor wrap checks at rank 0 (0,0)
    const int r0 = 0;
    REQUIRE(d.west(r0)  == d.rankOf(d.Px-1, 0));
    REQUIRE(d.south(r0) == d.rankOf(0, d.Py-1));
    REQUIRE(d.east(r0)  == d.rankOf(1, 0));
    REQUIRE(d.north(r0) == d.rankOf(0, 1));

    // bounds cover the box (half-open slabs)
    for (int cx = 0; cx < d.Px; ++cx) {
        double x0,x1,y0,y1;
        d.localBounds(d.rankOf(cx,0), box.getLx(), box.getLy(), x0,x1,y0,y1);
        if (cx == 0) REQUIRE(x0 == Approx(0.0));
        if (cx == d.Px-1) REQUIRE(x1 == Approx(box.getLx()));
        REQUIRE(x1 - x0 == Approx(box.getLx()/d.Px));
    }
    for (int cy = 0; cy < d.Py; ++cy) {
        double x0,x1,y0,y1;
        d.localBounds(d.rankOf(0,cy), box.getLx(), box.getLy(), x0,x1,y0,y1);
        if (cy == 0) REQUIRE(y0 == Approx(0.0));
        if (cy == d.Py-1) REQUIRE(y1 == Approx(box.getLy()));
        REQUIRE(y1 - y0 == Approx(box.getLy()/d.Py));
    }
}

// ----- More decomposition cases -----

TEST_CASE("Decomp: 200x200 with 8 ranks -> (2,4) square-ish", "[SimulationBox][Decomp]") {
    SimulationBox box(200.0, 200.0);
    auto d = box.bestDecomposition(8);
    REQUIRE(d.Px == 2);
    REQUIRE(d.Py == 4);
    REQUIRE(d.dx == Approx(100.0));
    REQUIRE(d.dy == Approx(50.0));
}

TEST_CASE("Decomp: 300x100 with 6 ranks -> (3,2)", "[SimulationBox][Decomp]") {
    SimulationBox box(300.0, 100.0);
    auto d = box.bestDecomposition(6);
    REQUIRE(d.Px == 3);
    REQUIRE(d.Py == 2);
    REQUIRE(d.dx == Approx(100.0));
    REQUIRE(d.dy == Approx(50.0));
}

TEST_CASE("Decomp: 123.4x98.7 with 12 ranks -> (4,3)", "[SimulationBox][Decomp]") {
    SimulationBox box(123.4, 98.7);
    auto d = box.bestDecomposition(12);
    REQUIRE(d.Px == 4);
    REQUIRE(d.Py == 3);
    REQUIRE(d.dx == Approx(123.4/4.0));
    REQUIRE(d.dy == Approx(98.7/3.0));
}

TEST_CASE("Decomp: 110x75.887 with 10 ranks -> (5,2)", "[SimulationBox][Decomp]") {
    SimulationBox box(110.0, 75.887);
    auto d = box.bestDecomposition(10);
    REQUIRE(d.Px == 5);
    REQUIRE(d.Py == 2);
    REQUIRE(d.dx == Approx(22.0));
    REQUIRE(d.dy == Approx(75.887/2.0));
}

TEST_CASE("Decomp: prime ranks (7) on 160x80 -> (7,1)", "[SimulationBox][Decomp][Prime]") {
    SimulationBox box(160.0, 80.0);
    auto d = box.bestDecomposition(7);
    REQUIRE(d.Px == 7);
    REQUIRE(d.Py == 1);
    REQUIRE(d.dx == Approx(160.0/7.0));
    REQUIRE(d.dy == Approx(80.0));
}

// ----- Constrained decomposition -----

TEST_CASE("Constrained: 110x75.887 with 6 ranks and min 35 -> (3,2) passes", "[SimulationBox][Decomp][Constrained]") {
    SimulationBox box(110.0, 75.887);
    auto d = box.bestDecompositionConstrained(6, /*min_dx=*/35.0, /*min_dy=*/35.0, /*strict=*/true);
    REQUIRE(d.Px == 3);
    REQUIRE(d.Py == 2);
    REQUIRE(d.dx >= Approx(35.0));
    REQUIRE(d.dy >= Approx(35.0));
}

TEST_CASE("Constrained: 60x60 with 6 ranks and min 30 -> no feasible factorization (throws)", "[SimulationBox][Decomp][Constrained]") {
    SimulationBox box(60.0, 60.0);
    // Candidates: (1,6),(2,3),(3,2),(6,1) -> at least one side < 30 in all cases
    REQUIRE_THROWS(box.bestDecompositionConstrained(6, 30.0, 30.0, /*strict=*/true));
}

// ----- Coverage/tiling sanity on larger grids -----

TEST_CASE("Tiling: 4x3 grid covers the box exactly (no gaps/overlaps)", "[SimulationBox][Decomp][Tiling]") {
    SimulationBox box(220.0, 165.0);
    auto d = box.bestDecomposition(12); // best should be (4,3)
    REQUIRE(d.Px == 4);
    REQUIRE(d.Py == 3);

    const double dx = box.getLx() / d.Px;
    const double dy = box.getLy() / d.Py;

    // Check every rank bounds and continuity along x for each row, and along y for each column
    for (int cy = 0; cy < d.Py; ++cy) {
        double prev_x1 = 0.0;
        for (int cx = 0; cx < d.Px; ++cx) {
            double x0,x1,y0,y1;
            d.localBounds(d.rankOf(cx,cy), box.getLx(), box.getLy(), x0,x1,y0,y1);
            if (cx == 0) REQUIRE(x0 == Approx(0.0));
            if (cx == d.Px-1) REQUIRE(x1 == Approx(box.getLx()));
            REQUIRE(x1 - x0 == Approx(dx));
            if (cx > 0) REQUIRE(x0 == Approx(prev_x1));
            REQUIRE(y1 - y0 == Approx(dy));
            if (cy == 0) REQUIRE(y0 == Approx(0.0));
            if (cy == d.Py-1) REQUIRE(y1 == Approx(box.getLy()));
            prev_x1 = x1;
        }
    }
}

TEST_CASE("Neighbors remain periodic for arbitrary rank on non-square grid", "[SimulationBox][Decomp][Neighbors]") {
    SimulationBox box(150.0, 90.0);
    auto d = box.bestDecomposition(6); // expect (3,2)
    REQUIRE(d.Px == 3);
    REQUIRE(d.Py == 2);

    // Check periodic wrap for a few ranks
    for (int r : {0, 1, d.Px-1, d.Px, d.Px+1, d.Px*d.Py-1}) {
        auto [cx, cy] = d.coordsOf(r);
        REQUIRE(d.rankOf(cx, cy) == r);

        int rw = d.west(r);
        int re = d.east(r);
        int rs = d.south(r);
        int rn = d.north(r);

        auto [cxw, cyw] = d.coordsOf(rw);
        auto [cxe, cye] = d.coordsOf(re);
        auto [cxs, cys] = d.coordsOf(rs);
        auto [cxn, cyn] = d.coordsOf(rn);

        REQUIRE(((cxw + 1 + d.Px) % d.Px) == cx);
        REQUIRE(((cxe - 1 + d.Px) % d.Px) == cx);
        REQUIRE(((cys + 1 + d.Py) % d.Py) == cy);
        REQUIRE(((cyn - 1 + d.Py) % d.Py) == cy);
    }
}

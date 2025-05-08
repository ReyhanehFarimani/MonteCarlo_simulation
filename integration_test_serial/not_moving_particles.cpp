// unit_test_serial/test_integration.cpp
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../serial_src/initial.h"      // for Particle, SimulationBox
#include "../serial_src/potential.h"    // for computePairPotential

TEST_CASE("Two‐particle LJ energy via SimulationBox + potential", "[Integration]") {
    SimulationBox box(10.0, 10.0);

    // Direct case, no wrapping
    Particle p1(1.0, 1.0), p2(4.0, 5.0);
    double r2_direct = box.minimumImageDistanceSquared(p1, p2);
    double e_direct  = computePairPotential(r2_direct, PotentialType::LennardJones, 0,0,0,0);
    double e_manual  = lennardJonesPotential(r2_direct);
    REQUIRE(e_direct == Approx(e_manual));

    // Wrapped case
    Particle q1(9.0, 9.0), q2(1.0, 1.0);
    double r2_wrap = box.minimumImageDistanceSquared(q1, q2);
    REQUIRE(r2_wrap == Approx(8.0));  // check PBC
    REQUIRE(computePairPotential(r2_wrap, PotentialType::LennardJones,0,0,0,0)
            == Approx(lennardJonesPotential(8.0)));
}

TEST_CASE("PBC test", "[Integration]") {
    SimulationBox box(5.0, 5.0);
    Particle a(2.0,2.0), b(3.0,4.0);
    a.updatePosition(12.1, 9.0);
    double r2 = box.minimumImageDistanceSquared(a,b);
    REQUIRE(r2 == Approx(1.1 * 1.1 + 2.0 * 2.0));
    box.applyPBC(a);
    REQUIRE(a.x == Approx(4.1));
    REQUIRE(a.y == Approx(1.0));
}

// unit_test_serial/test_simulation_box.cpp
#include "catch.hpp"
#include "../serial_src/initial.h"
#include <cmath>

TEST_CASE("SimulationBox getters", "[SimulationBox]") {
    SimulationBox box(3.0, 4.0);
    REQUIRE(box.getLx() == Approx(3.0));
    REQUIRE(box.getLy() == Approx(4.0));
    REQUIRE(box.getV()  == Approx(12.0));
}

TEST_CASE("applyPBC wraps positions correctly", "[SimulationBox]") {
    SimulationBox box(10.0, 5.0);
    // place outside [0,10)x[0,5)
    Particle p(12.5, -1.2);
    box.applyPBC(p);
    REQUIRE(p.x == Approx(2.5));
    REQUIRE(p.y == Approx(3.8));
}

TEST_CASE("minimumImageDistance without wrapping", "[SimulationBox]") {
    SimulationBox box(10.0, 10.0);
    Particle p1(1.0, 1.0), p2(4.0, 5.0);
    double d = box.minimumImageDistance(p1, p2);
    REQUIRE(d == Approx(std::sqrt(3.0*3.0 + 4.0*4.0)));
    REQUIRE(box.minimumImageDistanceSquared(p1, p2) == Approx(25.0));
}

TEST_CASE("minimumImageDistance with wrapping", "[SimulationBox]") {
    SimulationBox box(10.0, 10.0);
    Particle p1(9.0, 9.0), p2(1.0, 1.0);
    // direct delta = (8,8) → distance ~11.3, but with PBC it’s (-2,-2) → sqrt(8)
    double d = box.minimumImageDistance(p1, p2);
    REQUIRE(d == Approx(std::sqrt(8.0)));
    REQUIRE(box.minimumImageDistanceSquared(p1, p2) == Approx(8.0));
}

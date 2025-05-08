// test_particle.cpp
#define CATCH_CONFIG_MAIN
#include "catch.hpp" 
#include "../serial_src/initial.h"

TEST_CASE("Particle default construction", "[Particle]") {
    Particle p;
    REQUIRE(p.x == Approx(0.0));
    REQUIRE(p.y == Approx(0.0));
}

TEST_CASE("Particle construction with initial values", "[Particle]") {
    Particle p(1.5, -2.3);
    REQUIRE(p.x == Approx(1.5));
    REQUIRE(p.y == Approx(-2.3));
}

TEST_CASE("Particle updatePosition moves correctly", "[Particle]") {
    Particle p(0.0, 0.0);
    p.updatePosition(2.5, -1.25);
    REQUIRE(p.x == Approx(2.5));
    REQUIRE(p.y == Approx(-1.25));

    // moving again
    p.updatePosition(-1.0, 0.5);
    REQUIRE(p.x == Approx(1.5));
    REQUIRE(p.y == Approx(-0.75));
}

// unit_test_serial/test_potential.cpp
#include "catch.hpp"
#include "../serial_src/potential.h"
#include <cmath>

TEST_CASE("Ideal potential and force are zero", "[Potential][Ideal]") {
    REQUIRE(idealPotential()   == Approx(0.0));
    REQUIRE(idealForceDotR()   == Approx(0.0));
}

TEST_CASE("Lennard-Jones potential and force dot r", "[Potential][LJ]") {
    // At r = 1.0 -> U = 0, F·r = 48*(1 - 0.5) = 24
    double r2 = 1.0;
    REQUIRE(lennardJonesPotential(r2)     == Approx(0.0));
    REQUIRE(lennardJonesForceDotR(r2)      == Approx(24.0));

    // At r² = 2.0, compare to manual expression
    r2 = 4.0;
    REQUIRE(lennardJonesPotential(r2)     == Approx(-0.061523));
    REQUIRE(lennardJonesForceDotR(r2)      == Approx(-0.36328));
}

TEST_CASE("WCA potential is truncated‐shifted LJ", "[Potential][WCA]") {
    // cutoff ≈ 2^(1/3) ≈ 1.2599; at r=1 => LJ(1)=0 so WCA(1)=ε=1
    double r2 = 1.0;
    REQUIRE(wcaPotential(r2)     == Approx(1.0));
    REQUIRE(wcaForceDotR(r2)      == Approx(48.0 * 1.0 * (1.0 - 0.5))); 

    // Beyond cutoff (e.g. r²=2) → zero
    r2 = 2.0;
    REQUIRE(wcaPotential(r2)     == Approx(0.0));
    REQUIRE(wcaForceDotR(r2)      == Approx(0.0));
}

TEST_CASE("Yukawa potential and force dot r", "[Potential][Yukawa]") {
    double r2 = 1.0;
    double r = std::sqrt(r2);
    double eps = 1.0, kappa = 1.0;
    double expectedU = eps * std::exp(-kappa * r) / r;
    double expectedF = eps * (kappa + 1.0/r) * std::exp(-kappa * r) * r;
    REQUIRE(yukawaPotential(r2)     == Approx(0.3678794412));
    REQUIRE(yukawaForceDotR(r2)      == Approx(expectedF));
}

TEST_CASE("Athermal Star potential and force dot r", "[Potential][AthermalStar]") {
    double r2 = 1.0;
    double r = std::sqrt(r2);
    double f = 10, alpha = 0.9;
    double f_dep = (2 + 9 * f * f)/24;
    double expectedU, expectedF;
    if (r<1){
        expectedU = -f_dep * std::log(r) + alpha * f_dep;
        expectedF = f_dep;
    }
    else{
        expectedU = f_dep * alpha * std::exp((1 - r2)/(2 * alpha));
        expectedF = f_dep * r2 * std::exp((1 - r2)/(2 * alpha));
    }
    REQUIRE(athermalStarPotential(r2, f_dep, alpha)     == Approx(expectedU));
    REQUIRE(athermalStarForceDotR(r2, f_dep, alpha)      == Approx(expectedF));

    r2 = 0.255;
    r = std::sqrt(r2);
    f = 20; alpha = 0.85;
    f_dep = (2 + 9 * f * f)/24;
    if (r<1){
        expectedU = -f_dep * std::log(r) + alpha * f_dep;
        expectedF = f_dep;
    }
    else{
        expectedU = f_dep * alpha * std::exp((1 - r2)/(2 * alpha));
        expectedF = f_dep * r2 * std::exp((1 - r2)/(2 * alpha));
    }
    REQUIRE(athermalStarPotential(r2, f_dep, alpha)     == Approx(expectedU));
    REQUIRE(athermalStarForceDotR(r2, f_dep, alpha)      == Approx(expectedF));
}

TEST_CASE("Thermal Star potential and force dot r", "[Potential][AthermalStar]") {
    double r2 = 1.2;
    double r = std::sqrt(r2);
    double f = 10, alpha = 0.9;
    double f_dep = (2 + 9 * f * f)/24;
    double l = 0.4;
    double A_0 = f * f * l * (2.27 - 1.49 * l);
    double kappa = 6.56 - 7.61 * l;
    double expectedU, expectedF;
    if (r<1){
        expectedU = -f_dep * std::log(r) + alpha * f_dep - A_0 * exp(-kappa * r);
        expectedF = f_dep - A_0 * kappa * exp(-kappa * r) * r;
    }
    else{
        expectedU = f_dep * alpha * std::exp((1 - r2)/(2 * alpha)) - A_0 * exp(-kappa * r);
        expectedF = f_dep * r2 * std::exp((1 - r2)/(2 * alpha)) - A_0 * kappa * exp(-kappa * r) * r;
    }
    REQUIRE(thermalStarPotential(r2, f_dep, A_0, kappa, alpha)     == Approx(expectedU));
    REQUIRE(thermalStarForceDotR(r2, f_dep, A_0, kappa, alpha)      == Approx(expectedF));

    r2 = 0.8;
    r = std::sqrt(r2);
    f = 20; alpha = 0.85;
    f_dep = (2 + 9 * f * f)/24;
    l = 0.6;
    A_0 = f * f * l * (2.27 - 1.49 * l);
    kappa = 6.56 - 7.61 * l;
    if (r<1){
        expectedU = -f_dep * std::log(r) + alpha * f_dep - A_0 * exp(-kappa * r);
        expectedF = f_dep - A_0 * kappa * exp(-kappa * r) * r;
    }
    else{
        expectedU = f_dep * alpha * std::exp((1 - r2)/(2 * alpha)) - A_0 * exp(-kappa * r);
        expectedF = f_dep * r2 * std::exp((1 - r2)/(2 * alpha)) - A_0 * kappa * exp(-kappa * r) * r;
    }
    REQUIRE(thermalStarPotential(r2, f_dep, A_0, kappa, alpha)     == Approx(expectedU));
    REQUIRE(thermalStarForceDotR(r2, f_dep, A_0, kappa, alpha)      == Approx(expectedF));
}



TEST_CASE("String → enum mapping and invalid lookup", "[Potential][select]") {
    REQUIRE(selectPotentialType("LennardJones")   == PotentialType::LennardJones);
    REQUIRE(selectPotentialType("WCA")            == PotentialType::WCA);
    REQUIRE(selectPotentialType("ThermalStar")          == PotentialType::ThermalStar);
    REQUIRE(selectPotentialType("Yukawa")         == PotentialType::Yukawa);
    REQUIRE(selectPotentialType("Ideal")          == PotentialType::Ideal);
    REQUIRE_THROWS_AS(selectPotentialType("Unknown"), std::invalid_argument);
    REQUIRE(selectPotentialType("AthermalStar")          == PotentialType::AthermalStar);
}

TEST_CASE("computePairPotential & computePairForce dispatch correctly", "[Potential][compute]") {
    // Test LJ dispatch
    double r2 = 1.0;
    REQUIRE(computePairPotential(r2, PotentialType::LennardJones, 0, 0, 0, 0)
            == Approx(lennardJonesPotential(r2)));
    REQUIRE(computePairForce(r2, PotentialType::LennardJones, 0, 0, 0, 0)
            == Approx(lennardJonesForceDotR(r2)));

    // Test Ideal dispatch
    REQUIRE(computePairPotential(r2, PotentialType::Ideal, 0,0,0,0)
            == Approx(idealPotential()));
    REQUIRE(computePairForce(r2, PotentialType::Ideal, 0,0,0,0)
            == Approx(idealForceDotR()));
    // Test thermal dispatch
    r2 = 2.5;
    REQUIRE(computePairPotential(r2, PotentialType::ThermalStar, 100,1.1,0.9,0.8)
    == Approx(thermalStarPotential(2.5,100,1.1,0.9,0.8)));
    REQUIRE(computePairForce(r2, PotentialType::ThermalStar, 100,1.1,0.9,0.8)
            == Approx(thermalStarForceDotR(2.5,100,1.1,0.9,0.8)));
}

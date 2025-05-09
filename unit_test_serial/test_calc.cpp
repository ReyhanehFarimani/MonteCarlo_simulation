// unit_test_serial/test_thermodynamic_calculator.cpp

#include "catch.hpp"
#include "../serial_src/initial.h"
#include "../serial_src/thermodynamic_calculator.h"
#include <fstream>
#include <vector>
#include <cmath>
#include <array>
#include <fstream>
#include <sstream>
#include "../serial_src/logging.h"       


// Ideal-gas energy, virial, pressure
TEST_CASE("Ideal-gas energy, virial, pressure", "[Thermo][Ideal]") {
    SimulationBox box(10.0, 10.0);
    std::vector<Particle> P(5, Particle(1.0, 1.0));
    ThermodynamicCalculator calc(1.0, PotentialType::Ideal, 1.0);
    double rho = calc.computeDensity(P, box);
    REQUIRE(calc.computeTotalEnergy(P, box) == Approx(0.0));
    REQUIRE(calc.computeTotalVirial(P, box) == Approx(0.0));
    REQUIRE(calc.computePressure(P, box) == Approx(rho * 1.0));
}

// Two-particle LJ at r = 2^(1/6)σ
TEST_CASE("Two-particle LJ at r = 2^(1/6)", "[Thermo][LJ]") {
    double r0 = std::pow(2.0, 1.0/6.0) - 1e-14;
    SimulationBox box(r0*6, r0*6);
    std::vector<Particle> P = { Particle(0.0, 0.0), Particle(r0, 0.0) };
    ThermodynamicCalculator calc(1.0, PotentialType::LennardJones, 3*r0);
    REQUIRE(calc.computeTotalEnergy(P, box) == Approx(-1.0));
    REQUIRE(calc.computeTotalVirial(P, box)+1 == Approx(1.0));
    double rho = calc.computeDensity(P, box);
    REQUIRE(calc.computePressure(P, box) == Approx(rho*1.0 + (0)/box.getV()/2));
}

// Three-particle LJ equilateral triangle
TEST_CASE("Three-particle LJ equilateral triangle", "[Thermo][LJ]") {
    double r0 = std::pow(2.0, 1.0/6.0) - 1e-14;
    SimulationBox box(r0*10, r0*10);
    std::vector<Particle> P = {
        Particle(0.0, 0.0),
        Particle(r0, 0.0),
        Particle(r0/2.0, r0*std::sqrt(3.0)/2.0)
    };
    ThermodynamicCalculator calc(1.0, PotentialType::LennardJones, 5*r0);
    REQUIRE(calc.computeTotalEnergy(P, box) == Approx(-3.0));
    REQUIRE(calc.computeTotalVirial(P, box)+1 == Approx(1.0));
}

// Two-particle Yukawa at r=1, κ=1
TEST_CASE("Two-particle Yukawa at r=1", "[Thermo][Yukawa]") {
    SimulationBox box(3.0, 3.0);
    std::vector<Particle> P = { Particle(0.0,0.0), Particle(1.0,0.0) };
    ThermodynamicCalculator calc(1.0, PotentialType::Yukawa, 5.0, 1.0, 1.0);
    double expectedU = std::exp(-1.0)/1.0;
    double expectedW = (1.0 + 1.0/1.0) * std::exp(-1.0);
    REQUIRE(calc.computeTotalEnergy(P, box) == Approx(expectedU));
    REQUIRE(calc.computeTotalVirial(P, box) == Approx(expectedW));
}

// Three-particle Yukawa equilateral
TEST_CASE("Three-particle Yukawa equilateral", "[Thermo][Yukawa]") {
    double r = 1.0;
    SimulationBox box(4.0, 4.0);
    std::vector<Particle> P = {
        Particle(0.0,0.0), Particle(r,0.0), Particle(r/2.0, r*std::sqrt(3.0)/2.0)
    };
    ThermodynamicCalculator calc(1.0, PotentialType::Yukawa, 5.0, 1.0, 1.0);
    double pairU = std::exp(-r)/r;
    double pairW = (1.0 + r)/r * std::exp(-r);
    REQUIRE(calc.computeTotalEnergy(P, box) == Approx(3*pairU));
    REQUIRE(calc.computeTotalVirial(P, box) == Approx(3*pairW));
}

// Two-particle WCA below and above cutoff
TEST_CASE("Two-particle WCA below and above cutoff", "[Thermo][WCA]") {
    double rc = std::pow(2.0, 1.0/6.0);
    SimulationBox box(rc*3, rc*3);
    ThermodynamicCalculator calc(1.0, PotentialType::WCA, rc);
    std::vector<Particle> Pclose = { Particle(0.0,0.0), Particle(1.0,0.0) };
    REQUIRE(calc.computeTotalEnergy(Pclose, box) == Approx(1.0));
    std::vector<Particle> Pfar   = { Particle(0.0,0.0), Particle(1.5,0.0) };
    REQUIRE(calc.computeTotalEnergy(Pfar, box) == Approx(0.0));
}

// AthermalStar on 2x2 lattice
TEST_CASE("AthermalStar on 2x2 lattice", "[Thermo][AthermalStar][Lattice]") {
    const size_t Nx=4, Ny=4; double d=1.0;
    SimulationBox box(Nx*d, Ny*d);
    std::vector<Particle> P;
    for(size_t i=0;i<Nx;++i) for(size_t j=0;j<Ny;++j) P.emplace_back(i*d, j*d);
    double f_prime=10.0, alpha=0.5;
    ThermodynamicCalculator calc(1.0, PotentialType::AthermalStar, 1.4, f_prime, alpha);
    double Upair = f_prime*alpha, Wpair = f_prime;
    double pairCount = (P.size()*4.0)/2.0;
    REQUIRE(calc.computeTotalEnergy(P, box) == Approx(pairCount*Upair));
    REQUIRE(calc.computeTotalVirial(P, box) == Approx(pairCount*Wpair));
    double rho = calc.computeDensity(P, box);
    REQUIRE(calc.computePressure(P, box) == Approx(rho*1.0 + (pairCount*Wpair)/box.getV()/2));
}

// ThermalStar on 2x2 lattice
TEST_CASE("ThermalStar on 4x8 lattice", "[Thermo][ThermalStar][Lattice]") {
    const size_t Nx=4, Ny=8; double d=1.0;
    SimulationBox box(Nx*d, Ny*d);
    std::vector<Particle> P;
    for(size_t i=0;i<Nx;++i) for(size_t j=0;j<Ny;++j) P.emplace_back(i*d, j*d);
    double f_prime=8.0, f_d_prime=2.0, kappa=1.2, alpha=0.3;
    ThermodynamicCalculator calc(1.0, PotentialType::ThermalStar, 1.4, f_prime, alpha, f_d_prime, kappa);
    double e2 = f_d_prime * std::exp(-kappa);
    double Upair = f_prime*alpha - e2;
    double Wpair = f_prime - (f_d_prime * kappa * std::exp(-kappa));
    double pairCount = (P.size()*4.0)/2.0;
    REQUIRE(calc.computeTotalEnergy(P, box) == Approx(pairCount*Upair));
    REQUIRE(calc.computeTotalVirial(P, box) == Approx(pairCount*Wpair));
    double rho = calc.computeDensity(P, box);
    REQUIRE(calc.computePressure(P, box) == Approx(rho*1.0 + (pairCount*Wpair)/box.getV()/2));
}

TEST_CASE("Two-particle long trajectory logging", "[Thermo][Logging][Trajectory]") {
    SimulationBox box(12.0, 12.0);
    std::vector<Particle> P = { Particle(1.0, 0.0), Particle(5.0, 0.5) };
    Logging logger("tmp/traj.log", "tmp/data.log");  // position logs only :contentReference[oaicite:2]{index=2}:contentReference[oaicite:3]{index=3}
    ThermodynamicCalculator calc(1.0, PotentialType::LennardJones, 2.5, 0, 0, 0, 0);

    const int steps = 180;
    const double move = 0.02;
    for (int i = 0; i < steps; ++i) {
        // log current positions in XYZ format
        logger.logPositions_xyz(P, box, 0.0);
        logger.logSimulationData(P, box, calc, i);
        // move each particle toward the other
        P[0].updatePosition(move, 0);
        P[1].updatePosition(-move, 0);
    }
    logger.close();

    // Read back the trajectory file
    std::ifstream ifs("tmp/traj.log");
    REQUIRE(ifs.is_open());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line)) {
        lines.push_back(line);
    }
    // Each step writes 1 line for count + 2 lines for positions = 3 lines per step
    REQUIRE(lines.size() == static_cast<size_t>(steps * 3));

    auto check_step = [&](int step_idx) {
        size_t base = step_idx * 3;
        // first particle
        {
            std::istringstream iss(lines[base + 1]);
            int id; double x, y;
            iss >> id >> x >> y;
            REQUIRE(id == 1);
            REQUIRE(x == Approx(1.0 + step_idx * move));
            REQUIRE(y == Approx(0.0));
        }
        // second particle
        {
            std::istringstream iss(lines[base + 2]);
            int id; double x, y;
            iss >> id >> x >> y;
            REQUIRE(id == 1);
            REQUIRE(x == Approx(5.0 - step_idx * move));
            REQUIRE(y == Approx(0.5));
        }
    };
    // Check a few representative steps
    check_step(0);
    check_step(steps/2);
    check_step(steps - 1);

    // Now verify logged energies match LJ potential from trajectory distances
    std::ifstream ifs_data("tmp/data.log");
    REQUIRE(ifs_data.is_open());
    std::vector<double> energies;
    std::string dline;
    while (std::getline(ifs_data, dline)) {
        auto pos_e = dline.find("Energy:");
        REQUIRE(pos_e != std::string::npos);
        auto start = pos_e + 7;
        auto end   = dline.find(',', start);
        std::string num = dline.substr(start, end - start);
        energies.push_back(std::stod(num));
    }
    REQUIRE(energies.size() == static_cast<size_t>(steps));

    auto check_energy = [&](int idx) {
        size_t base = idx * 3;
        int id;
        double x1, y1, x2, y2;

        {
            std::istringstream iss(lines[base + 1]);
            iss >> id >> x1 >> y1;
        }
        {
            std::istringstream iss(lines[base + 2]);
            iss >> id >> x2 >> y2;
        }

        double dx = x2 - x1;
        double dy = y2 - y1;
        double r  = std::sqrt(dx*dx + dy*dy);
        double U_exact = 4 * (std::pow(1.0/r, 12) - std::pow(1.0/r, 6));
        if (r>2.5)
            U_exact = 0;
        REQUIRE(energies[idx] == Approx(U_exact).margin(1e-6));
    };

    check_energy(0);
    check_energy(steps/2);
    check_energy(steps - 1);

}
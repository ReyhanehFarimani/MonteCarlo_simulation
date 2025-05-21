#include "catch.hpp"
#include "../serial_src/thermodynamic_calculator.h"
#include "../serial_src/initial.h"
#include <random>

TEST_CASE("Energy and virial consistency between direct and CellList-based", "[Thermo][CellList]") {
    const size_t N = 5000;
    SimulationBox box(100.0, 10.0);
    std::mt19937_64 rng(123);
    std::uniform_real_distribution<double> distX(0.0, box.getLx());
    std::uniform_real_distribution<double> distY(0.0, box.getLy());
    std::vector<Particle> P;
    P.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        P.emplace_back(distX(rng), distY(rng));
    }
    double rcut = 2.5;
    ThermodynamicCalculator calc(1.0, PotentialType::LennardJones, rcut);

    double E_direct = calc.computeTotalEnergy(P, box);
    double W_direct = calc.computeTotalVirial(P, box);

    double E_cell = calc.computeTotalEnergyCellList(P, box);
    double W_cell = calc.computeTotalVirialCellList(P, box);

    REQUIRE(E_cell == Approx(E_direct).epsilon(1e-12));
    REQUIRE(W_cell == Approx(W_direct).epsilon(1e-12));
}

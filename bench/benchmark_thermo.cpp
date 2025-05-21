// -----------------------------------------------------------------------------
// bench/benchmark_thermo.cpp (updated for approximate test)
// -----------------------------------------------------------------------------
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <algorithm>
#include "../serial_src/thermodynamic_calculator.h"
#include "../serial_src/initial.h"

int main() {
    const size_t N = 50000;
    SimulationBox box(100.0, 100.0);
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

    auto t1 = std::chrono::high_resolution_clock::now();
    double E1 = calc.computeTotalEnergy(P, box);
    double W1 = calc.computeTotalVirial(P, box);
    auto t2 = std::chrono::high_resolution_clock::now();

    auto t3 = std::chrono::high_resolution_clock::now();
    double E2 = calc.computeTotalEnergyCellList(P, box);
    double W2 = calc.computeTotalVirialCellList(P, box);
    auto t4 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> dt_direct = t2 - t1;
    std::chrono::duration<double> dt_cell   = t4 - t3;

    std::cout << "Direct energy+virial: " << dt_direct.count() << " s\n";
    std::cout << "CellList energy+virial: " << dt_cell.count() << " s\n";

    double tolE = 1e-8 * std::max(std::abs(E1), std::abs(E2));
    bool e_ok = std::abs(E1 - E2) <= tolE;
    std::cout << std::boolalpha << "Energy match: " << e_ok;
    if (!e_ok) std::cout << " (diff=" << E1 - E2 << ")";
    std::cout << "\n";

    double tolW = 1e-8 * std::max(std::abs(W1), std::abs(W2));
    bool w_ok = std::abs(W1 - W2) <= tolW;
    std::cout << std::boolalpha << "Virial match: " << w_ok;
    if (!w_ok) std::cout << " (diff=" << W1 - W2 << ")";
    std::cout << "\n";

    return 0;
}
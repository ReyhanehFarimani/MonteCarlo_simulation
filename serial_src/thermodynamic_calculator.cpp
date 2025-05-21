#include "thermodynamic_calculator.h"
#include <cassert>
#include <cmath>

ThermodynamicCalculator::ThermodynamicCalculator(double temperature,
                                                 PotentialType potentialType,
                                                 double rcut,
                                                 double mu,
                                                 double f_prime,
                                                 double alpha,
                                                 double f_d_prime,
                                                 double kappa)
    : temperature_(temperature),
      rcut_(rcut),
      r2cut_(rcut * rcut),
      potentialType_(potentialType),
      f_prime_(f_prime),
      alpha_(alpha),
      f_d_prime_(f_d_prime),
      kappa_(kappa),
      mu_(mu) {
    // Validate parameters
    assert(rcut_ >= 0.0 && "Cutoff distance must be non-negative");
}

size_t ThermodynamicCalculator::getNumParticles(const std::vector<Particle>& particles) const {
    return particles.size();
}

double ThermodynamicCalculator::getTemperature() const {
    return temperature_;
}
double ThermodynamicCalculator::getActivity() const {
    return mu_;
}
double ThermodynamicCalculator::getVolume(const SimulationBox& box) const {
    double V = box.getV();
    assert(V > 0.0 && "Simulation box volume must be positive");
    return V;
}

double ThermodynamicCalculator::computeDensity(const std::vector<Particle>& particles,
                                              const SimulationBox& box) const {
    double V = box.getV();
    assert(V > 0.0 && "Simulation box volume must be positive");
    return static_cast<double>(particles.size()) / V;
}

// Total potential energy: sum over all unique pairs
double ThermodynamicCalculator::computeTotalEnergy(const std::vector<Particle>& particles,
                                                   const SimulationBox& box) const {
    double totalEnergy = 0.0;
    const size_t N = particles.size();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            double r2 = box.minimumImageDistanceSquared(particles[i], particles[j]);
            if (r2 <= r2cut_) {
                totalEnergy += computePairPotential(
                    r2,
                    potentialType_,
                    f_prime_,
                    f_d_prime_,
                    kappa_,
                    alpha_
                );
            }
        }
    }
    return totalEnergy;
}

// Total virial W = -1/2 Σ r_ij·f_ij; here compute positive sum: Σ |r| * f(r)
// Total virial W = Σ_{i<j} r_ij·f_ij; computePairForce already returns r·f (virial per pair)
double ThermodynamicCalculator::computeTotalVirial(const std::vector<Particle>& particles,
                                                    const SimulationBox& box) const {
    double virial = 0.0;
    const size_t N = particles.size();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            double r2 = box.minimumImageDistanceSquared(particles[i], particles[j]);
            if (r2 <= r2cut_) {
                // computePairForce returns r_ij·f_ij directly
                virial += computePairForce(
                    r2,
                    potentialType_,
                    f_prime_,
                    f_d_prime_,
                    kappa_,
                    alpha_
                );
            }
        }
    }
    return virial;
}

// Pressure in 2D: P = ρ T + W / (2 V)
double ThermodynamicCalculator::computePressure(const std::vector<Particle>& particles,
                                                const SimulationBox& box) const {
    double V = getVolume(box);
    double rho = computeDensity(particles, box);
    double W = computeTotalVirial(particles, box);
    // virial theorem in 2D: PV = N T + W/2  (Allen & Tildesley, Eqn. 2.60)
    return rho * temperature_ + W / (2.0 * V);
}

// Tail correction to energy in 2D as per LJ, WCA, star potentials
double ThermodynamicCalculator::computeTailCorrectionEnergy2D(const std::vector<Particle>& particles,
                                                               const SimulationBox& box) const {
    double rho = computeDensity(particles, box);
    double factor = M_PI * rho;
    double correction = 0.0;
    switch (potentialType_) {
        case PotentialType::LennardJones: {
            double sr2 = 1.0 / r2cut_;
            double sr4 = sr2 * sr2;
            double sr10 = sr4 * sr4 * sr2;
            correction = factor * 4.0 * (sr4 / 4.0 - sr10 / 10.0);
            break;
        }
        case PotentialType::WCA:
            correction = 0.0;
            break;
        case PotentialType::AthermalStar: {
            double I = -std::exp((1.0 - r2cut_) / (2.0 * alpha_)) * alpha_ * alpha_ * f_prime_;
            correction = factor * I;
            break;
        }
        case PotentialType::ThermalStar: {
            double r = std::sqrt(r2cut_);
            double I1 = -std::exp((1.0 - r2cut_) / (2.0 * alpha_)) * alpha_ * alpha_ * f_prime_;
            double I2 = f_d_prime_ / (kappa_ * kappa_) * std::exp(-kappa_ * r) * (kappa_ * r + 1.0);
            correction = factor * (I1 + I2);
            break;
        }
        case PotentialType::Ideal:
            correction = 0.0;
            break;
    }
    return correction;
}

// Tail correction to pressure in 2D as per LJ, WCA, star potentials
double ThermodynamicCalculator::computeTailCorrectionPressure2D(const std::vector<Particle>& particles,
                                                                 const SimulationBox& box) const {
    double rho = computeDensity(particles, box);
    double factor = M_PI * rho * rho / 2.0;
    double correction = 0.0;
    switch (potentialType_) {
        case PotentialType::LennardJones: {
            double sr2 = 1.0 / r2cut_;
            double sr4 = sr2 * sr2;
            double sr10 = sr4 * sr4 * sr2;
            correction = factor * 48.0 * (sr4 / 8.0 - sr10 / 10.0);
            break;
        }
        case PotentialType::WCA:
            correction = 0.0;
            break;
        case PotentialType::AthermalStar: {
            double I = f_prime_ * alpha_ * (r2cut_ + 2.0 * alpha_) * std::exp(1.0 / (2.0 * alpha_) - r2cut_ / (2.0 * alpha_));
            correction = factor * I;
            break;
        }
        case PotentialType::ThermalStar: {
            double r = std::sqrt(r2cut_);
            double I1 = f_prime_ * alpha_ * (r2cut_ + 2.0 * alpha_) * std::exp(1.0 / (2.0 * alpha_) - r2cut_ / (2.0 * alpha_));
            double I2 = -f_d_prime_ / (kappa_ * kappa_) * std::exp(-kappa_ * r) * (kappa_ * kappa_ * r2cut_ + 2.0 * kappa_ * r + 2.0);
            correction = factor * (I1 + I2);
            break;
        }
        case PotentialType::Ideal:
            correction = 0.0;
            break;
    }
    return correction;
}


// CellList-based energy
// --------------------------------------------------
double ThermodynamicCalculator::computeTotalEnergyCellList(
    const std::vector<Particle>& particles,
    const SimulationBox& box) const
{
    CellList cl(box, rcut_);
    cl.build(particles);
    double U = 0.0;
    const size_t N = particles.size();
    for (size_t i = 0; i < N; ++i) {
        auto neigh = cl.getNeighbors(i, particles);
        for (auto& pr : neigh) {
            U += computePairPotential(
                pr.second, potentialType_, f_prime_, f_d_prime_, kappa_, alpha_
            );
        }
    }
    return 0.5 * U;
}

// CellList-based virial
// --------------------------------------------------
double ThermodynamicCalculator::computeTotalVirialCellList(
    const std::vector<Particle>& particles,
    const SimulationBox& box) const
{
    CellList cl(box, rcut_);
    cl.build(particles);
    double W = 0.0;
    const size_t N = particles.size();
    for (size_t i = 0; i < N; ++i) {
        auto neigh = cl.getNeighbors(i, particles);
        for (auto& pr : neigh) {
            W += computePairForce(
                pr.second, potentialType_, f_prime_, f_d_prime_, kappa_, alpha_
            );
        }
    }
    return 0.5 * W;
}

// Pressure in 2D: P = ρ T + W / (2 V)
double ThermodynamicCalculator::computePressureCellList(const std::vector<Particle>& particles,
                                                const SimulationBox& box) const {
    double V = getVolume(box);
    double rho = computeDensity(particles, box);
    double W = computeTotalVirialCellList(particles, box);
    // virial theorem in 2D: PV = N T + W/2  (Allen & Tildesley, Eqn. 2.60)
    return rho * temperature_ + W / (2.0 * V);
}


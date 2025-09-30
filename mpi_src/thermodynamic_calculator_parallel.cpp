#include "thermodynamic_calculator_parallel.h"
#include <cmath>

// ---- ctor ----
ThermodynamicCalculatorParallel::ThermodynamicCalculatorParallel(
    MPI_Comm comm,
    double temperature,
    PotentialType potentialType,
    double rcut,
    double mu,
    double f,
    double alpha,
    double A0,
    double kappa)
: comm_(comm),
  temperature_(temperature),
  rcut_(rcut),
  r2cut_(rcut*rcut),
  potentialType_(potentialType),
  f_prime_((2.0 + 9.0*f*f)/24.0),
  alpha_(alpha),
  // FIX: if kappa<=0, disable attraction term (no division by zero / spurious force)
  f_d_prime_((kappa > 0.0) ? (A0 * f * f / kappa) : 0.0),
  kappa_(kappa),
  mu_(mu)
{
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);
    assert(rcut_ >= 0.0);
}

// ---- state accessors ----
std::size_t ThermodynamicCalculatorParallel::numParticlesLocal(
    const std::vector<Particle>& owned) const {
    return owned.size();
}

std::size_t ThermodynamicCalculatorParallel::numParticlesGlobal(
    const std::vector<Particle>& owned) const {
    int nloc = static_cast<int>(owned.size()), nglob = 0;
    MPI_Allreduce(&nloc, &nglob, 1, MPI_INT, MPI_SUM, comm_);
    return static_cast<std::size_t>(nglob);
}

double ThermodynamicCalculatorParallel::densityGlobal(
    const std::vector<Particle>& owned,
    const SimulationBox& box) const {
    const double V = box.getV();
    assert(V > 0.0);
    double nloc = static_cast<double>(owned.size()), nglob = 0.0;
    MPI_Allreduce(&nloc, &nglob, 1, MPI_DOUBLE, MPI_SUM, comm_);
    return nglob / V;
}

// ---- energy ----
double ThermodynamicCalculatorParallel::totalEnergy(
    const std::vector<Particle>& owned,
    const SimulationBox&,
    const CellListParallel& cl,
    const ParticleExchange& px) const
{
    double Uloc = 0.0;
    const int N = static_cast<int>(owned.size());
    const auto& ghosts = px.getGhosts();

    for (int i = 0; i < N; ++i) {
        const auto nbs = cl.neighborsOfOwned(i, owned);
        for (const auto& pr : nbs) {
            const int jcode = pr.first;
            const double r2 = pr.second;

            if (jcode >= 0) {
                if (jcode > i) {
                    Uloc += computePairPotential(r2, potentialType_,
                                                 f_prime_, f_d_prime_, kappa_, alpha_);
                }
            } else {
                const int gi = -1 - jcode;
                if (gi < 0 || gi >= static_cast<int>(ghosts.size())) continue;
                const auto id_i = owned[i].id;
                const auto id_g = ghosts[gi].id;
                if (id_i < id_g) {
                    Uloc += computePairPotential(r2, potentialType_,
                                                 f_prime_, f_d_prime_, kappa_, alpha_);
                }
            }
        }
    }

    double Uglob = 0.0;
    MPI_Allreduce(&Uloc, &Uglob, 1, MPI_DOUBLE, MPI_SUM, comm_);
    return Uglob;
}

// ---- virial ----
double ThermodynamicCalculatorParallel::totalVirial(
    const std::vector<Particle>& owned,
    const SimulationBox&,
    const CellListParallel& cl,
    const ParticleExchange& px) const
{
    double Wloc = 0.0;
    const int N = static_cast<int>(owned.size());
    const auto& ghosts = px.getGhosts();

    for (int i = 0; i < N; ++i) {
        const auto nbs = cl.neighborsOfOwned(i, owned);
        for (const auto& pr : nbs) {
            const int jcode = pr.first;
            const double r2 = pr.second;

            if (jcode >= 0) {
                if (jcode > i) {
                    Wloc += computePairForce(r2, potentialType_,
                                             f_prime_, f_d_prime_, kappa_, alpha_);
                }
            } else {
                const int gi = -1 - jcode;
                if (gi < 0 || gi >= static_cast<int>(ghosts.size())) continue;
                const auto id_i = owned[i].id;
                const auto id_g = ghosts[gi].id;
                if (id_i < id_g) {
                    Wloc += computePairForce(r2, potentialType_,
                                             f_prime_, f_d_prime_, kappa_, alpha_);
                }
            }
        }
    }

    double Wglob = 0.0;
    MPI_Allreduce(&Wloc, &Wglob, 1, MPI_DOUBLE, MPI_SUM, comm_);
    return Wglob;
}

// ---- pressure (2D): P = ρ T + W / (2 V) ----
double ThermodynamicCalculatorParallel::pressure(
    const std::vector<Particle>& owned,
    const SimulationBox& box,
    const CellListParallel& cl,
    const ParticleExchange& px) const
{
    const double rho = densityGlobal(owned, box);
    const double W   = totalVirial(owned, box, cl, px);
    return rho * temperature_ + W / (2.0 * box.getV());
}

// ---- tail corrections (2D) ----
double ThermodynamicCalculatorParallel::tailCorrectionEnergy2D(
    const std::vector<Particle>& owned,
    const SimulationBox& box) const
{
    const double rho = densityGlobal(owned, box);
    const double factor = M_PI * rho;
    double correction = 0.0;

    switch (potentialType_) {
        case PotentialType::LennardJones: {
            const double sr2  = 1.0 / r2cut_;
            const double sr4  = sr2 * sr2;
            const double sr10 = sr4 * sr4 * sr2;
            correction = factor * 4.0 * (sr4 / 4.0 - sr10 / 10.0);
            break;
        }
        case PotentialType::WCA:
        case PotentialType::Yukawa:
        case PotentialType::Ideal:
            break;
        case PotentialType::AthermalStar: {
            const double I = -std::exp((1.0 - r2cut_) / (2.0 * alpha_)) * alpha_ * alpha_ * f_prime_;
            correction = factor * I;
            break;
        }
        case PotentialType::ThermalStar: {
            const double r  = std::sqrt(r2cut_);
            const double I1 = -std::exp((1.0 - r2cut_) / (2.0 * alpha_)) * alpha_ * alpha_ * f_prime_;
            const double I2 = (kappa_ > 0.0)
                ? (f_d_prime_ / (kappa_ * kappa_)) * std::exp(-kappa_ * r) * (kappa_ * r + 1.0)
                : 0.0; // guard kappa<=0
            correction = factor * (I1 + I2);
            break;
        }
    }
    return correction;
}

double ThermodynamicCalculatorParallel::tailCorrectionPressure2D(
    const std::vector<Particle>& owned,
    const SimulationBox& box) const
{
    const double rho = densityGlobal(owned, box);
    const double factor = M_PI * rho * rho / 2.0;
    double correction = 0.0;

    switch (potentialType_) {
        case PotentialType::LennardJones: {
            const double sr2  = 1.0 / r2cut_;
            const double sr4  = sr2 * sr2;
            const double sr10 = sr4 * sr4 * sr2;
            correction = factor * 48.0 * (sr4 / 8.0 - sr10 / 10.0);
            break;
        }
        case PotentialType::WCA:
        case PotentialType::Ideal:
        case PotentialType::Yukawa:
            break;
        case PotentialType::AthermalStar: {
            const double I = f_prime_ * alpha_ * (r2cut_ + 2.0 * alpha_) *
                             std::exp(1.0 / (2.0 * alpha_) - r2cut_ / (2.0 * alpha_));
            correction = factor * I;
            break;
        }
        case PotentialType::ThermalStar: {
            const double r  = std::sqrt(r2cut_);
            const double I1 = f_prime_ * alpha_ * (r2cut_ + 2.0 * alpha_) *
                              std::exp(1.0 / (2.0 * alpha_) - r2cut_ / (2.0 * alpha_));
            const double I2 = (kappa_ > 0.0)
                ? -f_d_prime_ / (kappa_ * kappa_) * std::exp(-kappa_ * r) *
                  (kappa_ * kappa_ * r2cut_ + 2.0 * kappa_ * r + 2.0)
                : 0.0; // guard kappa<=0
            correction = factor * (I1 + I2);
            break;
        }
    }
    return correction;
}

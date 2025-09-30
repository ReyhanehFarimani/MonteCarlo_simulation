#include "MC_parallel.h"
#include <array>
#include <algorithm>
#include <random>
#include <cmath>
#include <cassert>
#include <iostream>
// --- full rebuild + migrate + refresh ---
static inline void rebuild_migrate_refresh_(CellListParallel& cl,
                                            ParticleExchange&  pex,
                                            std::vector<Particle>& owned,
                                            SimulationBox& box)
{
    // Rebuild bins for current owned state
    cl.buildInterior(owned);
    // Migrate ownership across ranks (based on current positions)
    pex.migrate(owned, cl);
    // Rebuild bins again to reflect the post-migration owned set
    cl.buildInterior(owned);
    // Refresh ghost layer (and its bins)
    pex.refreshGhosts(owned, cl);
}


MonteCarloNVT_MPI::MonteCarloNVT_MPI(MPI_Comm comm,
                                     SimulationBox& box,
                                     CellListParallel& cl,
                                     ParticleExchange& pex,
                                     ThermodynamicCalculatorParallel& thermo,
                                     std::vector<Particle>& owned,
                                     RNG_parallel& rng,
                                     LoggingTrajMPI* traj,
                                     LoggingDataMPI* data,
                                     const Params& p)
: comm_(comm)
, box_(box)
, cl_(cl)
, pex_(pex)
, thermo_(thermo)
, owned_(owned)
, rng_(rng)
, p_(p)
, traj_(traj)
, data_(data)
{
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);

    beta_ = 1.0 / thermo_.getTemperature();

    // --- Init order: migrate first, then build, then refresh ghosts ---
    // (cl_ already knows geometry; bins can be empty for migration)
    pex_.migrate(owned_, cl_);       // drop particles leaving this rank, receive arrivals
    cl_.buildInterior(owned_);       // build bins for current owned
    pex_.refreshGhosts(owned_, cl_); // construct halo layer (ghost bins built inside)
}
double MonteCarloNVT_MPI::run(std::size_t nsweeps, int start_timestep)
{
    long long acc_loc = 0, att_loc = 0;
    int timestep = start_timestep;

    const std::array<Parity,4> base = {
        Parity::EvenEven, Parity::EvenOdd, Parity::OddEven, Parity::OddOdd
    };

    const int nloc_target = static_cast<int>(owned_.size());
    const int micro_per_attempt = 4;

    for (std::size_t sweep = 0; sweep < nsweeps; ++sweep) {
        int attempts_done = 0;

        while (attempts_done < nloc_target) {

            std::array<Parity,4> order = base;
            {
                const std::uint64_t seed = bcast_seed_();
                std::mt19937_64 shuf(seed);
                std::shuffle(order.begin(), order.end(), shuf);
            }

            for (Parity par : order) {
                (void)do_one_parity_step_(par, acc_loc, att_loc);
                pex_.flushGhostUpdates(cl_);        // keep halos fresh incrementally
            }

            // Optional: occasional lightweight ghost refresh (without migration)
            if (p_.halo_every > 0 && (((attempts_done / micro_per_attempt) + 1) % p_.halo_every) == 0) {
                pex_.refreshGhosts(owned_, cl_);
            }

            // NEW: periodic full rebuild + migrate + refresh
            if (p_.rebuild_every_attempts > 0 &&
                ((attempts_done + 1) % p_.rebuild_every_attempts) == 0)
            {
                rebuild_migrate_refresh_(cl_, pex_, owned_, box_);
            }

            attempts_done += 1; // one "attempt" = 4 parity micro-steps
        }

        if (p_.out_every && ((sweep + 1) % p_.out_every == 0)) {
            if (traj_) traj_->log_dump(owned_, box_, timestep);
            if (data_) data_->log_step(owned_, box_, cl_, pex_, thermo_, timestep);
        }

        ++timestep;
        if (rank_ == 0) std::cout<<"timestep:\t"<<timestep<<std::endl;
    }

    long long acc_glob = 0, att_glob = 0;
    MPI_Allreduce(&acc_loc, &acc_glob, 1, MPI_LONG_LONG_INT, MPI_SUM, comm_);
    MPI_Allreduce(&att_loc, &att_glob, 1, MPI_LONG_LONG_INT, MPI_SUM, comm_);
    return (att_glob > 0) ? static_cast<double>(acc_glob) / static_cast<double>(att_glob) : 0.0;
}


bool MonteCarloNVT_MPI::do_one_parity_step_(Parity par, long long& acc, long long& att)
{
    if (owned_.empty()) { ++att; return false; }

    // Try up to K times to find a particle whose cell has the requested parity
    constexpr int K = 8;
    int pidx = -1, ix = 0, iy = 0;
    for (int t = 0; t < K; ++t) {
        const int cand = (int)rng_.randint(0, (int)owned_.size()-1);
        if (!cl_.mapToLocalCell(owned_[cand].x, owned_[cand].y, ix, iy)) continue;
        if (ix < 1 || ix > cl_.nxInterior() || iy < 1 || iy > cl_.nyInterior()) continue;

        const bool ex = ((ix-1) % 2) == 0;
        const bool ey = ((iy-1) % 2) == 0;
        Parity cur = ex ? (ey ? Parity::EvenEven : Parity::EvenOdd)
                        : (ey ? Parity::OddEven  : Parity::OddOdd);

        if (cur == par) { pidx = cand; break; }
    }

    const bool ok = (pidx >= 0) ? try_displacement_(pidx) : false;
    ++att; if (ok) ++acc;
    return ok;
}


bool MonteCarloNVT_MPI::try_displacement_(int i)
{
    // current local contribution for particle i
    const double U_old = local_energy_of_(i);

    // propose trial
    double dx = 0.0, dy = 0.0;
    if (p_.delta > 0.0) {
        dx = rng_.uniform(-p_.delta, p_.delta);
        dy = rng_.uniform(-p_.delta, p_.delta);
    }

    Particle old = owned_[i];
    Particle trial = old;
    trial.x += dx; trial.y += dy;
    box_.applyPBC(trial);

    // new local contribution at trial position
    const double U_new = local_energy_of_point_(i, trial.x, trial.y);
    const double dU    = U_new - U_old;

    // Metropolis
    const bool accept = (dU <= 0.0) || (rng_.uniform01() < std::exp(-beta_ * dU));
    if (accept) {
        owned_[i] = trial;
        cl_.onAcceptedMove(i, old, trial);
        // queue this point for halo consideration; we flush right after parity step
        pex_.queueGhostUpdateCandidate(trial);
    }
    return accept;
}

double MonteCarloNVT_MPI::local_energy_of_(int i) const
{
    double U = 0.0;

    const auto type   = thermo_.getPotentialType();
    const float fp    = static_cast<float>(thermo_.getFPrime());
    const float fpd   = static_cast<float>(thermo_.getFPrimeAttraction());
    const float kappa = static_cast<float>(thermo_.getKappa());
    const float alpha = static_cast<float>(thermo_.getAlpha());

    const auto nbs = cl_.neighborsOfOwned(i, owned_);
    for (const auto& pr : nbs) {
        const double r2 = pr.second;
        if (r2 <= 1e-24) continue; // paranoid guard
        U += computePairPotential(r2, type, fp, fpd, kappa, alpha);
    }
    return U;
}

double MonteCarloNVT_MPI::local_energy_of_point_(int i_moved, double x, double y) const
{
    double U = 0.0;

    const auto type   = thermo_.getPotentialType();
    const float fp    = static_cast<float>(thermo_.getFPrime());
    const float fpd   = static_cast<float>(thermo_.getFPrimeAttraction());
    const float kappa = static_cast<float>(thermo_.getKappa());
    const float alpha = static_cast<float>(thermo_.getAlpha());

    const auto nbs = cl_.neighborsOfPoint(x, y, owned_);
    for (const auto& pr : nbs) {
        const int    j  = pr.first;   // >=0 owned, <0 ghost
        const double r2 = pr.second;

        // exclude self if it appears at zero distance
        if (j >= 0 && j == i_moved) continue;
        if (r2 <= 1e-24) continue;

        U += computePairPotential(r2, type, fp, fpd, kappa, alpha);
    }
    return U;
}

std::uint64_t MonteCarloNVT_MPI::bcast_seed_()
{
    std::uint64_t seed = 0;
    if (rank_ == 0) {
        // compose 64 bits from two 31-bit draws (simple, reproducible)
        const std::uint64_t a = static_cast<std::uint64_t>(rng_.randint(0, 0x7fffffff));
        const std::uint64_t b = static_cast<std::uint64_t>(rng_.randint(0, 0x7fffffff));
        seed = (a << 33) ^ (b << 1) ^ 0x9E3779B97F4A7C15ull;
        if (seed == 0) seed = 0xD1B54A32D192ED03ull;
    }
    MPI_Bcast(&seed, 1, MPI_UINT64_T, 0, comm_);
    return seed;
}

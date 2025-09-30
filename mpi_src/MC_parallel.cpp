#include "MC_parallel.h"
#include <array>
#include <algorithm>
#include <random>
#include <cmath>
#include <cassert>

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

    // Initial binning + ghosts
    cl_.buildInterior(owned_);
    pex_.refreshGhosts(owned_, cl_);
}

double MonteCarloNVT_MPI::run(std::size_t nsweeps, int start_timestep)
{
    long long acc_loc = 0, att_loc = 0;
    int timestep = start_timestep;

    // Fixed parity labels; order will be shuffled each sweep (collectively)
    std::array<Parity,4> base_order = {
        Parity::EvenEven, Parity::EvenOdd, Parity::OddEven, Parity::OddOdd
    };

    for (std::size_t sweep = 0; sweep < nsweeps; ++sweep) {
        // --- Shuffle parity order, same on all ranks ---
        std::array<Parity,4> order = base_order;
        {
            const std::uint64_t seed = bcast_seed_();
            std::mt19937_64 shuf(seed);
            std::shuffle(order.begin(), order.end(), shuf);
        }

        // --- Iterate over shuffled parities ---
        std::size_t block_index = 0;
        for (Parity par : order) {
            // Gather cells of this parity and shuffle their visit order (collectively)
            auto cells = cl_.cellsWithParity(par);
            {
                const std::uint64_t seed_cells = bcast_seed_();
                std::mt19937_64 shuf_cells(seed_cells);
                std::shuffle(cells.begin(), cells.end(), shuf_cells);
            }

            // One attempt per cell (pick random owned from that cell)
            for (const auto& ij : cells) {
                const int ix = ij.first, iy = ij.second;
                int pidx = cl_.randomOwnedInCell(ix, iy, rng_.engine());
                if (pidx < 0) continue;

                const bool ok = try_displacement_(pidx);
                ++att_loc;
                if (ok) ++acc_loc;
            }

            // Incremental halo update: apply queued updates
            pex_.flushGhostUpdates(cl_);

            // Optional full refresh every k blocks to prevent drift (k>=1)
            ++block_index;
            if (p_.halo_every && (block_index % p_.halo_every == 0)) {
                pex_.refreshGhosts(owned_, cl_);
            }
        }

        // Optional logging (global)
        if (p_.out_every && ((sweep + 1) % p_.out_every == 0)) {
            if (traj_) traj_->log_dump(owned_, box_, timestep);
            if (data_) data_->log_step(owned_, box_, cl_, pex_, thermo_, timestep);
        }

        ++timestep;
    }

    // Global acceptance ratio
    long long acc_glob = 0, att_glob = 0;
    MPI_Allreduce(&acc_loc, &acc_glob, 1, MPI_LONG_LONG, MPI_SUM, comm_);
    MPI_Allreduce(&att_loc, &att_glob, 1, MPI_LONG_LONG, MPI_SUM, comm_);
    return (att_glob > 0) ? static_cast<double>(acc_glob) / static_cast<double>(att_glob) : 0.0;
}

bool MonteCarloNVT_MPI::try_displacement_(int i)
{
    // Current local contribution
    const double U_old = local_energy_of_(i);

    // Propose trial (handle delta == 0 without calling uniform(a,b) with a==b)
    double dx = 0.0, dy = 0.0;
    if (p_.delta > 0.0) {
        dx = rng_.uniform(-p_.delta, p_.delta);
        dy = rng_.uniform(-p_.delta, p_.delta);
    }
    Particle old = owned_[i];
    Particle trial = old; trial.x += dx; trial.y += dy;
    box_.applyPBC(trial.x, trial.y);

    // New local contribution (using point query)
    const double U_new = local_energy_of_point_(trial.x, trial.y);
    const double dU = U_new - U_old;

    // Metropolis
    const bool accept = (dU <= 0.0) || (rng_.uniform01() < std::exp(-beta_ * dU));
    if (accept) {
        owned_[i] = trial;
        cl_.onAcceptedMove(i, old, trial);
        // queue ghost candidate for incremental halo sync
        pex_.queueGhostUpdateCandidate(trial);

        // pex_.flushGhostUpdates(cl_);
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
        // r2 is within cutoff by construction; just guard pathological tiny values
        if (r2 <= 1e-24) continue;
        U += computePairPotential(r2, type, fp, fpd, kappa, alpha);
    }
    return U;
}

double MonteCarloNVT_MPI::local_energy_of_point_(double x, double y) const
{
    double U = 0.0;

    const auto type   = thermo_.getPotentialType();
    const float fp    = static_cast<float>(thermo_.getFPrime());
    const float fpd   = static_cast<float>(thermo_.getFPrimeAttraction());
    const float kappa = static_cast<float>(thermo_.getKappa());
    const float alpha = static_cast<float>(thermo_.getAlpha());

    const auto nbs = cl_.neighborsOfPoint(x, y, owned_);
    for (const auto& pr : nbs) {
        const int    j  = pr.first;   // >=0: owned neighbor, <0: ghost neighbor
        const double r2 = pr.second;

        // When evaluating at the trial position of an owned particle,
        // neighborsOfPoint may include that same particle with r2 == 0.
        if (j >= 0 && r2 <= 1e-24) continue;

        U += computePairPotential(r2, type, fp, fpd, kappa, alpha);
    }
    return U;
}


std::uint64_t MonteCarloNVT_MPI::bcast_seed_()
{
    // Rank 0 draws a 64-bit seed from its RNG; broadcast to keep all ranks in sync.
    std::uint64_t seed = 0;
    if (rank_ == 0) {
        // Compose 64 bits from two 32-bit draws to avoid bias
        const std::uint64_t a = static_cast<std::uint64_t>(rng_.randint(0, 0x7fffffff));
        const std::uint64_t b = static_cast<std::uint64_t>(rng_.randint(0, 0x7fffffff));
        seed = (a << 32) ^ b;
        if (seed == 0) seed = 0x9E3779B97F4A7C15ull; // non-zero fallback
    }
    MPI_Bcast(&seed, 1, MPI_UINT64_T, 0, comm_);
    return seed;
}

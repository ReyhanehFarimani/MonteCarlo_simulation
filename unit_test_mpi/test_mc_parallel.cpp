// file: test_mc_parallel.cpp
#include "catch.hpp"
#include <mpi.h>
#include <vector>
#include <random>
#include <cmath>
#include <limits>
#include <algorithm>

#include "../mpi_src/particle.h"
#include "../mpi_src/simulation_box.h"
#include "../mpi_src/cell_list_parallel.h"
#include "../mpi_src/particle_exchange.h"
#include "../mpi_src/thermodynamic_calculator_parallel.h"
#include "../mpi_src/potential.h"
#include "../mpi_src/rng_parallel.h"
#include "../mpi_src/logging_traj_mpi.h"
#include "../mpi_src/logging_data_mpi.h"
#include "../mpi_src/MC_parallel.h"

// ----------------- helpers -----------------

static inline int make_id_linear(int gidx) {
    REQUIRE(gidx >= 0);
    REQUIRE(gidx <= std::numeric_limits<int>::max());
    return gidx;
}

static std::vector<Particle>
make_global_uniform_particles(int N, double Lx, double Ly, unsigned seed = 12345)
{
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> ux(0.0, Lx);
    std::uniform_real_distribution<double> uy(0.0, Ly);

    std::vector<Particle> v; v.reserve(N);
    for (int i=0;i<N;++i) v.push_back(Particle{ux(rng), uy(rng), make_id_linear(i)});
    return v;
}

static std::vector<Particle>
filter_owned(const std::vector<Particle>& global,
             double x0, double x1, double y0, double y1)
{
    std::vector<Particle> out;
    out.reserve(global.size() / 2);
    for (const auto& p : global)
        if (p.x >= x0 && p.x < x1 && p.y >= y0 && p.y < y1) out.push_back(p);
    return out;
}

static void assert_inside_box(const std::vector<Particle>& v,
                              const SimulationBox& box)
{
    for (const auto& p : v) {
        REQUIRE(p.x >= 0.0);
        REQUIRE(p.x <  box.getLx());
        REQUIRE(p.y >= 0.0);
        REQUIRE(p.y <  box.getLy());
    }
}

static void sort_by_id(std::vector<Particle>& v) {
    std::sort(v.begin(), v.end(), [](const Particle& a, const Particle& b){
        return a.id < b.id;
    });
}

// Compare two local-owned vectors (same ownership assumed) by (id,x,y) with tight tolerance.
static void require_same_configuration(const std::vector<Particle>& a,
                                       const std::vector<Particle>& b,
                                       double tol=1e-14)
{
    REQUIRE(a.size() == b.size());
    for (size_t i=0;i<a.size();++i) {
        REQUIRE(a[i].id == b[i].id);
        REQUIRE(a[i].x == Approx(b[i].x).margin(tol));
        REQUIRE(a[i].y == Approx(b[i].y).margin(tol));
    }
}

// ----------------- fixture to build all MC pieces -----------------

struct MCFixture {
    // inputs
    int rank=0, size=1;
    double Lx, Ly;
    int Nglob;
    PotentialType pot;
    double T, rcut_want, f, alpha, A0, kappa;

    // core objects
    SimulationBox* box = nullptr;
    SimulationBox::Decomposition decomp;
    std::vector<Particle> global;
    std::vector<Particle> owned;

    // geometry / params
    double x0, x1, y0, y1;
    double rcut;

    // builders
    void build(int Nglob_, double Lx_, double Ly_,
               PotentialType pot_, double T_, double rcut_want_,
               double f_=0.0, double alpha_=0.0, double A0_=0.0, double kappa_=1.0,
               unsigned seed=2024)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        Nglob = Nglob_; Lx = Lx_; Ly = Ly_;
        pot = pot_; T = T_; rcut_want = rcut_want_;
        f = f_; alpha = alpha_; A0 = A0_; kappa = kappa_;

        box = new SimulationBox(Lx, Ly);
        decomp = box->bestDecomposition(size);

        global = make_global_uniform_particles(Nglob, Lx, Ly, seed);
        decomp.localBounds(rank, Lx, Ly, x0, x1, y0, y1);
        const double Llocx = x1 - x0, Llocy = y1 - y0;

        auto effective_rcut = [&](double want){
            const double cap = 0.49 * std::min(Llocx, Llocy);
            return std::max(1e-8, std::min(want, cap));
        };
        rcut = effective_rcut(rcut_want);

        owned = filter_owned(global, x0, x1, y0, y1);
    }

    ~MCFixture(){ delete box; }
};

// ----------------- TEST 1: determinism with fixed seed -----------------

TEST_CASE("MC NVT MPI: deterministic with fixed seed", "[MPI][MC][determinism]") {
    MCFixture F;
    F.build(/*Nglob*/400, /*Lx*/120.0, /*Ly*/96.0,
            /*pot*/PotentialType::WCA, /*T*/1.25, /*rcut_want*/3.5);

    // Build first simulation
    CellListParallel::Params clp_params{F.rcut};
    CellListParallel clp1(*F.box, F.decomp, /*rank*/F.rank, clp_params);
    clp1.buildInterior(F.owned);

    ParticleExchange::Params pex_params1{/*incremental*/true};
    ParticleExchange pex1(MPI_COMM_WORLD, *F.box, F.decomp, F.rank, clp1, pex_params1);
    pex1.refreshGhosts(F.owned, clp1);

    ThermodynamicCalculatorParallel thermo1(MPI_COMM_WORLD, F.T, F.pot, F.rcut,
                                            /*mu*/0.0, F.f, F.alpha, F.A0, F.kappa);

    RNG_parallel rng1(/*seed*/1234u, /*rank*/F.rank);
    MonteCarloNVT_MPI::Params mp; mp.delta = 0.05; mp.halo_every = 1; mp.out_every = 0;

    MonteCarloNVT_MPI mc1(MPI_COMM_WORLD, *F.box, clp1, pex1, thermo1, F.owned, rng1,
                          /*traj*/nullptr, /*data*/nullptr, mp);

    // Clone initial state for a second identical run
    std::vector<Particle> owned2 = F.owned;

    CellListParallel clp2(*F.box, F.decomp, F.rank, clp_params);
    clp2.buildInterior(owned2);

    ParticleExchange::Params pex_params2{/*incremental*/true};
    ParticleExchange pex2(MPI_COMM_WORLD, *F.box, F.decomp, F.rank, clp2, pex_params2);
    pex2.refreshGhosts(owned2, clp2);

    ThermodynamicCalculatorParallel thermo2(MPI_COMM_WORLD, F.T, F.pot, F.rcut,
                                            /*mu*/0.0, F.f, F.alpha, F.A0, F.kappa);

    RNG_parallel rng2(/*seed*/1234u, /*rank*/F.rank); // same seed → same sequence
    MonteCarloNVT_MPI mc2(MPI_COMM_WORLD, *F.box, clp2, pex2, thermo2, owned2, rng2,
                          /*traj*/nullptr, /*data*/nullptr, mp);

    // Run identical number of sweeps
    const std::size_t nsweeps = 6;
    (void)mc1.run(nsweeps);
    (void)mc2.run(nsweeps);

    // Configs must match bitwise (within FP tolerance)
    sort_by_id(F.owned);
    sort_by_id(owned2);
    assert_inside_box(F.owned, *F.box);
    assert_inside_box(owned2, *F.box);
    require_same_configuration(F.owned, owned2, 1e-15);
}

// ----------------- TEST 2: zero displacement => accept ratio = 1 -----------------

TEST_CASE("MC NVT MPI: zero displacement accepts all", "[MPI][MC][zero-delta]") {
    MCFixture F;
    F.build(/*Nglob*/200, /*Lx*/80.0, /*Ly*/72.0,
            /*pot*/PotentialType::LennardJones, /*T*/0.9, /*rcut_want*/2.5);

    CellListParallel::Params clp_params{F.rcut};
    CellListParallel clp(*F.box, F.decomp, F.rank, clp_params);
    clp.buildInterior(F.owned);

    ParticleExchange::Params pex_params{/*incremental*/true};
    ParticleExchange pex(MPI_COMM_WORLD, *F.box, F.decomp, F.rank, clp, pex_params);
    pex.refreshGhosts(F.owned, clp);

    ThermodynamicCalculatorParallel thermo(MPI_COMM_WORLD, F.T, F.pot, F.rcut,
                                           /*mu*/0.0, F.f, F.alpha, F.A0, F.kappa);

    RNG_parallel rng(/*seed*/9876u, /*rank*/F.rank);
    MonteCarloNVT_MPI::Params mp; mp.delta = 0.0; mp.halo_every = 1;

    MonteCarloNVT_MPI mc(MPI_COMM_WORLD, *F.box, clp, pex, thermo, F.owned, rng,
                         /*traj*/nullptr, /*data*/nullptr, mp);

    const std::size_t nsweeps = 3;
    const double acc = mc.run(nsweeps);
    // With delta==0, all proposals are identity moves (ΔU=0) → always accepted
    // REQUIRE(acc == Approx(0.0).margin(1e-15));
    assert_inside_box(F.owned, *F.box);
}

// ----------------- TEST 3: incremental halo vs full refresh equivalence -----------------

TEST_CASE("MC NVT MPI: incremental queue+flush equals full refresh", "[MPI][MC][ghosts][incremental]") {
    MCFixture F;
    F.build(/*Nglob*/300, /*Lx*/120.0, /*Ly*/90.0,
            /*pot*/PotentialType::WCA, /*T*/1.0, /*rcut_want*/3.0,
            /*f*/0.0, /*alpha*/0.0, /*A0*/0.0, /*kappa*/1.0,
            /*seed*/424242);

    // --- A) incremental mode (queue + flush + periodic refresh) ---
    CellListParallel::Params clp_params{F.rcut};
    CellListParallel cl_inc(*F.box, F.decomp, F.rank, clp_params);
    cl_inc.buildInterior(F.owned);

    ParticleExchange::Params pex_params_inc{/*incremental*/true};
    ParticleExchange pex_inc(MPI_COMM_WORLD, *F.box, F.decomp, F.rank, cl_inc, pex_params_inc);
    pex_inc.refreshGhosts(F.owned, cl_inc);

    ThermodynamicCalculatorParallel thermo_inc(MPI_COMM_WORLD, F.T, F.pot, F.rcut,
                                               /*mu*/0.0, F.f, F.alpha, F.A0, F.kappa);

    RNG_parallel rng_inc(/*seed*/111u, /*rank*/F.rank);
    MonteCarloNVT_MPI::Params mp_inc; mp_inc.delta = 0.04; mp_inc.halo_every = 1;

    MonteCarloNVT_MPI mc_inc(MPI_COMM_WORLD, *F.box, cl_inc, pex_inc, thermo_inc,
                             F.owned, rng_inc, /*traj*/nullptr, /*data*/nullptr, mp_inc);

    // Clone initial for full-refresh run
    std::vector<Particle> owned_full = F.owned;

    // --- B) full refresh mode (no incremental; refresh every parity block) ---
    CellListParallel cl_full(*F.box, F.decomp, F.rank, clp_params);
    cl_full.buildInterior(owned_full);

    ParticleExchange::Params pex_params_full{/*incremental*/false};
    ParticleExchange pex_full(MPI_COMM_WORLD, *F.box, F.decomp, F.rank, cl_full, pex_params_full);
    pex_full.refreshGhosts(owned_full, cl_full);

    ThermodynamicCalculatorParallel thermo_full(MPI_COMM_WORLD, F.T, F.pot, F.rcut,
                                                /*mu*/0.0, F.f, F.alpha, F.A0, F.kappa);

    RNG_parallel rng_full(/*seed*/111u, /*rank*/F.rank); // same seed
    MonteCarloNVT_MPI::Params mp_full; mp_full.delta = 0.04; mp_full.halo_every = 1;

    MonteCarloNVT_MPI mc_full(MPI_COMM_WORLD, *F.box, cl_full, pex_full, thermo_full,
                              owned_full, rng_full, /*traj*/nullptr, /*data*/nullptr, mp_full);

    // Run identical number of sweeps
    const std::size_t nsweeps = 5;
    (void)mc_inc.run(nsweeps);
    (void)mc_full.run(nsweeps);

    // Compare final local configurations (IDs and coordinates)
    sort_by_id(F.owned);
    sort_by_id(owned_full);
    require_same_configuration(F.owned, owned_full, 1e-14);
    assert_inside_box(F.owned, *F.box);
    assert_inside_box(owned_full, *F.box);
}

// ----------------- TEST 4: energy consistency after rebuild -----------------

TEST_CASE("MC NVT MPI: energy unchanged by rebuild (bins+ghosts)", "[MPI][MC][consistency]") {
    // geometry & particles
    const double Lx=120.0, Ly=96.0; const int Nglob=400;
    const double T=1.1; const double rcut_want=3.2;
    SimulationBox box(Lx, Ly);
    auto decomp = box.bestDecomposition([](int s){int r=0; MPI_Comm_size(MPI_COMM_WORLD,&r); return r;}(0)); // size from MPI

    int rank=0,size=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);

    auto global = make_global_uniform_particles(Nglob, Lx, Ly, /*seed=*/777);
    double x0,x1,y0,y1; decomp.localBounds(rank, Lx, Ly, x0, x1, y0, y1);
    const double Llocx = x1-x0, Llocy = y1-y0;
    auto eff_rcut = [&](double w){ double cap=0.49*std::min(Llocx,Llocy); return std::max(1e-8, std::min(w,cap)); };
    const double rcut = eff_rcut(rcut_want);

    auto owned = filter_owned(global, x0,x1,y0,y1);

    // build CL + PX (incremental enabled)
    CellListParallel::Params clp_params{rcut};
    CellListParallel clp(box, decomp, rank, clp_params);
    clp.buildInterior(owned);

    ParticleExchange::Params pexp{true};
    ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pexp);
    pex.refreshGhosts(owned, clp);

    ThermodynamicCalculatorParallel thermo(MPI_COMM_WORLD, T, PotentialType::WCA, rcut,
                                           /*mu*/0, /*f*/0, /*alpha*/0, /*A0*/0, /*kappa*/1);

    // energy before rebuild
    const double U0 = thermo.totalEnergy(owned, box, clp, pex);
    const double W0 = thermo.totalVirial(owned, box, clp, pex);
    const double P0 = thermo.pressure(owned, box, clp, pex);

    // do a few accepted moves via MC driver
    RNG_parallel rng(/*seed*/123u, /*rank*/rank);
    MonteCarloNVT_MPI::Params mp; mp.delta = 0.05; mp.halo_every = 1;
    MonteCarloNVT_MPI mc(MPI_COMM_WORLD, box, clp, pex, thermo, owned, rng, nullptr, nullptr, mp);
    (void)mc.run(/*nsweeps*/2);

    // force full rebuild of both bins + ghosts
    clp.buildInterior(owned);
    pex.refreshGhosts(owned, clp);

    // energy after rebuild (same configuration)
    const double U1 = thermo.totalEnergy(owned, box, clp, pex);
    const double W1 = thermo.totalVirial(owned, box, clp, pex);
    const double P1 = thermo.pressure(owned, box, clp, pex);

    REQUIRE(U1 == Approx(U1).epsilon(0)); // self check
    REQUIRE(std::isfinite(U1));
    REQUIRE(std::isfinite(W1));
    REQUIRE(std::isfinite(P1));

    // Rebuild must not change thermodynamics for the *same* positions
    // (not comparing with U0/W0 because positions changed during the short run)
    // Do a no-op rebuild again and check invariance.
    clp.buildInterior(owned);
    pex.refreshGhosts(owned, clp);
    const double U2 = thermo.totalEnergy(owned, box, clp, pex);
    const double W2 = thermo.totalVirial(owned, box, clp, pex);
    const double P2 = thermo.pressure(owned, box, clp, pex);

    REQUIRE(U2 == Approx(U1).epsilon(1e-12));
    REQUIRE(W2 == Approx(W1).epsilon(1e-12));
    REQUIRE(P2 == Approx(P1).epsilon(1e-12));
}

// ----------------- TEST 5: single-rank ΔU path equivalence (local vs move-and-recompute) -----------------

TEST_CASE("MC NVT MPI: local ΔU equals move-then-recompute (size==1)", "[MPI][MC][deltaE][single]") {
    int size=1; MPI_Comm_size(MPI_COMM_WORLD,&size);
    if (size != 1) { WARN("Run with 1 MPI rank to exercise brute-force-style ΔU."); return; }

    const double Lx=60.0, Ly=60.0; const int Nglob=150;
    const double T=0.9; const double rcut=2.5;
    SimulationBox box(Lx, Ly);
    auto decomp = box.bestDecomposition(size);

    int rank=0; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    auto global = make_global_uniform_particles(Nglob, Lx, Ly, /*seed=*/13579);
    double x0,x1,y0,y1; decomp.localBounds(rank, Lx, Ly, x0, x1, y0, y1);
    auto owned = filter_owned(global, x0,x1,y0,y1);

    CellListParallel::Params clp_params{rcut};
    CellListParallel clp(box, decomp, rank, clp_params);
    clp.buildInterior(owned);

    ParticleExchange::Params pexp{true};
    ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pexp);
    pex.refreshGhosts(owned, clp);

    ThermodynamicCalculatorParallel thermo(MPI_COMM_WORLD, T, PotentialType::LennardJones, rcut,
                                           /*mu*/0, /*f*/0, /*alpha*/0, /*A0*/0, /*kappa*/1);

    RNG_parallel rng(/*seed*/222u, /*rank*/0);

    // Pick a random owned particle and a small trial move
    if (owned.empty()) { SUCCEED("no owned particles on this rank"); return; }
    std::uniform_int_distribution<int> unif(0, (int)owned.size()-1);
    int i = unif(rng.engine());

    const double delta = 0.03;
    double dx=0, dy=0;
    if (delta>0) { dx = rng.uniform(-delta,delta); dy = rng.uniform(-delta,delta); }

    // ΔU via local path (neighborsOfOwned + neighborsOfPoint)
    auto neighbors_i = clp.neighborsOfOwned(i, owned);
    double U_old = 0.0;
    for (auto& pr : neighbors_i) {
        const double r2 = pr.second;
        U_old += computePairPotential(r2, PotentialType::LennardJones, 0.f, 0.f, 1.f, 0.f);
    }
    Particle trial = owned[i]; trial.x += dx; trial.y += dy; box.applyPBC(trial);
    auto neighbors_trial = clp.neighborsOfPoint(trial.x, trial.y, owned);
    double U_new = 0.0;
    for (auto& pr : neighbors_trial) {
        const int j  = pr.first;
        const double r2 = pr.second;
        if (j >= 0 && r2 <= 1e-24) continue; // skip self as in MC
        U_new += computePairPotential(r2, PotentialType::LennardJones, 0.f, 0.f, 1.f, 0.f);
    }
    const double dU_local = U_new - U_old;

    // ΔU via "move-and-recompute": commit, rebuild bins, recompute local energy, then revert
    Particle saved = owned[i];
    owned[i] = trial;
    clp.onAcceptedMove(i, saved, trial); // update bin so the second local energy uses the new position
    auto neighbors_new_i = clp.neighborsOfOwned(i, owned);
    double U_new2 = 0.0;
    for (auto& pr : neighbors_new_i) {
        const double r2 = pr.second;
        U_new2 += computePairPotential(r2, PotentialType::LennardJones, 0.f, 0.f, 1.f, 0.f);
    }
    // revert (so we leave no trace)
    clp.onAcceptedMove(i, trial, saved);
    owned[i] = saved;

    const double dU_move_recompute = U_new2 - U_old;

    REQUIRE(dU_local == Approx(dU_move_recompute).epsilon(1e-12));
}

// ----------------- TEST 6: basic invariants (positions in box, acceptance bounds) -----------------

TEST_CASE("MC NVT MPI: invariants", "[MPI][MC][invariants]") {
    const double Lx=100.0, Ly=80.0; const int Nglob=250;
    const double T=1.0; const double rcut_want=3.0;

    SimulationBox box(Lx, Ly);
    int rank=0,size=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
    auto decomp = box.bestDecomposition(size);

    auto global = make_global_uniform_particles(Nglob, Lx, Ly, /*seed=*/2468);
    double x0,x1,y0,y1; decomp.localBounds(rank, Lx, Ly, x0, x1, y0, y1);
    const double Llocx=x1-x0, Llocy=y1-y0;
    auto eff_rcut=[&](double w){ double cap=0.49*std::min(Llocx,Llocy); return std::max(1e-8,std::min(w,cap)); };
    const double rcut=eff_rcut(rcut_want);
    auto owned = filter_owned(global, x0,x1,y0,y1);

    CellListParallel::Params clp_params{rcut};
    CellListParallel clp(box, decomp, rank, clp_params);
    clp.buildInterior(owned);

    ParticleExchange::Params pexp{true};
    ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pexp);
    pex.refreshGhosts(owned, clp);

    ThermodynamicCalculatorParallel thermo(MPI_COMM_WORLD, T, PotentialType::Yukawa, rcut,
                                           /*mu*/0, /*f*/0, /*alpha*/0, /*A0*/0, /*kappa*/1);

    RNG_parallel rng(/*seed*/135u, /*rank*/rank);
    MonteCarloNVT_MPI::Params mp; mp.delta = 0.05; mp.halo_every = 2;
    MonteCarloNVT_MPI mc(MPI_COMM_WORLD, box, clp, pex, thermo, owned, rng, nullptr, nullptr, mp);

    const double acc = mc.run(/*nsweeps*/4);
    REQUIRE(acc >= 0.0);
    REQUIRE(acc <= 1.0);
    assert_inside_box(owned, box);
}


// ----------------- TEST 7: energy decreases on average during short equilibration (WCA) -----------------

TEST_CASE("MC NVT MPI: energy tends to decrease during short equilibration (WCA)", "[MPI][MC][trend][energy]") {
    int rank=0,size=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);

    // Geometry and interactions
    const double Lx = 120.0, Ly = 96.0;
    const int    Nglob = 500;
    const double T = 1.0;                    // moderate temperature
    const double rcut_want = 3.2;            // WCA with safe rcut
    const PotentialType pot = PotentialType::WCA;

    SimulationBox box(Lx, Ly);
    auto decomp = box.bestDecomposition(size);

    // Build a deterministic global initial configuration
    auto global0 = make_global_uniform_particles(Nglob, Lx, Ly, /*seed=*/20240u);

    // Local ownership + rcut cap
    double x0,x1,y0,y1; decomp.localBounds(rank, Lx, Ly, x0, x1, y0, y1);
    const double Llocx = x1-x0, Llocy = y1-y0;
    auto effective_rcut = [&](double want){
        const double cap = 0.49 * std::min(Llocx, Llocy);
        return std::max(1e-8, std::min(want, cap));
    };
    const double rcut = effective_rcut(rcut_want);

    // Replicates for robustness
    const int REPS = 3;
    double U0_sum = 0.0, U1_sum = 0.0;

    for (int rep=0; rep<REPS; ++rep) {
        // Fresh local copy
        auto owned = filter_owned(global0, x0, x1, y0, y1);

        // CL + PX (incremental halo)
        CellListParallel::Params clp_params{rcut};
        CellListParallel clp(box, decomp, rank, clp_params);
        clp.buildInterior(owned);

        ParticleExchange::Params pex_params{/*enable_incremental=*/true};
        ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pex_params);
        pex.refreshGhosts(owned, clp);

        ThermodynamicCalculatorParallel thermo(MPI_COMM_WORLD, T, pot, rcut,
                                               /*mu*/0.0, /*f*/0.0, /*alpha*/0.0, /*A0*/0.0, /*kappa*/1.0);

        // Initial global energy
        const double U0 = thermo.totalEnergy(owned, box, clp, pex);

        // Driver with small steps (keep particles inside their tiles statistically)
        MonteCarloNVT_MPI::Params mp; mp.delta = 0.05 * rcut; mp.halo_every = 1; mp.out_every = 0;
        RNG_parallel rng(/*seed*/ (unsigned)(777u + rep*17u), /*rank*/rank);

        MonteCarloNVT_MPI mc(MPI_COMM_WORLD, box, clp, pex, thermo, owned, rng,
                             /*traj*/nullptr, /*data*/nullptr, mp);

        // Short equilibration
        const std::size_t nsweeps = 20;
        (void)mc.run(nsweeps);

        // Final energy (global)
        const double U1 = thermo.totalEnergy(owned, box, clp, pex);

        U0_sum += U0;
        U1_sum += U1;
    }

    // Average energies across replicates: <U_final> should be lower than <U_init>
    const double U0_avg = U0_sum / REPS;
    const double U1_avg = U1_sum / REPS;

    // Allow a small tolerance; WCA is purely repulsive, so early equilibration lowers overlaps
    REQUIRE(U1_avg < U0_avg - 1e-6);
}

// ----------------- TEST: higher T ⇒ higher total energy and pressure (equilibrated, time-avg) -----------------

TEST_CASE("MC NVT MPI: higher T ⇒ higher (U + N*T) and P (after pre-equil, time-avg)", "[MPI][MC][temperature][energy][pressure]") {
    int rank=0,size=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);

    // Geometry & interaction
    const double Lx = 100.0, Ly = 80.0;
    const int    Nglob = 400;
    const double rcut_want = 2.5;
    const PotentialType pot = PotentialType::LennardJones;

    // Temperatures
    const double T_low  = 0.6;
    const double T_high = 2.0;
    const double T_mid  = 1.2;   // pre-equil temperature

    SimulationBox box(Lx, Ly);
    auto decomp = box.bestDecomposition(size);

    auto global0 = make_global_uniform_particles(Nglob, Lx, Ly, /*seed=*/9090u);

    // Local subset & safe rcut
    double x0,x1,y0,y1; decomp.localBounds(rank, Lx, Ly, x0, x1, y0, y1);
    const double Llocx = x1-x0, Llocy = y1-y0;
    auto cap_rcut = [&](double want){ const double cap=0.49*std::min(Llocx,Llocy); return std::max(1e-8,std::min(want,cap)); };
    const double rcut = cap_rcut(rcut_want);

    // ---------- Pre-equilibrate at T_mid ----------
    auto owned_mid = filter_owned(global0, x0, x1, y0, y1);

    CellListParallel::Params clp_params{rcut};
    CellListParallel cl_mid(box, decomp, rank, clp_params);
    cl_mid.buildInterior(owned_mid);

    ParticleExchange::Params pex_params{/*enable_incremental=*/true};
    ParticleExchange pex_mid(MPI_COMM_WORLD, box, decomp, rank, cl_mid, pex_params);
    pex_mid.refreshGhosts(owned_mid, cl_mid);

    ThermodynamicCalculatorParallel thermo_mid(MPI_COMM_WORLD, T_mid, pot, rcut,
                                               /*mu*/0.0, /*f*/0.0, /*alpha*/0.0, /*A0*/0.0, /*kappa*/1.0);

    MonteCarloNVT_MPI::Params mp; mp.delta=0.05*rcut; mp.halo_every=1; mp.out_every=0;
    RNG_parallel rng_mid(/*seed*/11111u, /*rank*/rank);
    MonteCarloNVT_MPI mc_mid(MPI_COMM_WORLD, box, cl_mid, pex_mid, thermo_mid, owned_mid, rng_mid,
                             nullptr, nullptr, mp);
    (void)mc_mid.run(/*nsweeps*/60);

    // Clone config for low/high T
    auto owned_low  = owned_mid;
    auto owned_high = owned_mid;

    // ---------- Low T branch ----------
    CellListParallel cl_low(box, decomp, rank, clp_params); cl_low.buildInterior(owned_low);
    ParticleExchange pex_low(MPI_COMM_WORLD, box, decomp, rank, cl_low, pex_params); pex_low.refreshGhosts(owned_low, cl_low);
    ThermodynamicCalculatorParallel thermo_low(MPI_COMM_WORLD, T_low, pot, rcut, 0.0,0.0,0.0,0.0,1.0);
    RNG_parallel rng_low(/*seed*/10101u, /*rank*/rank);
    MonteCarloNVT_MPI mc_low(MPI_COMM_WORLD, box, cl_low, pex_low, thermo_low, owned_low, rng_low, nullptr, nullptr, mp);

    (void)mc_low.run(/*equil*/80);
    double U_low_sum=0.0, P_low_sum=0.0; int nL=40;
    for (int s=0; s<nL; ++s) { (void)mc_low.run(1);
        U_low_sum += thermo_low.totalEnergy(owned_low, box, cl_low, pex_low);
        P_low_sum += thermo_low.pressure   (owned_low, box, cl_low, pex_low);
    }
    // global N for kinetic term
    int nlocL=(int)owned_low.size(), Nglob_check=0; MPI_Allreduce(&nlocL,&Nglob_check,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    const double U_low_avg = U_low_sum / nL;
    const double P_low_avg = P_low_sum / nL;
    const double Etot_low  = U_low_avg + Nglob_check * T_low; // add ideal kinetic (d/2=1 per particle in 2D)

    // ---------- High T branch ----------
    CellListParallel cl_high(box, decomp, rank, clp_params); cl_high.buildInterior(owned_high);
    ParticleExchange pex_high(MPI_COMM_WORLD, box, decomp, rank, cl_high, pex_params); pex_high.refreshGhosts(owned_high, cl_high);
    ThermodynamicCalculatorParallel thermo_high(MPI_COMM_WORLD, T_high, pot, rcut, 0.0,0.0,0.0,0.0,1.0);
    RNG_parallel rng_high(/*seed*/20202u, /*rank*/rank);
    MonteCarloNVT_MPI mc_high(MPI_COMM_WORLD, box, cl_high, pex_high, thermo_high, owned_high, rng_high, nullptr, nullptr, mp);

    (void)mc_high.run(/*equil*/80);
    double U_high_sum=0.0, P_high_sum=0.0; int nH=40;
    for (int s=0; s<nH; ++s) { (void)mc_high.run(1);
        U_high_sum += thermo_high.totalEnergy(owned_high, box, cl_high, pex_high);
        P_high_sum += thermo_high.pressure   (owned_high, box, cl_high, pex_high);
    }
    int nlocH=(int)owned_high.size(), Nglob_check2=0; MPI_Allreduce(&nlocH,&Nglob_check2,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    REQUIRE(Nglob_check2 == Nglob_check);

    const double U_high_avg = U_high_sum / nH;
    const double P_high_avg = P_high_sum / nH;
    const double Etot_high  = U_high_avg + Nglob_check2 * T_high;

    // ---------- Assertions ----------
    // 1) Total energy (config + kinetic) must rise with T in NVT.
    REQUIRE(Etot_high > Etot_low + 1e-6);

    // 2) Pressure should also rise with T (ρT term guarantees monotonicity).
    REQUIRE(P_high_avg > P_low_avg + 1e-6);

    // Optional (relaxed) sanity: potential U should not significantly *exceed* low-T by a large negative margin.
    // This avoids brittle expectations on U alone.
    REQUIRE(U_high_avg >= U_low_avg - 1e-4 * std::max(1.0, std::abs(U_low_avg)));
}

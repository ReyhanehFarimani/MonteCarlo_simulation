// file: test_thermodynamic_calculator_parallel.cpp
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
#include "../mpi_src/potential.h"
#include "../mpi_src/thermodynamic_calculator_parallel.h"

// ----------------- helpers -----------------

static inline int make_id_linear(int gidx) {
    REQUIRE(gidx >= 0);
    REQUIRE(gidx <= std::numeric_limits<int>::max());
    return gidx;
}

static inline int make_id_grid(int i, int nx, int j) {
    long long v = 1LL * (i + j*nx);
    REQUIRE(v <= std::numeric_limits<int>::max());
    return static_cast<int>(v);
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

static double serial_total_energy(const std::vector<Particle>& global,
                                  const SimulationBox& box,
                                  PotentialType pot,
                                  double rcut,
                                  double f_prime, double f_d_prime,
                                  double kappa, double alpha)
{
    const double r2cut = rcut * rcut;
    double U = 0.0;
    const int N = (int)global.size();
    for (int i=0;i<N;++i) {
        for (int j=i+1;j<N;++j) {
            double r2 = box.minimumImageDistanceSquared(global[i], global[j]);
            if (r2 <= r2cut) {
                U += computePairPotential(r2, pot,
                                          (float)f_prime, (float)f_d_prime,
                                          (float)kappa, (float)alpha);
            }
        }
    }
    return U;
}

static double serial_total_virial(const std::vector<Particle>& global,
                                  const SimulationBox& box,
                                  PotentialType pot,
                                  double rcut,
                                  double f_prime, double f_d_prime,
                                  double kappa, double alpha)
{
    const double r2cut = rcut * rcut;
    double W = 0.0;
    const int N = (int)global.size();
    for (int i=0;i<N;++i) {
        for (int j=i+1;j<N;++j) {
            double r2 = box.minimumImageDistanceSquared(global[i], global[j]);
            if (r2 <= r2cut) {
                W += computePairForce(r2, pot,
                                      (float)f_prime, (float)f_d_prime,
                                      (float)kappa, (float)alpha);
            }
        }
    }
    return W;
}

static void translate_all(std::vector<Particle>& global,
                          double dx, double dy,
                          const SimulationBox& box)
{
    for (auto& p : global) { p.updatePosition(dx, dy); box.applyPBC(p); }
}

// ----------------- TEST 1: WCA -----------------

TEST_CASE("Thermodynamics: parallel vs serial (WCA, batch halo)", "[MPI][thermo][parallel][WCA]") {
    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const double Lx = 120.0, Ly = 96.0;
    const double rcut_want = 3.5;
    const double T = 1.25;
    const PotentialType pot = PotentialType::WCA;
    const double f = 0.0, alpha = 0.0, A0 = 0.0, kappa = 1.0;
    const int Nglob = 400;

    SimulationBox box(Lx, Ly);
    auto decomp = box.bestDecomposition(size);

    auto global = make_global_uniform_particles(Nglob, Lx, Ly, /*seed=*/777);

    double x0,x1,y0,y1;
    decomp.localBounds(rank, Lx, Ly, x0, x1, y0, y1);
    const double Llocx = x1 - x0, Llocy = y1 - y0;
    auto effective_rcut = [&](double want){
        const double cap = 0.49 * std::min(Llocx, Llocy);
        return std::max(1e-8, std::min(want, cap));
    };
    const double rcut = effective_rcut(rcut_want);

    auto owned = filter_owned(global, x0,x1,y0,y1);

    CellListParallel::Params clp_params{rcut};
    CellListParallel clp(box, decomp, rank, clp_params);
    clp.buildInterior(owned);

    ParticleExchange::Params pex_params{/*enable_incremental=*/false};
    ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pex_params);
    pex.refreshGhosts(owned, clp);

    ThermodynamicCalculatorParallel tpar(MPI_COMM_WORLD, T, pot, rcut,
                                         /*mu*/0.0, f, alpha, A0, kappa);

    const double Upar = tpar.totalEnergy(owned, box, clp, pex);
    const double Wpar = tpar.totalVirial(owned, box, clp, pex);
    const double Ppar = tpar.pressure(owned, box, clp, pex);

    double User=0.0, Wser=0.0, Pser=0.0;
    if (rank == 0) {
        const double f_prime = (2.0 + 9.0*f*f)/24.0;
        const double f_d_prime = (kappa==0.0?0.0:(A0 * f * f / kappa));
        User = serial_total_energy(global, box, pot, rcut, f_prime, f_d_prime, kappa, alpha);
        Wser = serial_total_virial(global, box, pot, rcut, f_prime, f_d_prime, kappa, alpha);
        const double rho = (double)Nglob / box.getV();
        Pser = rho * T + Wser / (2.0 * box.getV());
    }
    MPI_Bcast(&User, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Wser, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Pser, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    REQUIRE(Upar == Approx(User).epsilon(1e-12));
    REQUIRE(Wpar == Approx(Wser).epsilon(1e-12));
    REQUIRE(Ppar == Approx(Pser).epsilon(1e-12));

    const double dx = 0.13 * rcut;
    const double dy = -0.07 * rcut;

    for (auto& p : owned) { p.updatePosition(dx, dy); box.applyPBC(p); }
    clp.buildInterior(owned);
    int sent = pex.migrate(owned, clp);
    REQUIRE(sent >= 0);
    pex.refreshGhosts(owned, clp);

    const double Upar2 = tpar.totalEnergy(owned, box, clp, pex);
    const double Wpar2 = tpar.totalVirial(owned, box, clp, pex);
    const double Ppar2 = tpar.pressure(owned, box, clp, pex);

    if (rank == 0) {
        translate_all(global, dx, dy, box);
        const double f_prime = (2.0 + 9.0*f*f)/24.0;
        const double f_d_prime = (kappa==0.0?0.0:(A0 * f * f / kappa));
        User = serial_total_energy(global, box, pot, rcut, f_prime, f_d_prime, kappa, alpha);
        Wser = serial_total_virial(global, box, pot, rcut, f_prime, f_d_prime, kappa, alpha);
        const double rho = (double)Nglob / box.getV();
        Pser = rho * T + Wser / (2.0 * box.getV());
    }
    MPI_Bcast(&User, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Wser, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Pser, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    REQUIRE(Upar2 == Approx(User).epsilon(1e-12));
    REQUIRE(Wpar2 == Approx(Wser).epsilon(1e-12));
    REQUIRE(Ppar2 == Approx(Pser).epsilon(1e-12));
}

// ----------------- TEST 2: LJ random -----------------

TEST_CASE("Thermodynamics: parallel vs serial (LJ, batch halo)", "[MPI][thermo][parallel][LJ]") {
    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const double Lx = 100.0, Ly = 90.0;
    const double rcut_want = 2.5;
    const double T = 0.9;
    const PotentialType pot = PotentialType::LennardJones;
    const double f = 0.0, alpha = 0.0, A0 = 0.0, kappa = 1.0;
    const int Nglob = 300;

    SimulationBox box(Lx, Ly);
    auto decomp = box.bestDecomposition(size);

    auto global = make_global_uniform_particles(Nglob, Lx, Ly, /*seed=*/2024);

    double x0,x1,y0,y1;
    decomp.localBounds(rank, Lx, Ly, x0, x1, y0, y1);
    const double Llocx = x1 - x0, Llocy = y1 - y0;
    auto effective_rcut = [&](double want){
        const double cap = 0.49 * std::min(Llocx, Llocy);
        return std::max(1e-8, std::min(want, cap));
    };
    const double rcut = effective_rcut(rcut_want);

    auto owned = filter_owned(global, x0,x1,y0,y1);

    CellListParallel::Params clp_params{rcut};
    CellListParallel clp(box, decomp, rank, clp_params);
    clp.buildInterior(owned);

    ParticleExchange::Params pex_params{/*enable_incremental=*/false};
    ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pex_params);
    pex.refreshGhosts(owned, clp);

    ThermodynamicCalculatorParallel tpar(MPI_COMM_WORLD, T, pot, rcut,
                                         /*mu*/0.0, f, alpha, A0, kappa);

    const double Upar = tpar.totalEnergy(owned, box, clp, pex);
    const double Wpar = tpar.totalVirial(owned, box, clp, pex);
    const double Ppar = tpar.pressure(owned, box, clp, pex);

    double User=0.0, Wser=0.0, Pser=0.0;
    if (rank == 0) {
        const double f_prime = (2.0 + 9.0*f*f)/24.0;
        const double f_d_prime = (kappa==0.0?0.0:(A0 * f * f / kappa));
        User = serial_total_energy(global, box, pot, rcut, f_prime, f_d_prime, kappa, alpha);
        Wser = serial_total_virial(global, box, pot, rcut, f_prime, f_d_prime, kappa, alpha);
        const double rho = (double)Nglob / box.getV();
        Pser = rho * T + Wser / (2.0 * box.getV());
    }
    MPI_Bcast(&User, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Wser, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Pser, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    REQUIRE(Upar == Approx(User).epsilon(1e-12));
    REQUIRE(Wpar == Approx(Wser).epsilon(1e-12));
    REQUIRE(Ppar == Approx(Pser).epsilon(1e-12));
}

// ----------------- TEST 3: LJ lattice (analytic) -----------------

static std::vector<Particle>
make_square_lattice(int nx, int ny, double a)
{
    std::vector<Particle> v;
    v.reserve(nx*ny);
    for (int j=0; j<ny; ++j)
        for (int i=0; i<nx; ++i)
            v.push_back(Particle{(i + 0.5) * a, (j + 0.5) * a, make_id_grid(i, nx, j)});
    return v;
}

static std::pair<double,double>
lj_square_lattice_per_particle(double a, double rcut)
{
    const double rc2 = rcut*rcut;
    const int mmax = static_cast<int>(std::floor(rcut / a));
    double U = 0.0, W = 0.0;

    for (int m = -mmax; m <= mmax; ++m) {
        for (int n = -mmax; n <= mmax; ++n) {
            if (m==0 && n==0) continue;
            const double r2 = a*a * (m*m + n*n);
            if (r2 > rc2) continue;
            const double inv_r2 = 1.0 / r2;
            const double r6inv  = inv_r2 * inv_r2 * inv_r2;
            U += 4.0 * r6inv * (r6inv - 1.0);
            W += 48.0 * r6inv * r6inv - 24.0 * r6inv;
        }
    }
    return {U, W};
}

TEST_CASE("LJ crystal: parallel vs lattice-sum (square lattice)", "[MPI][thermo][LJ][lattice]") {
    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int nx = 8, ny = 8;
    const double a = std::pow(2.0, 1.0/6.0);
    const double Lx = nx * a;
    const double Ly = ny * a;

    const double rcut_want = 2.5;
    const double T = 0.9;
    const PotentialType pot = PotentialType::LennardJones;
    const double f=0.0, alpha=0.0, A0=0.0, kappa=1.0;

    SimulationBox box(Lx, Ly);
    auto decomp = box.bestDecomposition(size);

    auto global = make_square_lattice(nx, ny, a);

    double x0,x1,y0,y1;
    decomp.localBounds(rank, Lx, Ly, x0, x1, y0, y1);
    const double Llocx = x1 - x0, Llocy = y1 - y0;
    auto effective_rcut = [&](double want){
        const double cap = 0.49 * std::min(Llocx, Llocy);
        return std::max(1e-8, std::min(want, cap));
    };
    const double rcut = effective_rcut(rcut_want);
    REQUIRE(rcut < 0.5 * std::min(Lx, Ly));

    auto owned = filter_owned(global, x0, x1, y0, y1);

    CellListParallel::Params clp_params{rcut};
    CellListParallel clp(box, decomp, rank, clp_params);
    clp.buildInterior(owned);

    ParticleExchange::Params pex_params{/*enable_incremental=*/false};
    ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pex_params);
    pex.refreshGhosts(owned, clp);

    ThermodynamicCalculatorParallel tpar(MPI_COMM_WORLD, T, pot, rcut,
                                         /*mu*/0.0, f, alpha, A0, kappa);

    const double Upar = tpar.totalEnergy(owned, box, clp, pex);
    const double Wpar = tpar.totalVirial(owned, box, clp, pex);
    const double Ppar = tpar.pressure(owned, box, clp, pex);

    double Ulat=0.0, Wlat=0.0, Plat=0.0;
    if (rank == 0) {
        auto [Upp, Wpp] = lj_square_lattice_per_particle(a, rcut);
        const int N = nx * ny;
        Ulat = 0.5 * N * Upp;
        Wlat = 0.5 * N * Wpp;
        const double rho = N / (Lx * Ly);
        Plat = rho * T + Wlat / (2.0 * (Lx * Ly));
    }
    MPI_Bcast(&Ulat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Wlat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Plat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    REQUIRE(Upar == Approx(Ulat).epsilon(1e-12));
    REQUIRE(Wpar == Approx(Wlat).epsilon(1e-12));
    REQUIRE(Ppar == Approx(Plat).epsilon(1e-12));

    for (auto& p : owned) { p.updatePosition(0.11*a, -0.07*a); box.applyPBC(p); }
    clp.buildInterior(owned);
    int sent = pex.migrate(owned, clp);
    REQUIRE(sent >= 0);
    pex.refreshGhosts(owned, clp);

    const double Upar2 = tpar.totalEnergy(owned, box, clp, pex);
    const double Wpar2 = tpar.totalVirial(owned, box, clp, pex);
    const double Ppar2 = tpar.pressure(owned, box, clp, pex);

    REQUIRE(Upar2 == Approx(Ulat).epsilon(1e-12));
    REQUIRE(Wpar2 == Approx(Wlat).epsilon(1e-12));
    REQUIRE(Ppar2 == Approx(Plat).epsilon(1e-12));
}


// --- test: exact-cutoff inclusion/exclusion across a rank boundary (LJ) ---
TEST_CASE("Rank-boundary pair at and above cutoff (LJ)", "[MPI][neighbors][cutoff][boundary]") {
    int rank=0, size=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
    if (size < 2) { WARN("Needs size >= 2; skipping."); return; }

    // Geometry and decomposition
    const double Lx = 80.0, Ly = 40.0;
    SimulationBox box(Lx, Ly);
    auto decomp = box.bestDecomposition(size);

    // Identify rank 0 and its East neighbor
    auto coords0 = decomp.coordsOf(0);
    const int rx0 = coords0.first, ry0 = coords0.second;
    auto wrap = [](int i, int n){ i%=n; if (i<0) i+=n; return i; };
    const int Px = decomp.Px, Py = decomp.Py;
    auto id = [&](int ix, int iy){ ix=wrap(ix,Px); iy=wrap(iy,Py); return iy*Px + ix; };
    const int rankE = id(rx0+1, ry0);

    // Local bounds and per-rank effective rcut
    double x0,x1,y0,y1;            decomp.localBounds(rank, Lx, Ly, x0, x1, y0, y1);
    double x0_0,x1_0,y0_0,y1_0;    decomp.localBounds(0,    Lx, Ly, x0_0, x1_0, y0_0, y1_0);
    double x0_E,x1_E,y0_E,y1_E;    decomp.localBounds(rankE,Lx, Ly, x0_E, x1_E, y0_E, y1_E);

    const double Llocx0 = x1_0 - x0_0, Llocy0 = y1_0 - y0_0;
    auto effective_rcut = [&](double want, double Llocx, double Llocy){
        const double cap = 0.49 * std::min(Llocx, Llocy);   // ensure Lloc >= 2*rcut
        return std::max(1e-8, std::min(want, cap));
    };

    // Desired rcut; ensure it's valid for BOTH tiles so both ranks can build
    const double rcut_desired = 3.0;
    const double rcut0 = effective_rcut(rcut_desired, Llocx0, Llocy0);

    // We'll build clp on each rank with its *own* cap; for symmetry, also cap on E
    const double LlocxE = x1_E - x0_E, LlocyE = y1_E - y0_E;
    const double rcutE = effective_rcut(rcut_desired, LlocxE, LlocyE);

    // Use the MIN to define the test target distance so both sides are happy
    const double rcut = std::min(rcut0, rcutE);
    REQUIRE(rcut > 0.0);

    // Helper to run one scenario (distance = rcut or rcut + eps)
    auto run_case = [&](double extra)->void {
        // Build owned per rank: rank 0 has A, rankE has B
        std::vector<Particle> owned;

        // Place both at same y (mid), straddling the vertical face x = x1_0 == x0_E
        const double ymid = 0.5*(y0_0 + y1_0);
        const double eps_inside = 0.25 * rcut; // keep well inside border but within halo
        // Particle A on rank 0 (just inside its right face)
        Particle A{x1_0 - eps_inside, ymid, 1001};
        // Particle B on rankE (just inside its left face) so that separation is rcut + extra
        Particle B{x0_E + (rcut + extra) - eps_inside, ymid, 2002};

        if (rank == 0)      owned.push_back(A);
        if (rank == rankE)  owned.push_back(B);

        // Build cell lists with per-rank rcut caps
        CellListParallel::Params clp_params{ (rank==rankE)? rcutE : rcut0 };
        CellListParallel clp(box, decomp, rank, clp_params);
        clp.buildInterior(owned);

        ParticleExchange::Params pex_params{false};
        ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pex_params);
        pex.refreshGhosts(owned, clp);

        // Check neighbors on rank 0 for its owned A (if present)
        bool found_inclusion = false;
        double found_r2 = -1.0;
        if (rank == 0 && !owned.empty()) {
            auto nb = clp.neighborsOfOwned(0, owned);
            for (auto& pr : nb) {
                if (pr.first < 0) { // ghost
                    found_inclusion = true;
                    found_r2 = pr.second;
                    break;
                }
            }
        }
        // Reduce boolean to root for clarity
        int loc_has = (rank==0 && !owned.empty() && found_inclusion) ? 1 : 0;
        int any_has = 0;
        MPI_Allreduce(&loc_has, &any_has, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        // Parallel thermo (LJ) vs serial reference on the two-particle global config
        ThermodynamicCalculatorParallel tpar(MPI_COMM_WORLD, /*T*/1.0, PotentialType::LennardJones,
                                             /*rcut*/rcut, /*mu*/0, /*f*/0, /*alpha*/0, /*A0*/0, /*kappa*/1);
        const double Upar = tpar.totalEnergy(owned, box, clp, pex);
        const double Wpar = tpar.totalVirial(owned, box, clp, pex);
        const double Ppar = tpar.pressure(owned, box, clp, pex);

        // Build the *global* two-particle set (same coordinates) for serial reference on rank 0
        double Uref=0.0, Wref=0.0, Pref=0.0;
        if (rank == 0) {
            std::vector<Particle> global2 = {A, B};
            const double fprime=0.0, fdprime=0.0, kappa=1.0, alpha=0.0;
            Uref = serial_total_energy(global2, box, PotentialType::LennardJones, rcut, fprime, fdprime, kappa, alpha);
            Wref = serial_total_virial (global2, box, PotentialType::LennardJones, rcut, fprime, fdprime, kappa, alpha);
            const double rho = (double)global2.size()/box.getV();
            Pref = rho*1.0 + Wref/(2.0*box.getV());
        }
        MPI_Bcast(&Uref,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&Wref,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&Pref,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

        if (extra <= 0.0) {
            // Exactly at rcut → included (neighbor present; energies match ref > 0)
            REQUIRE(any_has == 1);
            if (rank == 0 && found_inclusion) {
                REQUIRE(found_r2 == Approx(rcut*rcut).margin(1e-10));
            }
            REQUIRE(Upar == Approx(Uref).epsilon(1e-12));
            REQUIRE(Wpar == Approx(Wref).epsilon(1e-12));
            REQUIRE(Ppar == Approx(Pref).epsilon(1e-12));
        } else {
            // Slightly beyond rcut → excluded (no neighbor; energy=0)
            // Slightly beyond rcut → excluded (no neighbor; no pair contributions)
            REQUIRE(any_has == 0);
            REQUIRE(Uref == Approx(0.0));
            REQUIRE(Wref == Approx(0.0));

            // Pressure should be the ideal-gas term: P = rho * T
            const double T = 1.0;                  // same T used to build tpar
            const double rho = 2.0 / box.getV();   // 2 particles in the global config
            REQUIRE(Pref == Approx(rho * T).margin(1e-16));

            // Parallel must match serial
            REQUIRE(Upar == Approx(Uref).epsilon(1e-12));
            REQUIRE(Wpar == Approx(Wref).epsilon(1e-12));
            REQUIRE(Ppar == Approx(Pref).epsilon(1e-12));
        }
    };

    // Case A: r = rcut (inclusion)
    run_case(/*extra=*/0.0);
    // Case B: r = rcut + 1e-9 (exclusion)
    run_case(/*extra=*/1e-9);
}


// --- test: incremental ghost updates vs full refresh equivalence ---
TEST_CASE("Incremental halo flush equals full refresh (no extra motion)", "[MPI][ghosts][incremental]") {
    int rank=0, size=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);

    // Box and random particles
    const double Lx = 150.0, Ly = 120.0;
    const int    Nglob = 400;
    const double rcut_want = 3.2;
    const double T = 1.0;
    const PotentialType pot = PotentialType::WCA; // any is fine; choose WCA for variety
    const double f=0.0, alpha=0.0, A0=0.0, kappa=1.0;

    SimulationBox box(Lx, Ly);
    auto decomp = box.bestDecomposition(size);

    // Global random (deterministic)
    auto make_global_uniform = [&](int N, unsigned seed){
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> ux(0.0,Lx), uy(0.0,Ly);
        std::vector<Particle> v; v.reserve(N);
        for (int i=0;i<N;++i) v.push_back(Particle{ux(rng), uy(rng), i});
        return v;
    };
    auto global = make_global_uniform(Nglob, /*seed=*/4242);

    // Local owned subset + per-rank rcut cap
    double x0,x1,y0,y1; decomp.localBounds(rank, Lx, Ly, x0, x1, y0, y1);
    const double Llocx = x1-x0, Llocy = y1-y0;
    auto effective_rcut = [&](double want){
        const double cap = 0.49 * std::min(Llocx, Llocy);
        return std::max(1e-8, std::min(want, cap));
    };
    const double rcut = effective_rcut(rcut_want);

    std::vector<Particle> owned;
    owned.reserve(global.size());
    for (const auto& p : global) if (p.x>=x0 && p.x<x1 && p.y>=y0 && p.y<y1) owned.push_back(p);

    // Build cell list and particle exchange (incremental enabled)
    CellListParallel::Params clp_params{rcut};
    CellListParallel clp(box, decomp, rank, clp_params);
    clp.buildInterior(owned);

    ParticleExchange::Params pex_params{/*enable_incremental=*/true};
    ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pex_params);
    pex.refreshGhosts(owned, clp); // seed

    ThermodynamicCalculatorParallel tpar(MPI_COMM_WORLD, T, pot, rcut,
                                         /*mu*/0.0, f, alpha, A0, kappa);

    // Move all owned slightly and queue incremental updates
    const double dx = 0.07 * rcut;
    const double dy = -0.05 * rcut;
    for (auto& p : owned) {
        p.updatePosition(dx, dy);
        box.applyPBC(p);
        pex.queueGhostUpdateCandidate(p); // queue every moved particle
    }
    clp.buildInterior(owned);

    // Flush incremental halo updates
    pex.flushGhostUpdates(clp);

    // Compute incremental thermodynamics
    const double U_inc = tpar.totalEnergy(owned, box, clp, pex);
    const double W_inc = tpar.totalVirial(owned, box, clp, pex);
    const double P_inc = tpar.pressure(owned, box, clp, pex);

    // Now, WITHOUT moving again, do a full rebuild of ghosts
    pex.refreshGhosts(owned, clp);

    // Compute full-refresh thermodynamics
    const double U_full = tpar.totalEnergy(owned, box, clp, pex);
    const double W_full = tpar.totalVirial(owned, box, clp, pex);
    const double P_full = tpar.pressure(owned, box, clp, pex);

    // Must match tightly (identical configuration; just different update paths)
    REQUIRE(U_inc == Approx(U_full).epsilon(1e-12));
    REQUIRE(W_inc == Approx(W_full).epsilon(1e-12));
    REQUIRE(P_inc == Approx(P_full).epsilon(1e-12));
}
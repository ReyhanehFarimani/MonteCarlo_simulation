#include "catch.hpp"
#include <mpi.h>
#include <vector>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <algorithm>
#include <iomanip>

#include "../mpi_src/particle.h"
#include "../mpi_src/simulation_box.h"
#include "../mpi_src/cell_list_parallel.h"
#include "../mpi_src/particle_exchange.h"
#include "../mpi_src/logging_traj_mpi.h"

// ---------- helpers (same style as yours) ----------

// 9-per-cell layout (same as your current builder)
static std::vector<Particle> build_nine_per_cell(const CellListParallel& clp, int rank)
{
    std::vector<Particle> out;
    out.reserve(clp.nxInterior()*clp.nyInterior()*9);
    int shift = clp.nxInterior()*clp.nyInterior()*10;
    const double x0 = clp.x0(), y0 = clp.y0();
    const double dx = clp.dxCell(), dy = clp.dyCell();
    const double epsx = 0.1 * dx, epsy = 0.1 * dy;

    int id_seed = rank * shift;
    for (int jy=1; jy<=clp.nyInterior(); ++jy) {
        for (int ix=1; ix<=clp.nxInterior(); ++ix) {
            const double cx = x0 + (ix - 0.5)*dx;
            const double cy = y0 + (jy - 0.5)*dy;

            out.push_back(Particle{cx, cy, id_seed++});                 // center
            out.push_back(Particle{cx + (dx/2 - epsx), cy, id_seed++}); // right
            out.push_back(Particle{cx - (dx/2 - epsx), cy, id_seed++}); // left
            out.push_back(Particle{cx, cy + (dy/2 - epsy), id_seed++}); // top
            out.push_back(Particle{cx, cy - (dy/2 - epsy), id_seed++}); // bottom
            out.push_back(Particle{cx + (dx/2 - epsx), cy + (dy/2 - epsy), id_seed++}); // NE
            out.push_back(Particle{cx - (dx/2 - epsx), cy + (dy/2 - epsy), id_seed++}); // NW
            out.push_back(Particle{cx + (dx/2 - epsx), cy - (dy/2 - epsy), id_seed++}); // SE
            out.push_back(Particle{cx - (dx/2 - epsx), cy - (dy/2 - epsy), id_seed++}); // SW
        }
    }
    return out;
}

// Move all owners by (dx,dy) and apply PBC (REAL positions update)
static void move_all(std::vector<Particle>& owned, double dx, double dy, const SimulationBox& box)
{
    for (auto& p : owned) { p.updatePosition(dx, dy); box.applyPBC(p); }
}

// Wrap like your SimulationBox does
static inline double wrap(double a, double L) {
    a -= L * std::floor(a / L);
    if (a >= L) a -= L;
    if (a < 0)  a += L;
    return a;
}

// minimum-image squared distance using box
static inline double min_image_sq(double x1, double y1, double x2, double y2, double Lx, double Ly) {
    double dx = x1 - x2;
    double dy = y1 - y2;
    dx -= Lx * std::round(dx / Lx);
    dy -= Ly * std::round(dy / Ly);
    return dx*dx + dy*dy;
}

// allgather (x,y,id,owner_rank) of owned particles — gives REAL owner
struct OwnedXYZR { double x,y; int id; int owner; };
static std::vector<OwnedXYZR> allgather_owned_real(const std::vector<Particle>& owned)
{
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<OwnedXYZR> local(owned.size());
    for (size_t i=0;i<owned.size();++i) local[i] = OwnedXYZR{owned[i].x, owned[i].y, owned[i].id, rank};

    int nloc = (int)local.size();
    std::vector<int> counts(size,0);
    MPI_Allgather(&nloc, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> displs(size,0);
    int total = 0;
    for (int r=0;r<size;++r){ displs[r]=total; total += counts[r]; }

    std::vector<OwnedXYZR> global(total);

    MPI_Datatype T;
    {
        int bl[4] = {1,1,1,1};
        MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT};
        MPI_Aint disp[4], base;
        OwnedXYZR probe{};
        MPI_Get_address(&probe, &base);
        MPI_Get_address(&probe.x, &disp[0]);
        MPI_Get_address(&probe.y, &disp[1]);
        MPI_Get_address(&probe.id,&disp[2]);
        MPI_Get_address(&probe.owner,&disp[3]);
        for (int k=0;k<4;++k) disp[k] -= base;
        MPI_Type_create_struct(4, bl, disp, types, &T);
        MPI_Type_commit(&T);
    }

    std::vector<int> displs_el = displs;
    MPI_Allgatherv(local.data(), nloc, T,
                   global.data(), counts.data(), displs_el.data(), T,
                   MPI_COMM_WORLD);
    MPI_Type_free(&T);
    return global;
}

// compute neighbor id set (within rcut) of a point (x,y) against global snapshot
static std::unordered_set<int> neighbor_ids_of(double x, double y,
                                               const std::vector<OwnedXYZR>& global,
                                               double Lx, double Ly, double rcutsq,
                                               int self_id)
{
    std::unordered_set<int> ids;
    for (const auto& w : global) {
        if (w.id == self_id) continue;
        double r2 = min_image_sq(x, y, w.x, w.y, Lx, Ly);
        if (r2 <= rcutsq) ids.insert(w.id);
    }
    return ids;
}

// ---------- the test ----------

TEST_CASE("Track ONE REAL particle by ID: neighbors invariant; XYZ logs with REAL owner rank per frame",
          "[MPI][ParticleExchange][track][neighbors][logging][real]")
{
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Geometry
    const double Lx = 100.0, Ly = 100.0;
    const double rcut = 4.0;
    const double step_x = 0.10 * rcut;   // you can change these two to move in y as well
    const double step_y = 0.05 * rcut;

    SimulationBox box(Lx, Ly);
    auto decomp = box.bestDecomposition(size);
    CellListParallel::Params cl_params{rcut};
    CellListParallel clp(box, decomp, rank, cl_params);

    // Build owners and exchange
    auto owned = build_nine_per_cell(clp, rank);
    clp.buildInterior(owned);

    ParticleExchange::Params pex_params{/*enable_incremental=*/false};
    ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pex_params);

    // Logging (owners; rank-coded color)
    LoggingTrajMPI logger(
        /*basename=*/"track_one",
        /*mode=*/LoggingTrajMPI::Mode::MPIIO,
        /*comm=*/MPI_COMM_WORLD,
        /*append=*/false,
        /*rank_as_type=*/true
    );

    // Initial halo + log
    pex.refreshGhosts(owned, clp);
    logger.log_dump(owned, box, /*timestep=*/0);

    // ---- choose the tracked REAL particle (existing id from builder) ----
    // Pick a concrete id that certainly exists on all runs: for rank 0, it's 0;
    // to make it portable across any size, pick id 10 (first cell's "right" on rank 0).
    const int tracked_id = 10;

    // Gather global snapshot at t=0 and find the REAL tracked particle (x0,y0,owner0)
    auto global0 = allgather_owned_real(owned);

    double x0=0.0, y0=0.0;
    int owner0 = -1;
    bool found0 = false;
    for (const auto& w : global0) if (w.id == tracked_id) { x0=w.x; y0=w.y; owner0=w.owner; found0=true; break; }
    REQUIRE(found0);

    const double rcutsq = rcut * rcut;

    // Baseline neighbor id set at t=0 (global truth)
    auto base_neighbors = neighbor_ids_of(x0, y0, global0, Lx, Ly, rcutsq, tracked_id);

    // ---- multi-frame XYZ of the REAL tracked particle (rank 0 writes) ----
    // Each frame:
    //   line 1: "1"
    //   line 2: "step=... Lx=... Ly=... owner_rank=..."
    //   line 3: "<owner_rank_1based>  x  y  0"
    std::ofstream xyz;
    if (rank == 0) {
        xyz.open("track_one_traj.xyz", std::ios::out | std::ios::trunc);
        REQUIRE(!!xyz);
        xyz << std::setprecision(17);

        xyz << "1\n";
        xyz << "step=0 Lx=" << Lx << " Ly=" << Ly << " owner_rank=" << (owner0+1) << "\n";
        xyz << (owner0+1) << " " << x0 << " " << y0 << " " << 0.0 << "\n";
    }

    // march
    const int nsteps = 2000; // reasonably sized; adjust as you like
    for (int s=1; s<=nsteps; ++s) {
        // REAL movement of all particles
        move_all(owned, step_x, step_y, box);

        // rebin owners locally and refresh ghosts
        clp.buildInterior(owned);
        pex.refreshGhosts(owned, clp);

        // migrate periodically (REAL ownership changes)
        if (s % 20 == 0) {
            int sent = pex.migrate(owned, clp);
            REQUIRE(sent >= 0);
            pex.refreshGhosts(owned, clp);
        }

        // log owners (rank-colored)
        logger.log_dump(owned, box, /*timestep=*/s);

        // REAL global snapshot *after* movement/migration
        auto global_now = allgather_owned_real(owned);

        // find REAL current (x,y,owner) of the tracked particle by ID
        double xt=0.0, yt=0.0;
        int owner_now = -1;
        bool found_now = false;
        for (const auto& w : global_now) if (w.id == tracked_id) { xt=w.x; yt=w.y; owner_now=w.owner; found_now=true; break; }
        REQUIRE(found_now); // must exist somewhere

        // neighbor set now (global truth) — should be identical (uniform motion)
        auto now_neighbors = neighbor_ids_of(xt, yt, global_now, Lx, Ly, rcutsq, tracked_id);
        REQUIRE(now_neighbors.size() == base_neighbors.size());
        for (int id : base_neighbors) {
            REQUIRE(now_neighbors.find(id) != now_neighbors.end());
        }

        // append XYZ frame (rank 0 only) with REAL owner
        if (rank == 0) {
            xyz << "1\n";
            xyz << "step=" << s << " Lx=" << Lx << " Ly=" << Ly
                << " owner_rank=" << (owner_now+1) << "\n";
            xyz << (owner_now+1) << " " << xt << " " << yt << " " << 0.0 << "\n";
        }
    }

    if (rank == 0) {
        xyz.flush();
        xyz.close();
    }

    // global conservation at the end (optional safety)
    int loc_cells = clp.nxInterior() * clp.nyInterior();
    int glob_cells = 0;
    MPI_Allreduce(&loc_cells, &glob_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    int nloc = (int)owned.size(), nglob = 0;
    MPI_Allreduce(&nloc, &nglob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    REQUIRE(nglob == 9 * glob_cells);
}

#include "catch.hpp"
#include <mpi.h>
#include <vector>
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "../mpi_src/particle.h"
#include "../mpi_src/simulation_box.h"
#include "../mpi_src/cell_list_parallel.h"
#include "../mpi_src/particle_exchange.h"
#include "../mpi_src/logging_traj_mpi.h"

// ===================== Helpers =====================

// Build 9 particles per interior cell: center, 4 edges (inset), 4 corners (inset)
static std::vector<Particle> build_nine_per_cell(const CellListParallel& clp,int rank)
{
    std::vector<Particle> out;
    out.reserve(clp.nxInterior()*clp.nyInterior()*9);
    int shift = clp.nxInterior()*clp.nyInterior()*10;
    const double x0 = clp.x0(), y0 = clp.y0();
    const double dx = clp.dxCell(), dy = clp.dyCell();
    const double epsx = 0.1 * dx; // inset so they are strictly inside each cell
    const double epsy = 0.1 * dy;

    int id_seed = rank * shift;
    for (int jy=1; jy<=clp.nyInterior(); ++jy) {
        for (int ix=1; ix<=clp.nxInterior(); ++ix) {
            const double cx = x0 + (ix - 0.5)*dx;
            const double cy = y0 + (jy - 0.5)*dy;

            // center
            out.push_back(Particle{cx, cy, id_seed++});

            // edges (midpoints, inset)
            out.push_back(Particle{cx + (dx/2 - epsx), cy,                 id_seed++}); // right
            out.push_back(Particle{cx - (dx/2 - epsx), cy,                 id_seed++}); // left
            out.push_back(Particle{cx,                 cy + (dy/2 - epsy), id_seed++}); // top
            out.push_back(Particle{cx,                 cy - (dy/2 - epsy), id_seed++}); // bottom

            // corners (inset)
            out.push_back(Particle{cx + (dx/2 - epsx), cy + (dy/2 - epsy), id_seed++}); // NE
            out.push_back(Particle{cx - (dx/2 - epsx), cy + (dy/2 - epsy), id_seed++}); // NW
            out.push_back(Particle{cx + (dx/2 - epsx), cy - (dy/2 - epsy), id_seed++}); // SE
            out.push_back(Particle{cx - (dx/2 - epsx), cy - (dy/2 - epsy), id_seed++}); // SW
        }
    }
    return out;
}

// Move all owners by +step_x in x and apply PBC
static void move_all_plus_x(std::vector<Particle>& owned, double step_x, const SimulationBox& box)
{
    for (auto& p : owned) { p.updatePosition(step_x, 0.0); box.applyPBC(p); }
}

// Expected owner rank from (x,y), using pure math
static int expected_rank(double x, double y, double Lx, double Ly, int Px, int Py)
{
    auto wrap = [](double a, double L){ a = std::fmod(a, L); if (a < 0) a += L; return a; };
    x = wrap(x, Lx); y = wrap(y, Ly);
    const double w = Lx / Px, h = Ly / Py;
    int ix = std::min((int)std::floor(x / w), Px-1);
    int iy = std::min((int)std::floor(y / h), Py-1);
    return iy*Px + ix;
}

// Return true if (ix,iy) is on the one-cell ring
static inline bool is_ring_cell(int ix, int iy, int nx_in, int ny_in) {
    return (ix == 0 || ix == nx_in+1 || iy == 0 || iy == ny_in+1);
}

// Query whether there is a ghost exactly at (xg,yg)
static bool has_ghost_at_point(const CellListParallel& clp,
                               double xg, double yg,
                               const std::vector<Particle>& owned)
{
    auto v = clp.neighborsOfPoint(xg, yg, owned);
    for (const auto& pr : v) {
        int idx = pr.first;
        double r2 = pr.second;
        if (idx < 0 && r2 <= 1e-12) return true;
    }
    return false;
}

// Neighbor ranks of "rank" in a Px×Py periodic grid: E,W,N,S, NE,NW,SE,SW
static std::array<int,8> neighbor_ranks_of(int rank, int Px, int Py)
{
    int rx = rank % Px;
    int ry = rank / Px;
    auto wrap = [](int i, int n){ i%=n; if (i<0) i+=n; return i; };
    auto id = [&](int ix, int iy){ ix=wrap(ix,Px); iy=wrap(iy,Py); return iy*Px+ix; };
    return {
        id(rx+1, ry),   // E
        id(rx-1, ry),   // W
        id(rx,   ry+1), // N
        id(rx,   ry-1), // S
        id(rx+1, ry+1), // NE
        id(rx-1, ry+1), // NW
        id(rx+1, ry-1), // SE
        id(rx-1, ry-1)  // SW
    };
}

// All-gather owned positions (x,y) from all ranks
struct OwnedXY { double x,y; };
static std::vector<OwnedXY> allgather_owned_xy(const std::vector<Particle>& owned)
{
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<OwnedXY> local(owned.size());
    for (size_t i=0;i<owned.size();++i) local[i] = OwnedXY{owned[i].x, owned[i].y};

    int nloc = (int)local.size();
    std::vector<int> counts(size,0);
    MPI_Allgather(&nloc, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> displs(size,0);
    int total = 0;
    for (int r=0;r<size;++r){ displs[r]=total; total += counts[r]; }

    std::vector<OwnedXY> global(total);

    // MPI datatype for OwnedXY
    MPI_Datatype T;
    {
        int bl[2] = {1,1};
        MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};
        MPI_Aint disp[2], base;
        OwnedXY probe{};
        MPI_Get_address(&probe, &base);
        MPI_Get_address(&probe.x, &disp[0]);
        MPI_Get_address(&probe.y, &disp[1]);
        for (int k=0;k<2;++k) disp[k] -= base;
        MPI_Type_create_struct(2, bl, disp, types, &T);
        MPI_Type_commit(&T);
    }

    std::vector<int> displs_el = displs; // element displacements

    MPI_Allgatherv(local.data(), nloc, T,
                   global.data(), counts.data(), displs_el.data(), T,
                   MPI_COMM_WORLD);

    MPI_Type_free(&T);
    return global;
}

void writeGhostsDumpPerRank(const std::string& basename,
                                              const SimulationBox& box,
                                              int timestep, const std::vector<PackedParticle>& ghosts_, 
                                            int rank_)
{
    // Build per-rank filename: <basename>.ghosts.rank####.dump
    std::ostringstream name;
    name << basename << ".ghosts.rank"
         << std::setfill('0') << std::setw(4) << rank_ << ".dump";

    std::ofstream ofs(name.str(), std::ios::out | std::ios::trunc);
    if (!ofs) {
        throw std::runtime_error("ParticleExchange: cannot open ghost dump: " + name.str());
    }

    // LAMMPS-like header
    ofs << "ITEM: TIMESTEP\n" << timestep << "\n";
    ofs << "ITEM: NUMBER OF ATOMS\n" << ghosts_.size() << "\n";
    ofs << "ITEM: BOX BOUNDS pp pp ff\n";
    ofs << 0.0 << " " << box.getLx() << "\n";
    ofs << 0.0 << " " << box.getLy() << "\n";
    ofs << 0.0 << " " << 0.0 << "\n";
    ofs << "ITEM: ATOMS id x y z\n";

    // Emit ghosts (wrapped into [0,L))
    for (const auto& g : ghosts_) {
        const double xw = g.x;
        const double yw = g.y;
        ofs << g.id << " " << std::setprecision(10) << xw << " " << yw << " " << 0.0 << "\n";
    }
    ofs.flush();
}

// ===================== The test =====================

TEST_CASE("ParticleExchange: 9-per-cell march + ghost/migration + MPI-IO logging + ghost-only logging",
          "[MPI][ParticleExchange][march][ghosts][migration][logging]")
{
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Problem setup
    const double Lx = 300.0, Ly = 200.0;
    const double rcut = 10.0;
    const double step = 0.1 * rcut; // 0.2

    SimulationBox box(Lx, Ly);
    auto decomp = box.bestDecomposition(size);
    CellListParallel::Params cl_params{rcut};
    CellListParallel clp(box, decomp, rank, cl_params);

    // Logger: single MPI-IO dump, colored by rank (owners)
    LoggingTrajMPI logger(
        /*basename=*/"march9",
        /*mode=*/LoggingTrajMPI::Mode::MPIIO,
        /*comm=*/MPI_COMM_WORLD,
        /*append=*/false,
        /*rank_as_type=*/true
    );

    // NEW: Ghost-only logger (separate file)
    // LoggingTrajMPI ghost_logger(
    //     /*basename=*/"march9_ghosts",
    //     /*mode=*/LoggingTrajMPI::Mode::MPIIO,
    //     /*comm=*/MPI_COMM_WORLD,
    //     /*append=*/false,
    //     /*rank_as_type=*/true   // type is the logging (receiver) rank
    // );

    // Build owners: 9 per interior cell
    auto owned = build_nine_per_cell(clp, rank);
    clp.buildInterior(owned);

    // Wire ParticleExchange (batch halos)
    ParticleExchange::Params pex_params{/*enable_incremental=*/true};
    ParticleExchange pex(MPI_COMM_WORLD, box, decomp, rank, clp, pex_params);

    // Initial halo + dumps
    pex.refreshGhosts(owned, clp);
    logger.log_dump(owned, box, /*timestep=*/0);


    // March
    const int nsteps = 250; // migrations at 5,10,15,20
    const auto nbrs = neighbor_ranks_of(rank, decomp.Px, decomp.Py);

    for (int s=1; s<=nsteps; ++s) {
        // Move + rebuild + refresh
        move_all_plus_x(owned, step, box);
        clp.buildInterior(owned);
        pex.refreshGhosts(owned, clp);

        // Log owners this step
        logger.log_dump(owned, box, /*timestep=*/s);

        // Build and log ghost-only snapshot for this rank

        // const std::vector<PackedParticle>& g = pex.getGhosts();
        // writeGhostsDumpPerRank("march9", box, /*timestep=*/s, g, rank);

       

        // Ghost completeness (only for 8-neighbor owners that map into our ring)
        auto global = allgather_owned_xy(owned);


        // for (const auto& g : global) {
        //     const int owner = expected_rank(g.x, g.y, Lx, Ly, decomp.Px, decomp.Py);
        //     if (owner == rank) continue;

        //     bool owner_is_neighbor = false;
        //     for (int r : nbrs) if (r == owner) { owner_is_neighbor = true; break; }
        //     if (!owner_is_neighbor) continue;

        //     int ix, iy;
        //     bool ok = clp.mapToLocalCell(g.x, g.y, ix, iy);
        //     REQUIRE(ok);
        //     if (!is_ring_cell(ix, iy, clp.nxInterior(), clp.nyInterior())) continue;

        //     REQUIRE(has_ghost_at_point(clp, g.x, g.y, owned));
        // }

        // Migration every 5 steps
        if (s % 5 == 0) {
            int sent = pex.migrate(owned, clp);
            REQUIRE(sent >= 0);
            // After migration, all locals must belong to this rank by tiling
            for (const auto& p : owned) {
                int exp_r = expected_rank(p.x, p.y, Lx, Ly, decomp.Px, decomp.Py);
                REQUIRE(exp_r == rank);
            }
            pex.refreshGhosts(owned, clp);

            const std::vector<PackedParticle>& g = pex.getGhosts();
            writeGhostsDumpPerRank("march9", box, /*timestep=*/s, g, rank);
            // Log both streams again at the same timestep to visualize ownership jumps
            // logger.log_dump(owned, box, /*timestep=*/s);
            // auto global2 = allgather_owned_xy(owned);
            // auto ghosts_local2 = build_local_ghost_snapshot(clp, owned, global2, box, rank, decomp.Px, decomp.Py);
            // ghost_logger.log_dump(ghosts_local2, box, /*timestep=*/s);
        }
    }

    // Global conservation: owners = 9 * sum_r (nx_in * ny_in)
    int loc_cells = clp.nxInterior() * clp.nyInterior();
    int glob_cells = 0;
    MPI_Allreduce(&loc_cells, &glob_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    int nloc = (int)owned.size(), nglob = 0;
    MPI_Allreduce(&nloc, &nglob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    REQUIRE(nglob == 9 * glob_cells);

    // MPI_Finalize();
}

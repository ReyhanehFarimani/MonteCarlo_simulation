#include "cell_list_parallel.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>

// If your project uses a different header name/path for SimulationBox,
// include it here instead of relying on the forward declarations:
#include "simulation_box.h" // adjust if your header name differs

// ---------------------- ctor / geometry ---------------------------------------

CellListParallel::CellListParallel(const SimulationBox& box,
                                   const SimulationBox::Decomposition& decomp,
                                   int rank,
                                   const Params& p)
: rcut_(p.rcut), rcutsq_(p.rcut * p.rcut)
{
    resetGeometry(box, decomp, rank);
}

void CellListParallel::resetGeometry(const SimulationBox& box,
                                     const SimulationBox::Decomposition& d,
                                     int rank)
{
    // Global box sizes and rank grid from decomposition/box.
    Lx_ = box.getLx();
    Ly_ = box.getLy();
    Px_ = d.Px;
    Py_ = d.Py;

    // My rank coordinates in the Px x Py grid.
    {
        auto cc = d.coordsOf(rank);
        cx_ = cc.first;
        cy_ = cc.second;
    }

    // Local owner subdomain bounds for this rank.
    computeLocalBounds_(box, d, rank);

    // Decide interior cell counts (even) and derived sizes (dx,dy) and totals.
    chooseEvenCells_();

    // Allocate/clear bin storage for current geometry.
    clearBins_();

    // Clear any ghosts carried from prior geometry.
    ghosts_.clear();
}

void CellListParallel::computeLocalBounds_(const SimulationBox& box,
                                           const SimulationBox::Decomposition& d,
                                           int rank)
{
    double x0,x1,y0,y1;
    d.localBounds(rank, box.getLx(), box.getLy(), x0, x1, y0, y1);
    x0_ = x0; x1_ = x1; y0_ = y0; y1_ = y1;
}

void CellListParallel::chooseEvenCells_()
{
    // Interior lengths (owner subdomain only).
    const double Lx_loc = x1_ - x0_;
    const double Ly_loc = y1_ - y0_;

    // Choose counts so that dx,dy >= rcut. Then force to nearest even >= 2.
    auto even_down = [](double Lloc, double r) -> int {
        // n <= floor(Lloc / r), then round DOWN to even, but minimum 2.
        int n = (int)std::floor(Lloc / std::max(r, 1e-16)); // guard divide-by-0
        if (n < 2) n = 2;
        if (n % 2) --n;
        if (n < 2) n = 2;
        return n;
    };

    nx_in_ = even_down(Lx_loc, rcut_);
    ny_in_ = even_down(Ly_loc, rcut_);

    // Total grid adds a one-cell ghost ring on each side in both directions.
    nx_tot_ = nx_in_ + 2;
    ny_tot_ = ny_in_ + 2;

    // Interior cell sizes.
    dx_cell_ = Lx_loc / nx_in_;
    dy_cell_ = Ly_loc / ny_in_;

    // Final safety check: enforce evenness (fail fast if violated).
    if ((nx_in_ % 2) != 0 || (ny_in_ % 2) != 0) {
        throw std::runtime_error("CellListParallel: nx_in_ and ny_in_ must be even.");
    }
}

void CellListParallel::clearBins_()
{
    bins_owned_.assign(nx_tot_ * ny_tot_, {});
    bins_ghost_.assign(nx_tot_ * ny_tot_, {});
}

// ---------------------- owned binning / updates --------------------------------

void CellListParallel::buildInterior(const std::vector<Particle>& owned)
{
    // Only owned bins are cleared here; ghost bins are independent.
    for (auto& v : bins_owned_) v.clear();

    // Bin owned particles into the interior indices [1..nx_in_] x [1..ny_in_].
    for (int i = 0; i < (int)owned.size(); ++i) {
        int ix, iy;
        localCellOfGlobal_(owned[i].x, owned[i].y, ix, iy);

        // Clamp to interior window for robustness against FP on upper edges.
        if (ix < 1) ix = 1; else if (ix > nx_in_) ix = nx_in_;
        if (iy < 1) iy = 1; else if (iy > ny_in_) iy = ny_in_;

        bins_owned_[cid(ix, iy)].push_back(i);
    }
}

void CellListParallel::onAcceptedMove(int pidx, const Particle& oldp, const Particle& newp)
{
    // Compute old/new interior bins (clamped to interior window).
    auto clamp_interior = [&](int& ix, int& iy) {
        if (ix < 1) ix = 1; else if (ix > nx_in_) ix = nx_in_;
        if (iy < 1) iy = 1; else if (iy > ny_in_) iy = ny_in_;
    };

    int ix_old, iy_old, ix_new, iy_new;
    localCellOfGlobal_(oldp.x, oldp.y, ix_old, iy_old);
    localCellOfGlobal_(newp.x, newp.y, ix_new, iy_new);
    clamp_interior(ix_old, iy_old);
    clamp_interior(ix_new, iy_new);

    if (ix_old == ix_new && iy_old == iy_new) return; // same bin → nothing to update

    // Remove from old cell (swap-erase).
    {
        auto& v = bins_owned_[cid(ix_old, iy_old)];
        for (size_t k = 0; k < v.size(); ++k) {
            if (v[k] == pidx) { v[k] = v.back(); v.pop_back(); break; }
        }
    }

    // Add to new cell.
    bins_owned_[cid(ix_new, iy_new)].push_back(pidx);
}

// ---------------------- parity & random selection ------------------------------

std::vector<std::pair<int,int>>
CellListParallel::cellsWithParity(Parity par) const
{
    std::vector<std::pair<int,int>> out;
    out.reserve((nx_in_ * ny_in_) / 4);
    // Interior indices are 1..nx_in_, 1..ny_in_. We map them to 0-based for parity:
    // (ix-1, iy-1) → parity in {EvenEven, EvenOdd, OddEven, OddOdd}.
    for (int iy = 1; iy <= ny_in_; ++iy) {
        for (int ix = 1; ix <= nx_in_; ++ix) {
            const bool ex = ((ix - 1) % 2) == 0;
            const bool ey = ((iy - 1) % 2) == 0;
            Parity cur = ex ? (ey ? Parity::EvenEven : Parity::EvenOdd)
                            : (ey ? Parity::OddEven  : Parity::OddOdd);
            if (cur == par) out.emplace_back(ix, iy);
        }
    }
    return out;
}

int CellListParallel::randomOwnedInCell(int ix, int iy, std::mt19937_64& rng) const
{
    const auto& v = bins_owned_[cid(ix, iy)];
    if (v.empty()) return -1;
    std::uniform_int_distribution<size_t> dist(0, v.size() - 1);
    return v[dist(rng)];
}

// ---------------------- ghost layer --------------------------------------------

void CellListParallel::clearGhosts()
{
    ghosts_.clear();
    for (auto& v : bins_ghost_) v.clear();
}

void CellListParallel::addGhost(const Particle& g)
{
    ghosts_.push_back(g);
}

void CellListParallel::rebuildGhostBins()
{
    for (auto& v : bins_ghost_) v.clear();

    // Place ghosts strictly in the ring cells:
    //   ix==0 or ix==nx_in_+1 or iy==0 or iy==ny_in_+1
    for (int gi = 0; gi < (int)ghosts_.size(); ++gi) {
        int ix, iy;
        localCellOfGlobal_(ghosts_[gi].x, ghosts_[gi].y, ix, iy);
        const bool on_ring = (ix == 0 || ix == nx_in_ + 1 ||
                              iy == 0 || iy == ny_in_ + 1);
        if (on_ring) bins_ghost_[cid(ix, iy)].push_back(gi);
        // If a ghost lands interior, it likely reflects a too-wide halo send.
    }
}

// ---------------------- neighbor queries ---------------------------------------

std::vector<std::pair<int,double>>
CellListParallel::neighborsOfOwned(int localOwnedIdx,
                                   const std::vector<Particle>& owned) const
{
    std::vector<std::pair<int,double>> out;
    if (localOwnedIdx < 0 || localOwnedIdx >= (int)owned.size()) return out;

    // Starting cell for the owned particle, clamped to interior.
    int ix, iy;
    localCellOfGlobal_(owned[localOwnedIdx].x, owned[localOwnedIdx].y, ix, iy);
    if (ix < 1) ix = 1; else if (ix > nx_in_) ix = nx_in_;
    if (iy < 1) iy = 1; else if (iy > ny_in_) iy = ny_in_;

    const double xi = owned[localOwnedIdx].x;
    const double yi = owned[localOwnedIdx].y;

    // Search 3x3 neighborhood around (ix,iy), including ring cells at the edges.
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            const int jx = ix + dx;
            const int jy = iy + dy;
            if (jx < 0 || jx >= nx_tot_ || jy < 0 || jy >= ny_tot_) continue;

            // Interior → owned neighbors in that bucket.
            if (jx >= 1 && jx <= nx_in_ && jy >= 1 && jy <= ny_in_) {
                const auto& v = bins_owned_[cid(jx, jy)];
                for (int jloc : v) {
                    if (jloc == localOwnedIdx) continue;
                    double ddx = owned[jloc].x - xi;
                    double ddy = owned[jloc].y - yi;
                    min_image(ddx, ddy, ddx, ddy);
                    const double r2 = ddx*ddx + ddy*ddy;
                    if (r2 <= rcutsq_) out.emplace_back(jloc, r2);
                }
            } else {
                // Ghost ring → check ghosts in this bucket.
                const auto& vg = bins_ghost_[cid(jx, jy)];
                for (int gi : vg) {
                    const auto& g = ghosts_[gi];
                    double ddx = g.x - xi;
                    double ddy = g.y - yi;
                    min_image(ddx, ddy, ddx, ddy);
                    const double r2 = ddx*ddx + ddy*ddy;
                    if (r2 <= rcutsq_) out.emplace_back(-1 - gi, r2); // encode ghosts as negative ids
                }
            }
        }
    }
    return out;
}

std::vector<std::pair<int,double>>
CellListParallel::neighborsOfPoint(double xg, double yg,
                                   const std::vector<Particle>& owned) const
{
    std::vector<std::pair<int,double>> out;

    // Locate the total-grid cell (including ring). Clamp into [0..nx_tot_-1],[0..ny_tot_-1].
    int ix, iy;
    localCellOfGlobal_(xg, yg, ix, iy);
    ix = std::min(std::max(ix, 0), nx_tot_ - 1);
    iy = std::min(std::max(iy, 0), ny_tot_ - 1);

    // 3x3 search around (ix,iy).
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            const int jx = ix + dx;
            const int jy = iy + dy;
            if (jx < 0 || jx >= nx_tot_ || jy < 0 || jy >= ny_tot_) continue;

            if (jx >= 1 && jx <= nx_in_ && jy >= 1 && jy <= ny_in_) {
                const auto& v = bins_owned_[cid(jx, jy)];
                for (int jloc : v) {
                    double ddx = owned[jloc].x - xg;
                    double ddy = owned[jloc].y - yg;
                    min_image(ddx, ddy, ddx, ddy);
                    const double r2 = ddx*ddx + ddy*ddy;
                    if (r2 <= rcutsq_) out.emplace_back(jloc, r2);
                }
            } else {
                const auto& vg = bins_ghost_[cid(jx, jy)];
                for (int gi : vg) {
                    const auto& g = ghosts_[gi];
                    double ddx = g.x - xg;
                    double ddy = g.y - yg;
                    min_image(ddx, ddy, ddx, ddy);
                    const double r2 = ddx*ddx + ddy*ddy;
                    if (r2 <= rcutsq_) out.emplace_back(-1 - gi, r2);
                }
            }
        }
    }
    return out;
}

// ---------------------- mapping & geometry helpers -----------------------------

bool CellListParallel::mapToLocalCell(double xg, double yg, int& ix, int& iy) const
{
    localCellOfGlobal_(xg, yg, ix, iy);
    return (ix >= 0 && ix < nx_tot_ && iy >= 0 && iy < ny_tot_);
}

void CellListParallel::min_image(double dx, double dy, double& rdx, double& rdy) const
{
    // Minimum-image displacement in global periodic box [0,Lx) x [0,Ly).
    rdx = dx - Lx_ * std::round(dx / Lx_);
    rdy = dy - Ly_ * std::round(dy / Ly_);
}

void CellListParallel::localCellOfGlobal_(double xg, double yg, int& ix, int& iy) const
{
    // We want to map a global point into the local *total* grid:
    //   X-cells: 0 .. nx_in_+1  (0 and nx_in_+1 are ghost ring in x)
    //   Y-cells: 0 .. ny_in_+1  (0 and ny_in_+1 are ghost ring in y)
    //
    // Step 1) Define the extended window covering owner + 1 ghost cell on each side:
    const double gxw = dx_cell_; // one-cell ring width in x
    const double gyw = dy_cell_; // one-cell ring width in y

    const double left   = x0_ - gxw;
    const double right  = x1_ + gxw;
    const double bottom = y0_ - gyw;
    const double top    = y1_ + gyw;

    // Step 2) Wrap (xg,yg) into that window with global periodicity.
    auto wrap_into = [](double val, double lo, double hi, double L) {
        const double span = hi - lo;
        // Move by +/-L until in [lo, hi). Using while to be robust.
        while (val <  lo) val += L;
        while (val >= hi) val -= L;
        // Guard for rare FP: clamp if marginally outside.
        if (val < lo)  val = lo;
        if (val >= hi) val = std::nextafter(hi, lo);
        return val;
    };

    const double xw = wrap_into(xg, left,  right,  Lx_);
    const double yw = wrap_into(yg, bottom, top,    Ly_);

    // Step 3) Classify x into ring/interior cell index.
    if (xw < x0_)        ix = 0;
    else if (xw >= x1_)  ix = nx_in_ + 1;
    else {
        const double fx = (xw - x0_) / dx_cell_;      // in [0, nx_in_)
        int ixi = 1 + (int)std::floor(fx);            // map to 1..nx_in_
        if (ixi < 1) ixi = 1;
        if (ixi > nx_in_) ixi = nx_in_;
        ix = ixi;
    }

    // Step 4) Classify y similarly.
    if (yw < y0_)        iy = 0;
    else if (yw >= y1_)  iy = ny_in_ + 1;
    else {
        const double fy = (yw - y0_) / dy_cell_;
        int iyi = 1 + (int)std::floor(fy);
        if (iyi < 1) iyi = 1;
        if (iyi > ny_in_) iyi = ny_in_;
        iy = iyi;
    }

    // Final clamp (numerical safety).
    ix = std::max(0, std::min(nx_in_ + 1, ix));
    iy = std::max(0, std::min(ny_in_ + 1, iy));
}

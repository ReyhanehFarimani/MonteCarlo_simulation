#include "cell_list_parallel.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

static inline int pos_mod_int(int a, int m) {
    int r = a % m; return (r < 0) ? (r + m) : r;
}

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
    Lx_ = box.getLx();
    Ly_ = box.getLy();
    Px_ = d.Px;
    Py_ = d.Py;

    // get my coords
    auto cc = d.coordsOf(rank);
    cx_ = cc.first;
    cy_ = cc.second;

    computeLocalBounds_(box, d, rank);
    chooseEvenCells_();
    clearBins_();
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
    const double Lx_loc = x1_ - x0_;
    const double Ly_loc = y1_ - y0_;

    auto even_down = [](double Lloc, double r) -> int {
        // Max cells so that dx >= rcut  => n <= floor(L/r)
        int n = (int)std::floor(Lloc / r);
        // enforce even by rounding DOWN to nearest even
        if (n % 2) --n;
        // minimum 2 cells for robustness
        if (n < 2) n = 2;
        return n;
    };

    nx_in_ = even_down(Lx_loc, rcut_);
    ny_in_ = even_down(Ly_loc, rcut_);

    nx_tot_ = nx_in_ + 2;
    ny_tot_ = ny_in_ + 2;

    dx_cell_ = Lx_loc / nx_in_;
    dy_cell_ = Ly_loc / ny_in_;
}


void CellListParallel::clearBins_()
{
    bins_owned_.assign(nx_tot_ * ny_tot_, {});
    bins_ghost_.assign(nx_tot_ * ny_tot_, {});
}

void CellListParallel::buildInterior(const std::vector<Particle>& owned)
{
    // clear only owned bins
    for (auto& v : bins_owned_) v.clear();

    // assign owned particles to interior cells [1..nx_in_]x[1..ny_in_]
    for (int i = 0; i < (int)owned.size(); ++i) {
        int ix, iy;
        localCellOfGlobal_(owned[i].x, owned[i].y, ix, iy);
        // clamp into interior in case of round-off at upper edge
        if (ix < 1) ix = 1; else if (ix > nx_in_) ix = nx_in_;
        if (iy < 1) iy = 1; else if (iy > ny_in_) iy = ny_in_;
        bins_owned_[cid(ix, iy)].push_back(i);
    }
}

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
    for (int gi = 0; gi < (int)ghosts_.size(); ++gi) {
        int ix, iy;
        localCellOfGlobal_(ghosts_[gi].x, ghosts_[gi].y, ix, iy);
        // keep ghosts only if they fall in the ghost ring (0 or nx_in_+1, 0 or ny_in_+1)
        if (ix == 0 || ix == nx_in_ + 1 || iy == 0 || iy == ny_in_ + 1) {
            bins_ghost_[cid(ix, iy)].push_back(gi);
        }
        // else: it's interior; your halo exchange should normally not send those.
    }
}

std::vector<std::pair<int,double>>
CellListParallel::neighborsOfOwned(int localOwnedIdx,
                                   const std::vector<Particle>& owned) const
{
    std::vector<std::pair<int,double>> out;
    if (localOwnedIdx < 0 || localOwnedIdx >= (int)owned.size()) return out;

    // find the owned particle cell
    int ix, iy;
    localCellOfGlobal_(owned[localOwnedIdx].x, owned[localOwnedIdx].y, ix, iy);
    if (ix < 1) ix = 1; else if (ix > nx_in_) ix = nx_in_;
    if (iy < 1) iy = 1; else if (iy > ny_in_) iy = ny_in_;

    const double xi = owned[localOwnedIdx].x;
    const double yi = owned[localOwnedIdx].y;

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            const int jx = ix + dx;
            const int jy = iy + dy;
            if (jx < 0 || jx >= nx_tot_ || jy < 0 || jy >= ny_tot_) continue;

            // owned neighbors inside interior window
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
                // ghost ring
                const auto& vg = bins_ghost_[cid(jx, jy)];
                for (int gi : vg) {
                    const auto& g = ghosts_[gi];
                    double ddx = g.x - xi;
                    double ddy = g.y - yi;
                    min_image(ddx, ddy, ddx, ddy);
                    const double r2 = ddx*ddx + ddy*ddy;
                    if (r2 <= rcutsq_) out.emplace_back(-1 - gi, r2); // encode ghosts as negative ids if you like
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

    int ix, iy;
    localCellOfGlobal_(xg, yg, ix, iy);
    // if strictly outside the subdomain+ghost by numerical issues, clamp to nearest border cell
    ix = std::min(std::max(ix, 0), nx_tot_ - 1);
    iy = std::min(std::max(iy, 0), ny_tot_ - 1);

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

bool CellListParallel::mapToLocalCell(double xg, double yg, int& ix, int& iy) const
{
    localCellOfGlobal_(xg, yg, ix, iy);
    // permit ghost ring [0..nx_in+1], [0..ny_in+1]
    if (ix < 0 || ix >= nx_tot_ || iy < 0 || iy >= ny_tot_) return false;
    return true;
}

void CellListParallel::min_image(double dx, double dy, double& rdx, double& rdy) const
{
    // minimum image on the GLOBAL box (not only subdomain), consistent with your SimulationBox
    rdx = dx - Lx_ * std::round(dx / Lx_);
    rdy = dy - Ly_ * std::round(dy / Ly_);
}

void CellListParallel::wrapIntoSubdomain_(double& x) const
{
    // map global x into [x0_, x1_) by adding/subtracting multiples of Lx_
    const double w = x - x0_;
    const double Lloc = x1_ - x0_;
    double k = std::floor(w / Lloc);
    x -= k * Lloc;
    // bring back to global frame relative to x0_
    x = x0_ + (w - k * Lloc);
}

void CellListParallel::wrapIntoSubdomainY_(double& y) const
{
    const double w = y - y0_;
    const double Lloc = y1_ - y0_;
    double k = std::floor(w / Lloc);
    y = y0_ + (w - k * Lloc);
}
void CellListParallel::localCellOfGlobal_(double xg, double yg, int& ix, int& iy) const
{
    const double Llocx = x1_ - x0_;
    const double Llocy = y1_ - y0_;

    // Ghost ring widths (one cell)
    const double gx = dx_cell_;
    const double gy = dy_cell_;

    // Target window: [x0_-gx, x1_+gx) × [y0_-gy, y1_+gy)
    const double left   = x0_ - gx;
    const double right  = x1_ + gx;
    const double bottom = y0_ - gy;
    const double top    = y1_ + gy;

    // Wrap by the GLOBAL box so that (x,y) lands inside the window above
    auto wrap_into = [](double val, double lo, double hi, double L) {
        // bring into [lo, hi)
        const double span = hi - lo;           // here span = Lloc + 2*ghost_cell
        // shift by multiples of L until inside
        while (val <  lo) val += L;
        while (val >= hi) val -= L;
        // very rare: if span < L and val still escapes due to FP, clamp
        if (val < lo)  val = lo;
        if (val >= hi) val = std::nextafter(hi, lo);
        return val;
    };

    double xw = wrap_into(xg, left,  right,  Lx_);
    double yw = wrap_into(yg, bottom, top,    Ly_);

    // Classify X
    if (xw < x0_)        ix = 0;
    else if (xw >= x1_)  ix = nx_in_ + 1;
    else {
        const double fx = (xw - x0_) / dx_cell_;   // in [0, nx_in)
        int ixi = 1 + (int)std::floor(fx);
        // numeric guard (xw very close to x1_)
        if (ixi < 1) ixi = 1;
        if (ixi > nx_in_) ixi = nx_in_;
        ix = ixi;
    }

    // Classify Y
    if (yw < y0_)        iy = 0;
    else if (yw >= y1_)  iy = ny_in_ + 1;
    else {
        const double fy = (yw - y0_) / dy_cell_;
        int iyi = 1 + (int)std::floor(fy);
        if (iyi < 1) iyi = 1;
        if (iyi > ny_in_) iyi = ny_in_;
        iy = iyi;
    }

    // Final clamp to [0..nx_in+1], [0..ny_in+1]
    ix = std::max(0, std::min(nx_in_ + 1, ix));
    iy = std::max(0, std::min(ny_in_ + 1, iy));
}

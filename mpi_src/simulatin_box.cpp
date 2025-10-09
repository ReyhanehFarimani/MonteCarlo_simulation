#include "simulation_box.h"
#include <cmath>
#include <stdexcept>
#include <limits>

SimulationBox::SimulationBox(double Lx, double Ly)
    : Lx_(Lx), Ly_(Ly)
{
    if (Lx_ <= 0.0 || Ly_ <= 0.0) throw std::invalid_argument("Box lengths must be > 0");
    invLx_ = 1.0 / Lx_;
    invLy_ = 1.0 / Ly_;
    V_     = Lx_ * Ly_;
}

void SimulationBox::applyPBC(Particle &p) const {
    p.x -= Lx_ * std::floor(p.x * invLx_);
    p.y -= Ly_ * std::floor(p.y * invLy_);
    if (p.x >= Lx_) p.x -= Lx_;
    if (p.y >= Ly_) p.y -= Ly_;
    if (p.x <  0.0) p.x += Lx_;
    if (p.y <  0.0) p.y += Ly_;
}

void SimulationBox::applyPBC(double &x, double &y) const {
    x -= Lx_ * std::floor(x * invLx_);
    y -= Ly_ * std::floor(y * invLy_);
    if (x >= Lx_) x -= Lx_;
    if (y >= Ly_) y -= Ly_;
    if (x <  0.0) x += Lx_;
    if (y <  0.0) y += Ly_;
}

double SimulationBox::minimumImageDistanceSquared(const Particle &p1, const Particle &p2) const {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    dx -= Lx_ * std::round(dx * invLx_);
    dy -= Ly_ * std::round(dy * invLy_);
    return dx*dx + dy*dy;
}

double SimulationBox::minimumImageDistance(const Particle &p1, const Particle &p2) const {
    return std::sqrt(minimumImageDistanceSquared(p1, p2));
}

void SimulationBox::setLx(double lx) {
    if (lx <= 0.0) throw std::invalid_argument("Lx must be > 0");
    Lx_ = lx; invLx_ = 1.0 / Lx_; V_ = Lx_ * Ly_;
}

void SimulationBox::setLy(double ly) {
    if (ly <= 0.0) throw std::invalid_argument("Ly must be > 0");
    Ly_ = ly; invLy_ = 1.0 / Ly_; V_ = Lx_ * Ly_;
}

void SimulationBox::setV(double v) {
    if (v <= 0.0) throw std::invalid_argument("V must be > 0");
    const double aspect = Lx_ / Ly_;
    Ly_    = std::sqrt(v / aspect);
    Lx_    = aspect * Ly_;
    invLx_ = 1.0 / Lx_;
    invLy_ = 1.0 / Ly_;
    V_     = v;
}

void SimulationBox::recenter(Particle &p, double lx_old, double ly_old) const {
    if (lx_old <= 0.0 || ly_old <= 0.0) throw std::invalid_argument("old box lengths must be > 0");
    double fx = p.x / lx_old;
    double fy = p.y / ly_old;
    p.x = fx * Lx_;
    p.y = fy * Ly_;
    applyPBC(p);
}

// ---- Decomposition helpers ----

int SimulationBox::Decomposition::rankOf(int cx, int cy) const {
    cx = (cx % Px + Px) % Px;
    cy = (cy % Py + Py) % Py;
    return cy * Px + cx;
}

std::pair<int,int> SimulationBox::Decomposition::coordsOf(int rank) const {
    return { rank % Px, rank / Px };
}

int SimulationBox::Decomposition::west (int rank) const { auto [cx,cy]=coordsOf(rank); return rankOf(cx-1, cy  ); }
int SimulationBox::Decomposition::east (int rank) const { auto [cx,cy]=coordsOf(rank); return rankOf(cx+1, cy  ); }
int SimulationBox::Decomposition::south(int rank) const { auto [cx,cy]=coordsOf(rank); return rankOf(cx  , cy-1); }
int SimulationBox::Decomposition::north(int rank) const { auto [cx,cy]=coordsOf(rank); return rankOf(cx  , cy+1); }

void SimulationBox::Decomposition::localBounds(int rank, double Lx, double Ly,
                                               double &x0, double &x1, double &y0, double &y1) const {
    auto [cx, cy] = coordsOf(rank);
    const double dx = Lx / Px;
    const double dy = Ly / Py;
    x0 = cx * dx; x1 = (cx + 1) * dx;
    y0 = cy * dy; y1 = (cy + 1) * dy;
}

SimulationBox::Decomposition SimulationBox::bestDecomposition(int nranks) const {
    if (nranks <= 0) throw std::invalid_argument("nranks must be > 0");

    int bestPx = 1, bestPy = nranks;
    double bestScore = std::numeric_limits<double>::infinity();

    auto evaluate = [&](int Px, int Py) {
        const double cell_aspect = (Lx_ / Px) / (Ly_ / Py);
        const double aspect_dev  = std::abs(std::log(cell_aspect));
        const double balance     = std::abs(Px - Py) * 1e-3;
        return aspect_dev + balance;
    };

    for (int px = 1; px * px <= nranks; ++px) {
        if (nranks % px) continue;
        const int py = nranks / px;

        const double s1 = evaluate(px, py);
        if (s1 < bestScore) { bestScore = s1; bestPx = px; bestPy = py; }

        const double s2 = evaluate(py, px);
        if (s2 < bestScore) { bestScore = s2; bestPx = py; bestPy = px; }
    }

    Decomposition d;
    d.Px = bestPx; d.Py = bestPy;
    d.dx = Lx_ / d.Px; d.dy = Ly_ / d.Py;
    return d;
}

SimulationBox::Decomposition SimulationBox::bestDecompositionConstrained(
    int nranks, double min_dx, double min_dy, bool strict) const
{
    if (nranks <= 0) throw std::invalid_argument("nranks must be > 0");
    if (min_dx <= 0.0 || min_dy <= 0.0) throw std::invalid_argument("min_dx/min_dy must be > 0");

    int bestPx = 1, bestPy = nranks;
    double bestScore = std::numeric_limits<double>::infinity();
    bool found = false;

    auto evaluate = [&](int Px, int Py) {
        const double dx = Lx_ / Px, dy = Ly_ / Py;
        if (dx < min_dx || dy < min_dy) return std::numeric_limits<double>::infinity();
        const double cell_aspect = (Lx_ / Px) / (Ly_ / Py);
        const double aspect_dev  = std::abs(std::log(cell_aspect));
        const double balance     = std::abs(Px - Py) * 1e-3;
        return aspect_dev + balance;
    };

    for (int px = 1; px * px <= nranks; ++px) {
        if (nranks % px) continue;
        const int py = nranks / px;

        const double s1 = evaluate(px,  py);
        if (s1 < bestScore) { bestScore = s1; bestPx = px; bestPy = py; found = std::isfinite(s1); }

        const double s2 = evaluate(py,  px);
        if (s2 < bestScore) { bestScore = s2; bestPx = py; bestPy = px; found = std::isfinite(s2); }
    }

    if (!found) {
        if (strict) {
            throw std::runtime_error("No (Px,Py) satisfies min_dx/min_dy.");
        } // else fall through with best unconstrained (which will be 1 x nranks)
    }

    Decomposition d;
    d.Px = bestPx; d.Py = bestPy;
    d.dx = Lx_ / d.Px; d.dy = Ly_ / d.Py;
    return d;
}

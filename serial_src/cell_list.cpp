#include "cell_list.h"
#include <cmath>

CellList::CellList(const SimulationBox& box, double rcut)
    : box_(box), rcut_(rcut), rcutsq_(rcut*rcut) {
    // choose number of cells such that cell size >= rcut
    nx_ = static_cast<int>(box_.getLx() / rcut_);
    ny_ = static_cast<int>(box_.getLy() / rcut_);
    if (nx_ < 1) nx_ = 1;
    if (ny_ < 1) ny_ = 1;
    cell_dx_ = box_.getLx() / nx_;
    cell_dy_ = box_.getLy() / ny_;
    cells_.resize(nx_ * ny_);
}

void CellList::build(const std::vector<Particle>& particles) {
    int N = particles.size();
    particle_cell_.assign(N, 0);
    // clear cells
    for (auto& cell : cells_) cell.clear();
    // assign particles to cells
    for (int i = 0; i < N; ++i) {
        double x = particles[i].x;
        double y = particles[i].y;
        int idx = cellIndex(x, y);
        particle_cell_[i] = idx;
        cells_[idx].push_back(i);
    }
}

std::vector<std::pair<int,double>> CellList::getNeighbors(int idx, const std::vector<Particle>& particles) const {
    std::vector<std::pair<int,double>> neighbors;
    int ci = particle_cell_[idx];
    int cx = ci % nx_;
    int cy = ci / nx_;

    // loop over this cell and 8 neighbors for periodic cells
    for (int ddx = -1; ddx <= 1; ++ddx) {
        for (int ddy = -1; ddy <= 1; ++ddy) {
            int ncx = (cx + ddx + nx_) % nx_;
            int ncy = (cy + ddy + ny_) % ny_;
            int nci = ncx + ncy * nx_;
            for (int j : cells_[nci]) {
                if (j == idx) continue;
                // compute dx, dy with PBC
                double dx = particles[j].x - particles[idx].x;
                double dy = particles[j].y - particles[idx].y;
                dx -= box_.getLx() * std::round(dx/box_.getLx());
                dy -= box_.getLy() * std::round(dy/box_.getLy());
                double r_sq = dx*dx + dy*dy;
                if (r_sq <= rcutsq_) {
                    neighbors.emplace_back(j, r_sq);
                }
            }
        }
    }
    return neighbors;
}

int CellList::cellIndex(double x, double y) const {
    int cx = static_cast<int>(std::floor(x / cell_dx_));
    int cy = static_cast<int>(std::floor(y / cell_dy_));
    cx = (cx % nx_ + nx_) % nx_;
    cy = (cy % ny_ + ny_) % ny_;
    return cx + cy * nx_;
}
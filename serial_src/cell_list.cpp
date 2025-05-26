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

std::vector<std::pair<int,double>> CellList::getNeighbors2(const Particle& p, const std::vector<Particle>& particles) const {
    std::vector<std::pair<int,double>> neighbors;
    double x = p.x;
    double y = p.y;
    int cx = static_cast<int>(std::floor(x / cell_dx_));
    int cy = static_cast<int>(std::floor(y / cell_dy_));
    cx = (cx % nx_ + nx_) % nx_;
    cy = (cy % ny_ + ny_) % ny_;

    // loop over this cell and 8 neighbors for periodic cells
    for (int ddx = -1; ddx <= 1; ++ddx) {
        for (int ddy = -1; ddy <= 1; ++ddy) {
            int ncx = (cx + ddx + nx_) % nx_;
            int ncy = (cy + ddy + ny_) % ny_;
            int nci = ncx + ncy * nx_;
            for (int j : cells_[nci]) {
                // compute dx, dy with PBC
                double dx = particles[j].x - x;
                double dy = particles[j].y - y;
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

void CellList::addParticle(const Particle& p, int index) {
    int cellIdx = cellIndex(p.x, p.y);
    if (index >= particle_cell_.size()) {
        particle_cell_.resize(index + 1);  // accommodate new index
    }
    particle_cell_[index] = cellIdx;
    cells_[cellIdx].push_back(index);
}

void CellList::removeParticle(int index) {
    if (index < 0 || index >= static_cast<int>(particle_cell_.size()))
        return;

    int cellIdx = particle_cell_[index];

    // 1. Remove from the appropriate cell
    auto& cell = cells_[cellIdx];
    cell.erase(std::remove(cell.begin(), cell.end(), index), cell.end());

    // 2. Fix remaining indices in all cells
    for (auto& c : cells_) {
        for (int& i : c) {
            if (i > index) {
                i -= 1;
            }
        }
    }

    // 3. Fix particle_cell_ entries
    for (size_t i = index + 1; i < particle_cell_.size(); ++i) {
        particle_cell_[i - 1] = particle_cell_[i];
    }
    particle_cell_.pop_back();  // remove the extra entry at the end
}

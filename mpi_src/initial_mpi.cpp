#include "initial_mpi.h"
#include <cmath>
#include <algorithm>

int compute_local_count(int N_global, const SimulationBox::Decomposition& d, int rank) {
    const int P = d.Px * d.Py;
    const int base = (P > 0) ? (N_global / P) : 0;
    const int rem  = (P > 0) ? (N_global % P) : 0;
    // Deterministic: first 'rem' ranks get +1
    return base + (rank < rem ? 1 : 0);
}

static void local_bounds(const SimulationBox& box,
                         const SimulationBox::Decomposition& d,
                         int rank,
                         double& x0, double& x1, double& y0, double& y1)
{
    d.localBounds(rank, box.getLx(), box.getLy(), x0, x1, y0, y1);
}

void initializeParticles_rank(std::vector<Particle>& particles,
                              const SimulationBox& box,
                              const SimulationBox::Decomposition& decomp,
                              int rank,
                              int N_local,
                              bool random,
                              unsigned seed,
                              long long id_base)
{
    particles.clear();
    if (N_local <= 0) return;
    particles.reserve(static_cast<size_t>(N_local));

    double x0, x1, y0, y1;
    local_bounds(box, decomp, rank, x0, x1, y0, y1);
    const double dx = x1 - x0;
    const double dy = y1 - y0;

    if (random) {
        // Reproducible, rank-unique RNG
        RNG_parallel rng(seed, rank);
        for (int i = 0; i < N_local; ++i) {
            const double xr = rng.uniform01(); // in [0,1)
            const double yr = rng.uniform01();
            Particle p;
            p.id = static_cast<int>(id_base + i);
            p.x  = x0 + xr * dx;
            p.y  = y0 + yr * dy;
            particles.emplace_back(p);
        }
    } else {
        // Grid tile inside this subdomain: as square as possible
        // Choose nx, ny such that nx*ny >= N_local and nx ~ dx/dy * sqrt(N_local)
        const double aspect = (dy > 0.0) ? dx / dy : 1.0;
        const double rootN  = std::sqrt(static_cast<double>(N_local));
        int nx = std::max(1, static_cast<int>(std::round(aspect * rootN)));
        int ny = std::max(1, (N_local + nx - 1) / nx); // ceil(N_local / nx)

        // spacing (half-open [x0,x1), [y0,y1)), center particles in cells
        const double sx = dx / nx;
        const double sy = dy / ny;

        int count = 0;
        for (int iy = 0; iy < ny && count < N_local; ++iy) {
            for (int ix = 0; ix < nx && count < N_local; ++ix) {
                Particle p;
                p.id = static_cast<int>(id_base + count);
                p.x  = x0 + (ix + 0.5) * sx; // center
                p.y  = y0 + (iy + 0.5) * sy;
                // Apply global PBC guard (it should already be in-range)
                // but ensure [0,L) if floating-point edge cases occur:
                double gx = p.x, gy = p.y;
                const_cast<SimulationBox&>(box).applyPBC(gx, gy);
                p.x = gx; p.y = gy;
                particles.emplace_back(p);
                ++count;
            }
        }
    }
}

void initializeParticles_globalN(std::vector<Particle>& particles,
                                 const SimulationBox& box,
                                 const SimulationBox::Decomposition& decomp,
                                 int rank,
                                 int N_global,
                                 bool random,
                                 unsigned seed,
                                 long long id_stride)
{
    const int n_local = compute_local_count(N_global, decomp, rank);
    // Give each rank a disjoint id range for clarity: [rank*id_stride, rank*id_stride + n_local)
    const long long id_base = static_cast<long long>(rank) * id_stride;
    initializeParticles_rank(particles, box, decomp, rank, n_local, random, seed, id_base);
}

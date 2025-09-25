#ifndef INITIAL_MPI_H
#define INITIAL_MPI_H

#include <vector>
#include "particle.h"
#include "simulation_box.h"
#include "rng_parallel.h"

/**
 * @brief Compute a deterministic local particle count for this rank
 *        from a global count N_global over a Px*Py grid.
 *
 * Rule: base = N_global / (Px*Py); remainder r = N_global % (Px*Py);
 *       ranks 0..r-1 get (base+1), others get base.
 */
int compute_local_count(int N_global, const SimulationBox::Decomposition& d, int rank);

/**
 * @brief Initialize N_local particles for this rank inside its subdomain (global coords).
 *        If random=true, uniform in [x0,x1)×[y0,y1) using RNG_parallel(seed,rank).
 *        Else, place on a grid that tiles the local subdomain.
 *
 * @param particles  output vector (cleared)
 * @param box        global simulation box
 * @param decomp     chosen decomposition (Px,Py, etc.)
 * @param rank       this rank id (row-major: rank=cy*Px+cx)
 * @param N_local    number of particles to place on this rank
 * @param random     random or grid
 * @param seed       global seed for RNG (used only if random==true)
 * @param id_base    base id for this rank (ids become id_base + i)
 */
void initializeParticles_rank(std::vector<Particle>& particles,
                              const SimulationBox& box,
                              const SimulationBox::Decomposition& decomp,
                              int rank,
                              int N_local,
                              bool random,
                              unsigned seed = 0,
                              long long id_base = 0LL);

/**
 * @brief Convenience: initialize from a global N, splitting deterministically across ranks.
 */
void initializeParticles_globalN(std::vector<Particle>& particles,
                                 const SimulationBox& box,
                                 const SimulationBox::Decomposition& decomp,
                                 int rank,
                                 int N_global,
                                 bool random,
                                 unsigned seed = 0,
                                 long long id_stride = 1000000LL);

#endif // INITIAL_MPI_H

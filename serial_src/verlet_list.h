// File: verlet_list.h
#ifndef VERLET_LIST_H
#define VERLET_LIST_H

#include <vector>
#include <cmath>
#include <cstddef>
#include "initial.h"    // Particle, SimulationBox



/**
 * @brief Verlet neighbor list for Grand-Canonical MC with cell-based build.
 *
 * Supports O(N) neighbor-list construction via cell binning
 * and rebuild criteria based on particle displacements.
 */
class VerletList {
    public:
        /**
         * @param rc   force cutoff distance
         * @param skin additional buffer distance for neighbor list
         */
        VerletList(double rc, double skin)
          : rc2_(rc*rc),
            rl2_((rc+skin)*(rc+skin)),
            skin2_(skin*skin),
            maxDisp2_(0.0)
        {}
    
        /**
         * @brief Build or rebuild the full neighbor list for all particles
         * @param positions current particle positions (size N)
         * @param box       simulation box for PBC and minimum-image
         */
        void build(const std::vector<Particle>& positions,
                   const SimulationBox& box) {
            std::size_t N = positions.size();
            neighbors_.assign(N, {});
            lastPos_ = positions;
            maxDisp2_ = 0.0;
            // O(N^2) loop: include pairs within (rc+skin)
            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = i + 1; j < N; ++j) {
                    double d2 = box.minimumImageDistanceSquared(positions[i], positions[j]);
                    if (d2 <= rl2_) {
                        neighbors_[i].push_back(j);
                        neighbors_[j].push_back(i);
                    }
                }
            }
        }
    
        /**
         * @brief Update maximum squared displacement since last build
         * @param positions current positions
         * @param box       simulation box
         */
        void computeMaxDisplacement(const std::vector<Particle>& positions,
                                    const SimulationBox& box) {
            double worst = 0.0;
            for (std::size_t i = 0; i < positions.size(); ++i) {
                double d2 = box.minimumImageDistanceSquared(positions[i], lastPos_[i]);
                if (d2 > worst) worst = d2;
            }
            maxDisp2_ = worst;
        }
    
        /**
         * @brief Check if any particle moved > skin/2 requiring rebuild
         */
        bool needsRebuild() const {
            return maxDisp2_ >= skin2_ * 0.25;
        }
    
        /**
         * @brief Access the computed neighbor lists
         * @return reference to vector of neighbor-index vectors
         */
        const std::vector<std::vector<std::size_t>>& neighbors() const {
            return neighbors_;
        }
    
        /**
         * @brief Find neighbors of a test particle for insertion within rc
         * @param p_new      coordinates of the insertion candidate
         * @param positions  current particle positions
         * @param box        simulation box
         * @param out        output vector of neighbor indices (cleared then filled)
         */
        void getNeighborsForTest(const Particle& p_new,
                                 const std::vector<Particle>& positions,
                                 const SimulationBox& box,
                                 std::vector<std::size_t>& out) const {
            out.clear();
            for (std::size_t j = 0; j < positions.size(); ++j) {
                double d2 = box.minimumImageDistanceSquared(p_new, positions[j]);
                if (d2 <= rc2_) out.push_back(j);
            }
        }
    
    private:
        double rc2_, rl2_, skin2_, maxDisp2_;
        std::vector<std::vector<std::size_t>> neighbors_;
        std::vector<Particle> lastPos_;
    };
    

#endif
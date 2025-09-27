#ifndef CELL_LIST_PARALLEL_H
#define CELL_LIST_PARALLEL_H
#include <utility>
#include "particle.h"
#include "simulation_box.h"
#include <vector>
#include <cstdint>
#include <random>

// 4-color checkerboard on the *interior* cell grid.
// Interior cells are indexed [1..nx_in_] x [1..ny_in_].
// Ghost ring is the frame cells with i==0 or nx_in_+1, j==0 or ny_in_+1.
enum class Parity : std::uint8_t { EvenEven=0, EvenOdd=1, OddEven=2, OddOdd=3 };

class CellListParallel {
public:
    // Static/geometry params provided by caller.
    struct Params {
        double rcut{3.0};  // pair cutoff; also used to size cells (dx,dy >= rcut)
    };

    // ---- Construction / geometry ------------------------------------------------

    // Build for a given global box + decomposition and the rank's id.
    explicit CellListParallel(const SimulationBox& box,
                              const SimulationBox::Decomposition& decomp,
                              int rank,
                              const Params& p);

    // Recompute local geometry (e.g., after Px,Py or domain changes).
    void resetGeometry(const SimulationBox& box,
                       const SimulationBox::Decomposition& decomp,
                       int rank);

    // ---- Binning: owned particles (interior grid only) --------------------------

    // Fill interior bins from current owned list (clears owned bins).
    void buildInterior(const std::vector<Particle>& owned);

    // Update bins after an accepted move for particle index pidx
    // (indices refer to 'owned' vector provided to buildInterior()).
    void onAcceptedMove(int pidx, const Particle& oldp, const Particle& newp);

    // Return interior cells matching a parity, as (ix,iy) with 1<=ix<=nx_in_, 1<=iy<=ny_in_.
    std::vector<std::pair<int,int>> cellsWithParity(Parity par) const;

    // Pick a random owned particle index from cell (ix,iy). Returns -1 if empty.
    int randomOwnedInCell(int ix, int iy, std::mt19937_64& rng) const;

    // ---- Ghost layer (one-cell ring) --------------------------------------------

    // Drop all ghosts and clear ghost bins.
    void clearGhosts();

    // Append one ghost particle (caller ensures it’s relevant for the halo).
    void addGhost(const Particle& g);

    // Re-bin all ghosts into the ghost ring buckets (clears ghost bins).
    void rebuildGhostBins();

    // ---- Neighbor queries -------------------------------------------------------

    // Return neighbors of an owned particle (local index), within rcut.
    // Owned neighbors return (index >= 0), ghost neighbors return (index < 0)
    // encoded as  -1 - ghostIndex  (i.e., ghostIndex = -1 - returnedIndex).
    // The second value is squared distance r^2 (minimum-image).
    std::vector<std::pair<int,double>>
    neighborsOfOwned(int localOwnedIdx, const std::vector<Particle>& owned) const;

    // Same as above, but for an arbitrary global point (xg,yg).
    std::vector<std::pair<int,double>>
    neighborsOfPoint(double xg, double yg, const std::vector<Particle>& owned) const;

    // Map a global point to a local cell index (including the ghost ring).
    // Returns false if the point falls outside the [0..nx_tot_-1]x[0..ny_tot_-1] frame
    // after wrapping (should be rare; guards extreme FP).
    bool mapToLocalCell(double xg, double yg, int& ix, int& iy) const;

    // ---- Accessors (useful for wiring/tests) ------------------------------------

    // Interior grid (must be even in both directions).
    int nxInterior() const { return nx_in_; }
    int nyInterior() const { return ny_in_; }

    // Total grid includes the 1-cell ghost ring on each side.
    int nxTotal() const { return nx_tot_; }
    int nyTotal() const { return ny_tot_; }

    // Cell sizes and local owner bounds.
    double dxCell() const { return dx_cell_; }
    double dyCell() const { return dy_cell_; }
    double x0() const { return x0_; }
    double x1() const { return x1_; }
    double y0() const { return y0_; }
    double y1() const { return y1_; }

    // Global box and rank grid info.
    double Lx() const { return Lx_; }
    double Ly() const { return Ly_; }
    int Px() const { return Px_; }
    int Py() const { return Py_; }
    int rx() const { return cx_; }
    int ry() const { return cy_; }

    double rcut()   const { return rcut_;   }
    double rcutsq() const { return rcutsq_; }

private:
    // ---- geometry / parameters ----
    double rcut_{0.0};
    double rcutsq_{0.0};

    // Global box size (periodic).
    double Lx_{0.0}, Ly_{0.0};

    // Rank grid Px x Py and my coordinates (cx_, cy_).
    int Px_{1}, Py_{1};
    int cx_{0}, cy_{0};

    // Local owner bounds [x0_,x1_) x [y0_,y1_) for this rank.
    double x0_{0.0}, x1_{0.0}, y0_{0.0}, y1_{0.0};

    // Interior cell grid (even counts) and total grid including ghost ring.
    int nx_in_{0}, ny_in_{0}; // interior counts (even)
    int nx_tot_{0}, ny_tot_{0}; // = nx_in_ + 2, ny_in_ + 2
    double dx_cell_{0.0}, dy_cell_{0.0}; // interior cell sizes

    // Bins:
    //  - bins_owned_  : size nx_tot_ * ny_tot_ ; we only use interior indices [1..nx_in_]x[1..ny_in_]
    //  - bins_ghost_  : size nx_tot_ * ny_tot_ ; we only use ring cells (i==0 or nx_in_+1, j==0 or ny_in_+1)
    // Each bucket stores indices into 'owned' (for owned bins) or into ghosts_ (for ghost bins).
    std::vector<std::vector<int>> bins_owned_;
    std::vector<std::vector<int>> bins_ghost_;
    std::vector<Particle> ghosts_; // flat storage for ghost coordinates

    // ---- helpers (internal) ----

    // Compute local owner bounds.
    void computeLocalBounds_(const SimulationBox& box,
                             const SimulationBox::Decomposition& d,
                             int rank);

    // Choose even interior cell counts so that dx,dy >= rcut (robust lower bound).
    void chooseEvenCells_();

    // Clear/recreate bin arrays to current sizes (empties buckets).
    void clearBins_();

    // Convert (ix,iy) in total grid to linear index.
    inline int cid(int ix, int iy) const { return iy * nx_tot_ + ix; }

    // Map global position (xg,yg) to local total-grid cell (including ring).
    // Result ix,iy in [0..nx_tot_-1],[0..ny_tot_-1].
    void localCellOfGlobal_(double xg, double yg, int& ix, int& iy) const;

    // Minimum-image displacement in the *global* periodic box.
    void min_image(double dx, double dy, double& rdx, double& rdy) const;
};

#endif // CELL_LIST_PARALLEL_H

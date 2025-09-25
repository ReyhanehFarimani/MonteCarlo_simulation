#ifndef CELL_LIST_PARALLEL_H
#define CELL_LIST_PARALLEL_H

#include <vector>
#include <utility>
#include "particle.h"
#include "simulation_box.h"

/**
 * @brief Per-rank 2D cell list with a 1-cell ghost layer on all sides.
 *        Local cell ids start at 0 (row-major): cid = ix + iy*nx_tot.
 *
 * Interior cells: ix=1..nx_in, iy=1..ny_in
 * Ghost ring:     ix=0 and ix=nx_in+1; iy=0 and iy=ny_in+1
 *
 * Positions are GLOBAL; mapping uses the rank's subdomain [x0,x1)×[y0,y1).
 */
class CellListParallel {
public:
    struct Params {
        double rcut;             // cutoff
        int    enforce_even = 1; // force even interior cell counts
    };

    CellListParallel(const SimulationBox& box,
                     const SimulationBox::Decomposition& decomp,
                     int rank,
                     const Params& p);

    // Recompute grid if box/decomp changes (does NOT rebuild content)
    void resetGeometry(const SimulationBox& box,
                       const SimulationBox::Decomposition& decomp,
                       int rank);

    // Build interior cells from owned particles (GLOBAL coords)
    void buildInterior(const std::vector<Particle>& owned);

    // Ghost management (you fill these after your halo exchange)
    void clearGhosts();
    void addGhost(const Particle& g);                  // push a ghost particle (GLOBAL coords)
    void rebuildGhostBins();                           // (re)bin current ghosts

    // Neighbor queries (returns local indices of neighbors and squared distance)
    // These search both interior + ghost bins.
    std::vector<std::pair<int,double>>
    neighborsOfOwned(int localOwnedIdx,
                     const std::vector<Particle>& owned) const;

    std::vector<std::pair<int,double>>
    neighborsOfPoint(double xg, double yg,
                     const std::vector<Particle>& owned) const;

    // Accessors
    int nx_in()   const { return nx_in_; }
    int ny_in()   const { return ny_in_; }
    int nx_tot()  const { return nx_tot_; }
    int ny_tot()  const { return ny_tot_; }
    double dx()   const { return dx_cell_; }
    double dy()   const { return dy_cell_; }
    double x0()   const { return x0_; }
    double x1()   const { return x1_; }
    double y0()   const { return y0_; }
    double y1()   const { return y1_; }

    // Map GLOBAL (x,y) to local cell (ix,iy) with ghosts. Returns false if numerically outside by >1 cell.
    bool mapToLocalCell(double xg, double yg, int& ix, int& iy) const;

    // Flat local cell id (row-major)
    inline int cid(int ix, int iy) const { return ix + iy * nx_tot_; }

    // Interior predicate (ix in 1..nx_in, iy in 1..ny_in)
    static inline bool isInterior(int ix, int iy, int nx_in, int ny_in) {
        return (ix >= 1 && ix <= nx_in && iy >= 1 && iy <= ny_in);
    }

private:
    // geometry
    double Lx_, Ly_;
    double x0_, x1_, y0_, y1_;   // subdomain bounds for this rank
    int Px_, Py_;
    int cx_, cy_;                // this rank's coords in the grid
    double dx_cell_, dy_cell_;   // cell size
    int nx_in_, ny_in_;          // interior cells (even)
    int nx_tot_, ny_tot_;        // with ghosts = nx_in+2, ny_in+2
    double rcut_, rcutsq_;

    // bins: owned bins + ghost bins share the same grid; we keep separate index lists
    std::vector<std::vector<int>> bins_owned_; // local owned indices per cell
    std::vector<std::vector<int>> bins_ghost_; // ghost array indices per cell

    // ghost storage (global coords)
    std::vector<Particle> ghosts_;

    // helpers
    void computeLocalBounds_(const SimulationBox& box,
                             const SimulationBox::Decomposition& d,
                             int rank);
    void chooseEvenCells_(); // set nx_in_, ny_in_ and nx_tot_, ny_tot_, dx_cell_, dy_cell_
    void clearBins_();

    inline void min_image(double dx, double dy, double& rdx, double& rdy) const;
    inline void wrapIntoSubdomain_(double& x) const; // wrap global → subdomain periodic image (x)
    inline void wrapIntoSubdomainY_(double& y) const; // (y)

    // Cell index from GLOBAL (x,y); returns (ix,iy) possibly in ghost ring.
    void localCellOfGlobal_(double xg, double yg, int& ix, int& iy) const;
};

#endif // CELL_LIST_PARALLEL_H

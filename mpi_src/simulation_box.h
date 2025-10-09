#ifndef SIMULATION_BOX_H
#define SIMULATION_BOX_H

#include <utility>
#include "particle.h"

class SimulationBox {
private:
    double Lx_, Ly_, invLx_, invLy_, V_;

public:
    explicit SimulationBox(double Lx, double Ly);

    // PBC
    void applyPBC(Particle &p) const;
    void applyPBC(double &x, double &y) const;

    // Minimum-image
    double minimumImageDistanceSquared(const Particle &p1, const Particle &p2) const;
    double minimumImageDistance(const Particle &p1, const Particle &p2) const;

    // Getters
    double getLx() const { return Lx_; }
    double getLy() const { return Ly_; }
    double getV()  const { return V_;  }

    // Setters (keep invariants)
    void setLx(double lx);
    void setLy(double ly);
    void setV(double v); // preserves aspect ratio

    // Recenter after box resize from (lx_old, ly_old)
    void recenter(Particle &p, double lx_old, double ly_old) const;

    // Decomposition (MPI-agnostic)
    struct Decomposition {
        int Px{1}, Py{1};
        double dx{0.0}, dy{0.0};

        int rankOf(int cx, int cy) const;
        std::pair<int,int> coordsOf(int rank) const;

        int west (int rank) const;
        int east (int rank) const;
        int south(int rank) const;
        int north(int rank) const;

        void localBounds(int rank, double Lx, double Ly,
                         double &x0, double &x1, double &y0, double &y1) const;
    };

    Decomposition bestDecomposition(int nranks) const;

    // Constrained version (skip layouts with dx<min_dx or dy<min_dy)
    Decomposition bestDecompositionConstrained(int nranks,
                                               double min_dx, double min_dy,
                                               bool strict = true) const;
};

#endif // SIMULATION_BOX_H

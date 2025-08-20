#ifndef GIBBSMC_H
#define GIBBSMC_H

#include <vector>
#include "initial.h"
#include "cell_list.h"
#include "thermodynamic_calculator.h"
#include "logging.h"
#include "rng.h"
#include "MC.h"

/**
 * @brief Gibbs Monte Carlo driver orchestrator: handles moves, cell-list updates,
 *        sampling, and logging. similar to Montecarlo class but with two simulation box!
 */
class GibbsMonteCarlo {
    public:
    /**
     * @param box1 Simulation box (dimensions & PBC)
     * @param box2 Simulation box (dimensions & PBC)
     * @param particles1 Particle positions of part 1
     * @param particles2 Particle positions of part 2
     * @param calc1 Thermodynamic calculator for energy/virial
     * @param calc2 Thermodynamic calculator for energy/virial
     * @param rcut Interaction cutoff distance
     * @param logger1 Logging instance for trajectory & data of particles 1
     * @param logger2 Logging instance for trajectory & data of particles 2
     * @param rng Seeded random number generator
     */
    GibbsMonteCarlo(SimulationBox& box1,
                    SimulationBox& box2,
                    std::vector<Particle> particles1,
                    std::vector<Particle> particles2,
                    ThermodynamicCalculator& calc1,
                    ThermodynamicCalculator& calc2,
                    double rcut,
                    double delta,
                    double delta_V,
                    Logging& logger1,
                    Logging& logger2,
                    RNG& rng);
    
    /**
     * @brief Run MC simulation for a number of cycles.
     * @param nSteps Number of MC cycles (each cycle = N displacements + possible GC move)
     * @param fOutputStep Frequecy of writing outputs
     * @param fUpdateCell Frequency of updating cell list
     */
    int run(size_t nSteps, size_t fOutputStep, size_t fUpdateCell);

    private:
    SimulationBox box_1;
    SimulationBox box_2;
    std::vector<Particle> particles_1;
    std::vector<Particle> particles_2;
    ThermodynamicCalculator& calc_1;
    ThermodynamicCalculator& calc_2;
    CellList cellList_1;
    CellList cellList_2;
    Logging& logger_1;
    Logging& logger_2;
    RNG& rng_;
    double rcut_;    ///< interaction cutoff
    double delta_;   ///< max displacement magnitude
    double delta_V;
    double beta;
    double energy_1;
    double energy_2;
    int simulation_step_time = 0;

    // Internal moves
    void updateCellList();
    
    bool particle_displacement_1();
    bool particle_displacement_2();
    bool volume_change();
    bool particle_exchange();
    void recordObservables(size_t step);


};
#endif

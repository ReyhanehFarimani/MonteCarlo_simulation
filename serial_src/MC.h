#ifndef MC_H
#define MC_H

#include <vector>
#include "initial.h"
#include "cell_list.h"
#include "thermodynamic_calculator.h"
#include "logging.h"
#include "rng.h"

/**
 * @brief Supported Monte Carlo ensembles.
 */
enum class Ensemble {
    NVT,    ///< Constant NVT ensemble: displacements only
    GCMC,    ///< Grand canonical (μVT): insertion/deletion + displacements
    NPT     ///< isobaric, isothermal ensemble: random walk in lnV + displacement
};

/**
 * @brief Monte Carlo driver orchestrator: handles moves, cell-list updates,
 *        sampling, and logging.
 */
class MonteCarlo {
public:
    /**
     * @param box Simulation box (dimensions & PBC)
     * @param particles Particle positions
     * @param calc Thermodynamic calculator for energy/virial
     * @param rcut Interaction cutoff distance
     * @param ensemble Ensemble type (NVT or GCMC)
     * @param logger Logging instance for trajectory & data
     * @param rng Seeded random number generator
     */
    MonteCarlo(const SimulationBox& box,
               std::vector<Particle>& particles,
               ThermodynamicCalculator& calc,
               double rcut,
               double delta,
               double delta_V,
               Ensemble ensemble,
               Logging& logger,
               RNG& rng);

    /**
     * @brief Run MC simulation for a number of cycles.
     * @param nSteps Number of MC cycles (each cycle = N displacements + possible GC move)
     * @param fOutputStep Frequecy of writing outputs
     * @param fUpdateCell Frequency of updating cell list
     */
    int run(size_t nSteps, size_t fOutputStep, size_t fUpdateCell);
    // void run_no_cell_list(size_t nStep, size_t fOutputStel);

    // --- Testing interfaces ---
    /**
     * @brief Perform one displacement using cell-list ΔE.
     */
    bool stepCellList();
    /**
     * @brief Perform one displacement using full-system energy.
     */
    bool stepBruteForce();

    /**
     * @brief Access current total energy (for testing).
     */
    double getEnergy() const;
    /**
     * @brief Access current particle positions (for testing).
     */
    const std::vector<Particle>& getParticles() const;
    void updateCellList();

private:
    SimulationBox box_;
    std::vector<Particle>& particles_;
    ThermodynamicCalculator& calc_;
    CellList cellList_;
    RNG& rng_;
    double rcut_;    ///< interaction cutoff
    double delta_;   ///< max displacement magnitude
    double delta_V;
    Ensemble ensemble_;
    Logging& logger_;
    double beta;
    double press;
    double energy;
    double z;
    int simulation_step_time = 0;

    // Internal move implementations
    bool displacementMove_cell_list_dE();
    bool displacementMove_no_cell_list();
    bool npt_move_no_cell_list();
    bool grandCanonicalMove();
    void recordObservables(size_t step);
};

#endif // MC_H

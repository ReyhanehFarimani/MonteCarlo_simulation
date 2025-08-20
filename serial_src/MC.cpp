#include "MC.h"
#include <cmath>
#include <iostream>
MonteCarlo::MonteCarlo(const SimulationBox& box,
                       std::vector<Particle>& particles,
                       ThermodynamicCalculator& calc,
                       double rcut,
                       double delta,
                       double delta_V,
                       Ensemble ensemble,
                       Logging& logger,
                       RNG& rng)
    : box_(box),
      particles_(particles),
      calc_(calc),
      cellList_(box, rcut),
      rng_(rng),
      rcut_(rcut),
      delta_(delta),
      delta_V(delta_V),
      ensemble_(ensemble),
      logger_(logger),
      beta(1.0 / calc_.getTemperature()),
      press(calc_.getPressure()),
      z(calc_.getActivity()),
      energy(calc_.computeTotalEnergy(particles, box))
{
    cellList_.build(particles_);
}

int MonteCarlo::run(size_t nSteps, size_t fOutputStep, size_t fUpdateCell) {
        int accept_rate = 0;
        int total = 0;
        for (size_t step = 0; step < nSteps; ++step) {

            size_t N = particles_.size();

            for (size_t i = 0; i < N; ++i) {
                bool accept = displacementMove_cell_list_dE();
                simulation_step_time += 1;
                if (accept) accept_rate++;
                total += 1;
            }
            if (ensemble_ == Ensemble::GCMC) {
                bool accept = grandCanonicalMove();
                simulation_step_time += 1;
                if (accept) accept_rate++;
                total += 1;
                // updateCellList();
            }
            if (ensemble_ == Ensemble::NPT) {
                bool accept = npt_move_no_cell_list();
                simulation_step_time += 1;
                if (accept) updateCellList();
                if (accept) accept_rate++;
                total += 1;
            }
            if (step%fUpdateCell == 0){
                updateCellList();
                double e_tmp = calc_.computeTotalEnergyCellList(particles_, box_);
                if (abs(energy - e_tmp)>10){
                    // std::cout<<"energy_diff = "<<energy - e_tmp<<std::endl;
                    std::cout<<"cell list update frequency is too high!"<<std::endl;
                }
                energy = e_tmp;
                // std::cout<<"step:\t"<<step<<"\t, energy:"<<energy<<std::endl;
            }
            if ((step + 1) % fOutputStep == 0 || step == nSteps - 1)
                recordObservables(simulation_step_time);    
            
    }
    return accept_rate/total;
    // logger_.close();
}

bool MonteCarlo::stepCellList() {
    return displacementMove_cell_list_dE();
}

bool MonteCarlo::stepBruteForce() {
    return displacementMove_no_cell_list();
}

double MonteCarlo::getEnergy() const {
    return energy;
}

const std::vector<Particle>& MonteCarlo::getParticles() const {
    return particles_;
}

bool MonteCarlo::displacementMove_cell_list_dE() {
    size_t idx = rng_.uniformInt(0, static_cast<int>(particles_.size() - 1));
    Particle oldPos = particles_[idx];

    double dx = rng_.uniform(-delta_, delta_);
    double dy = rng_.uniform(-delta_, delta_);
    auto neighbors = cellList_.getNeighbors(idx, particles_);
    double U_old = calc_.computeLocalEnergy(idx, particles_, box_, neighbors);

    particles_[idx].updatePosition(dx, dy);
    box_.applyPBC(particles_[idx]);
    neighbors = cellList_.getNeighbors(idx, particles_);
    double U_new = calc_.computeLocalEnergy(idx, particles_, box_, neighbors);
    double dU = U_new - U_old;
    energy += dU;
    bool accept = false;
    if (dU <= 0.0) {
        accept = true;
    } else {
        double p = std::exp(-beta * dU);
        accept = (rng_.uniform01() < p);
    }
    if (!accept) {
        particles_[idx] = oldPos;
        energy -= dU;
    }
    return accept;
}

bool MonteCarlo::displacementMove_no_cell_list() {
    size_t idx = rng_.uniformInt(0, static_cast<int>(particles_.size() - 1));
    Particle oldPos = particles_[idx];

    double dx = rng_.uniform(-delta_, delta_);
    double dy = rng_.uniform(-delta_, delta_);
    double U_old = calc_.computeTotalEnergy(particles_, box_);

    particles_[idx].updatePosition(dx, dy);
    box_.applyPBC(particles_[idx]);

    double U_new = calc_.computeTotalEnergy(particles_, box_);
    double dU = U_new - U_old;
    energy += dU;
    bool accept = false;
    if (dU <= 0.0) {
        accept = true;
    } else {
        double p = std::exp(-beta * dU);
        accept = (rng_.uniform01() < p);
    }
    if (!accept) {
        particles_[idx] = oldPos;
        energy -= dU;
    }
    return accept;
}

bool MonteCarlo::npt_move_no_cell_list(){
    bool accept = true;
    double V_o = box_.getV();
    double lx_o = box_.getLx();
    double ly_o = box_.getLy();
    double U_old = calc_.computeTotalEnergy(particles_, box_);
    double lnvn = log(V_o) + rng_.uniform(-0.5, 0.5) * delta_V;
    double V_n = exp(lnvn);
    box_.setV(V_n);
    size_t N = particles_.size();
    for (size_t i = 0; i < N; ++i) {
        box_.recenter(particles_[i], lx_o, ly_o);
    }
    double U_new = calc_.computeTotalEnergy(particles_, box_);

    double arg = - beta * ((U_new - U_old) + press * (V_n - V_o) + (-N + 1) * log(V_n/V_o)/beta);
    if (rng_.uniform01() >= exp(arg)) {
        double lx_o = box_.getLx();
        double ly_o = box_.getLy();
        box_.setV(V_o);
        for (size_t i = 0; i < N; ++i) {
            box_.recenter(particles_[i], lx_o, ly_o);
        }
        accept = false;

    }
    return accept;
}
bool MonteCarlo::grandCanonicalMove() {
    double V = box_.getV();
    size_t N = particles_.size();
    bool accept = false;
//  50/50 chance to insert or delete
    double check = rng_.uniform01();
    if (check < 0.3) {
        simulation_step_time += 1;
        // Insertion attempt
        Particle pnew(rng_.uniform(0, box_.getLx()), rng_.uniform(0, box_.getLy()));
        auto new_neighbors = cellList_.getNeighbors2(pnew, particles_);  // neighbors for the new particle
        double dU = calc_.computeLocalEnergy(0, particles_, box_, new_neighbors);
        double acc =  (z * V) / (N + 1) * std::exp(- beta * dU);
        if (acc >= 1){
            energy += dU;
            // cellList_.addParticle(pnew, N);
            particles_.push_back(pnew);
            updateCellList();
            return true;
        }
        else if (rng_.uniform01() < acc) {
            energy += dU;
            // cellList_.addParticle(pnew, N);
            particles_.push_back(pnew);
            updateCellList();
            return true;
        } 
        else {
            return false;
        }
    } else if (check <0.6){
        simulation_step_time += 1;
        // Deletion attempt
        if (N == 0) return false;
        size_t idx = rng_.uniformInt(0, N - 1);
        auto neighbors = cellList_.getNeighbors(idx, particles_);
        double dU = -calc_.computeLocalEnergy(idx, particles_, box_, neighbors);
        double acc =  N / (z * V) * std::exp(- beta * dU);
        if (acc >= 1)
        {
            // cellList_.removeParticle(idx);
            particles_.erase(particles_.begin() + idx);
            updateCellList();
            energy += dU;
            return true;
        }
        else if (rng_.uniform01() < acc) {
            // cellList_.removeParticle(idx);
            particles_.erase(particles_.begin() + idx);
            updateCellList();
            energy += dU;
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }

    return false;
}

void MonteCarlo::updateCellList() {
    cellList_.build(particles_);
}

void MonteCarlo::recordObservables(size_t step) {
    logger_.logSimulationData(particles_, box_, calc_, static_cast<int>(step));
    logger_.logPositions_dump(particles_, box_, step);
}
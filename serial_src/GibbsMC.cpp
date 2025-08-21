#include "GibbsMC.h"
#include <cmath>
#include <iostream>

GibbsMonteCarlo::GibbsMonteCarlo(SimulationBox& box1,
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
                    RNG& rng)
                    : box_1(box1),
                    box_2(box2),
                    particles_1(particles1),
                    particles_2(particles2),
                    calc_1(calc1),
                    calc_2(calc2),
                    cellList_1(box1, rcut),
                    cellList_2(box2, rcut),
                    logger_1(logger1),
                    logger_2(logger2),
                    rng_(rng),
                    rcut_(rcut),
                    delta_(delta),
                    delta_V(delta_V),
                    beta(1./ calc_1.getTemperature())
                    
                    {
                        cellList_1.build(particles_1);
                        cellList_2.build(particles_2);
                        energy_1 = calc1.computeTotalEnergyCellList(particles_1, box_1);
                        energy_2 = calc1.computeTotalEnergyCellList(particles_2, box_2);
                    }

int GibbsMonteCarlo::run(size_t nSteps, size_t fOutputStep, size_t fUpdateCell){
    int accept_rate = 0;
    int total = 0;
    bool accept;
    for (size_t step = 0; step < nSteps; ++step) {

            size_t N1 = particles_1.size();
            for (size_t i = 0; i < N1; ++i) {
                accept = particle_displacement_1();
                if (accept)
                    accept_rate ++;
                total ++;
                simulation_step_time ++;
            }

            size_t N2 = particles_2.size();
            for (size_t i = 0; i < N2; ++i) {
                accept = particle_displacement_2();
                if (accept)
                    accept_rate ++;
                total ++;
                simulation_step_time ++;
            }

            accept = particle_exchange();
            if (accept)
                    accept_rate ++;
            total ++;
            simulation_step_time ++;

            accept = volume_change();
            if (accept)
                    accept_rate ++;
            total ++;
            simulation_step_time ++;

            if (step%fUpdateCell == 0){
                updateCellList();
                double e_tmp1 = calc_1.computeTotalEnergyCellList(particles_1, box_1);
                if (abs(energy_1 - e_tmp1)>10){
                    // std::cout<<"energy_diff = "<<energy - e_tmp<<std::endl;
                    std::cout<<"cell list update frequency is too high!"<<std::endl;
                }
                energy_1 = e_tmp1;
                double e_tmp2 = calc_2.computeTotalEnergyCellList(particles_2, box_2);
                if (abs(energy_2 - e_tmp2)>10){
                    // std::cout<<"energy_diff = "<<energy - e_tmp<<std::endl;
                    std::cout<<"cell list update frequency is too high!"<<std::endl;
                }
                energy_2 = e_tmp2;

                // std::cout<<"step:\t"<<step<<"\t, energy:"<<energy<<std::endl;
            }
            if ((step + 1) % fOutputStep == 0 || step == nSteps - 1)
                recordObservables(simulation_step_time);    
        
        }
    return accept_rate/total;
}
bool GibbsMonteCarlo::particle_displacement_1(){
    bool accept=false;

    size_t idx = rng_.uniformInt(0, static_cast<int>(particles_1.size() - 1));
    Particle oldPos = particles_1[idx];

    double dx = rng_.uniform(-delta_, delta_);
    double dy = rng_.uniform(-delta_, delta_);

    auto neighbors = cellList_1.getNeighbors(idx, particles_1);
    double U_old = calc_1.computeLocalEnergy(idx, particles_1, box_1, neighbors);
    
    particles_1[idx].updatePosition(dx, dy);
    box_1.applyPBC(particles_1[idx]);
    neighbors = cellList_1.getNeighbors(idx, particles_1);

    double U_new = calc_1.computeLocalEnergy(idx, particles_1, box_1, neighbors);
    double dU = U_new - U_old;
    energy_1 += dU;

    if (dU <= 0.0) {
        accept = true;
    } else {
        double p = std::exp(-beta * dU);
        accept = (rng_.uniform01() < p);
    }
    if (!accept) {
        particles_1[idx] = oldPos;
        energy_1 -= dU;
    }
    return accept;

}

bool GibbsMonteCarlo::particle_displacement_2(){
    bool accept=false;

    size_t idx = rng_.uniformInt(0, static_cast<int>(particles_2.size() - 1));
    Particle oldPos = particles_2[idx];

    double dx = rng_.uniform(-delta_, delta_);
    double dy = rng_.uniform(-delta_, delta_);

    auto neighbors = cellList_2.getNeighbors(idx, particles_2);
    double U_old = calc_2.computeLocalEnergy(idx, particles_2, box_2, neighbors);
    
    particles_2[idx].updatePosition(dx, dy);
    box_2.applyPBC(particles_2[idx]);
    neighbors = cellList_2.getNeighbors(idx, particles_2);

    double U_new = calc_1.computeLocalEnergy(idx, particles_2, box_2, neighbors);
    double dU = U_new - U_old;
    energy_2 += dU;
    
    if (dU <= 0.0) {
        accept = true;
    } else {
        double p = std::exp(-beta * dU);
        accept = (rng_.uniform01() < p);
    }
    if (!accept) {
        particles_2[idx] = oldPos;
        energy_2 -= dU;
    }
    return accept;
}


bool GibbsMonteCarlo::particle_exchange(){
    bool accept=false;

    return accept;
}

bool GibbsMonteCarlo::volume_change(){
    bool accept=true;
    double V_o_1 = box_1.getV();
    double lx_o_1 = box_1.getLx();
    double ly_o_1 = box_1.getLy();
    double V_o_2 = box_2.getV();
    double lx_o_2 = box_2.getLx();
    double ly_o_2 = box_2.getLy();
    double u1 = energy_1;
    double u2 = energy_2;
    double U_old = energy_1 + energy_2;
    double dlnV = rng_.uniform(-0.5, 0.5) * delta_V;
    double vn_1_v_2 = std::exp(std::log(V_o_1/V_o_2) + dlnV);
    double V_n_1 = (V_o_1 + V_o_2) * vn_1_v_2 / (1 + vn_1_v_2);
    double V_n_2 = (V_o_1 + V_o_2) - V_n_1;
    box_1.setV(V_n_1);
    box_2.setV(V_n_2);
    size_t N_1 = particles_1.size();
    for (size_t i = 0; i < N_1; ++i) {
        box_1.recenter(particles_1[i], lx_o_1, ly_o_1);
    }
    size_t N_2 = particles_2.size();
    for (size_t i = 0; i < N_2; ++i) {
        box_2.recenter(particles_2[i], lx_o_2, ly_o_2);
    }
    cellList_1.adjust_box(box_1);
    cellList_1.build(particles_1);
    cellList_2.adjust_box(box_2);
    cellList_2.build(particles_2);

    energy_1 = calc_1.computeTotalEnergyCellList(particles_1, box_1);
    energy_2 = calc_2.computeTotalEnergyCellList(particles_2, box_2);
    double U_new = energy_1 + energy_2;
    double arg = -beta * (U_new - U_old);
    arg = std::exp(arg) * std::pow(V_n_1/V_o_1, N_1+1) * std::pow(V_n_2/V_o_2, N_2+1);
    if (rng_.uniform01() >= arg){
        // std::cout<<1<<std::endl;
        // double V_o_1 = box_1.getV();
        double lx_n_1 = box_1.getLx();
        double ly_n_1 = box_1.getLy();
        // double V_o_2 = box_2.getV();
        double lx_n_2 = box_2.getLx();
        double ly_n_2 = box_2.getLy();
        box_1.setV(V_o_1);
        box_2.setV(V_o_2);
        size_t N_1 = particles_1.size();
        for (size_t i = 0; i < N_1; ++i) {
            box_1.recenter(particles_1[i], lx_n_1, ly_n_1);
        }
        size_t N_2 = particles_2.size();
        for (size_t i = 0; i < N_2; ++i) {
            box_2.recenter(particles_2[i], lx_n_2, ly_n_2);
        }
        energy_1 = u1;
        energy_2 = u2;
        cellList_1.adjust_box(box_1);
        cellList_1.build(particles_1);
        cellList_2.adjust_box(box_2);
        cellList_2.build(particles_2);
        accept = false;
    }

    
    return accept;
}

void GibbsMonteCarlo::updateCellList(){
    cellList_1.build(particles_1);
    cellList_2.build(particles_2);
}

void GibbsMonteCarlo::recordObservables(size_t step){
    logger_1.logSimulationData(particles_1, box_1, calc_1, static_cast<int>(step));
    logger_1.logPositions_dump(particles_1, box_1, step);

    logger_2.logSimulationData(particles_2, box_2, calc_2, static_cast<int>(step));
    logger_2.logPositions_dump(particles_2, box_2, step);
}
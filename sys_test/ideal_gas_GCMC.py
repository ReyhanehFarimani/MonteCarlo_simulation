import numpy as np
import os
import matplotlib.pyplot as plt

def theoretical_avgN(mu, temperature, box_area, lambda_square):
    # Calculate the theoretical average number of particles using the ideal gas equation in 2D
    return box_area * np.exp(mu / temperature) / lambda_square

def read_energy_data(filename):
    timesteps = []
    energies = []

    with open(filename, 'r') as file:  # Replace with your actual filename
        for line in file:
            parts = line.split(',')
            timestep = int(parts[0].split(':')[1].strip())
            energy = float(parts[1].split(':')[1].strip())
            timesteps.append(timestep)
            energies.append(energy)
    return timesteps, energies

def read_pressure_data(filename):
    timesteps = []
    pressures = []

    with open(filename, 'r') as file:  # Replace with your actual filename
        for line in file:
            parts = line.split(',')
            timestep = int(parts[0].split(':')[1].strip())
            pressure = float(parts[3].split(':')[1].strip())
            timesteps.append(timestep)
            pressures.append(pressure)
    return timesteps, pressures

def read_nparticles_data(filename):
    timesteps = []
    num_particles = []

    with open(filename, 'r') as file:  # Replace with your actual filename
        for line in file:
            parts = line.split(',')
            timestep = int(parts[0].split(':')[1].strip())
            number_of_particles = int(parts[4].split(':')[1].strip())
            timesteps.append(timestep)
            num_particles.append(number_of_particles)
    return timesteps, num_particles

def compare_results(simulated_data, mu_value, temperature, box_area, lambda_square):
    # Compare each mu value in the simulation data with the theoretical prediction
    sim_avgN = np.mean(simulated_data)
    theory_avgN = theoretical_avgN(mu_value, temperature, box_area, lambda_square)
    if not np.isclose(sim_avgN, theory_avgN, rtol=1e-1):
        print(f"Discrepancy found at mu = {mu_value}: Simulated N = {sim_avgN}, Theoretical N = {theory_avgN}")
        
    return sim_avgN, theory_avgN

if __name__ == "__main__":
    # Constants
    temperature = 1.0
    box_length_x = 20
    box_length_y = 20
    box_area = box_length_x * box_length_y
    lambda_square = 1.0  # Assuming lambda^2 = 1 for simplicity
    sim_data = []
    theo_data = []
    mu_list = np.arange(-2.0, 2.1, 0.1)

    print(r"Testing number of particles, and $\mu$:")

    for mu in mu_list:
        t_list, n_list = read_nparticles_data("data_" + str(mu) + ".txt")
        tmp1, tmp2 = compare_results(n_list, mu, temperature, box_area, lambda_square)
        sim_data.append(tmp1)
        theo_data.append(tmp2)

    ##### A nice plot for the relation of mu to N:

    plt.figure(figsize=(10, 6))
    plt.plot(mu_list, sim_data, 'bo-', label='Simulated Average N')
    plt.plot(mu_list, theo_data, 'r--', label='Theoretical Average N')
    plt.xlabel(r'Chemical Potential $\mu$', fontsize=14)
    plt.ylabel(r'Average Number of Particles $N$', fontsize=14)
    plt.title('Average Number of Particles vs. Chemical Potential', fontsize=16)
    plt.legend()
    plt.grid(True)
    plt.savefig("ideal_gas_GMCM_mu_N.pdf", dpi = 400)

# gcmc_post_processing/data_loader.py
'''
Reyhaneh 28 Aug 2024
'''
import numpy as np

def load_xyz(filename):
    """
    Loads particle positions from an .xyz file.

    Parameters
    ----------
    filename : str
        The path to the .xyz file.

    Returns
    -------
    np.ndarray
        Array of particle positions.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()
        n_particles = int(lines[0])
        positions = np.array([list(map(float, line.split()[1:4]))
                              for line in lines[2:2 + n_particles]])
    return positions


def load_txt_data(filename, maximum_energy_per_particle = 50):
    """
    Loads and processes simulation data from a .txt file.

    The data includes timesteps, number of particles, pressures, and energies.

    Parameters
    ----------
    filename : str
        The path to the .txt file.

    Returns
    -------
    timesteps : list of int
        List of timesteps from the simulation.
    num_particles : list of int
        List of particle counts at each timestep.
    pressures : list of float
        List of pressure values at each timestep.
    energies : list of float
        List of energy values at each timestep.
    """
    timesteps = []
    num_particles = []
    pressures = []
    energies = []

    with open(filename, 'r') as file:
        for line in file:
            parts = line.split(',')
            timestep = int(parts[0].split(':')[1].strip())
            number_of_particles = int(parts[4].split(':')[1].strip())
            pressure = float(parts[3].split(':')[1].strip())
            energy = float(parts[1].split(':')[1].strip())
            if energy < maximum_energy_per_particle * number_of_particles:  # Data cleaning as per your original script
                timesteps.append(timestep)
                num_particles.append(number_of_particles)
                pressures.append(pressure)
                energies.append(energy)
    
    return timesteps, num_particles, pressures, energies


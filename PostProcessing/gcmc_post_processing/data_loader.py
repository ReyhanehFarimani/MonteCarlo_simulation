
'''
Reyhaneh 28 Aug 2024
'''
import numpy as np
import os
import glob
def get_data_files(pattern, directory="."):
    """
    Retrieve a list of data files that match the given pattern in the specified directory.

    Parameters:
    - pattern: The file name pattern to search for (e.g., "NEW_LJ_data_*_*.txt").
    - directory: The directory to search in (defaults to the current directory).

    Returns:
    - List of file paths that match the pattern.
    """
    search_pattern = os.path.join(directory, pattern)
    return glob.glob(search_pattern)

def load_xyz(filename):
    """
    Loads particle positions from an .xyz file, handling variable particle counts per timestep.

    Parameters
    ----------
    filename : str
        The path to the .xyz file.

    Returns
    -------
    timesteps_positions : list of np.ndarray
        A list where each element is an array of particle positions for a timestep.
    """
    timesteps_positions = []
    
    with open(filename, 'r') as file:
        lines = file.readlines()
        i = 0
        while i < len(lines):
            n_particles = int(lines[i].strip())
            positions = []
            i += 1 # Skip the comment line and move to the next timestep
            for j in range(n_particles):
                i += 1
                parts = lines[i].strip().split()
                positions.append([float(parts[1]), float(parts[2]), 0.0])
            timesteps_positions.append(np.array(positions))
            i += 1  # Skip the comment line and move to the next timestep
    
    return timesteps_positions

def load_xyz_crpt_check(filename:str):
    """
    Loads particle positions from an .xyz file, handling variable particle counts per timestep.

    Parameters
    ----------
    filename : str
        The path to the .xyz file.
    n: int 
        Number of particles in case it is const.
    Returns
    -------
    timesteps_positions : list of np.ndarray
        A list where each element is an array of particle positions for a timestep.
    """

    
    timesteps_positions = []
    
    with open(filename, 'r') as file:
        lines = file.readlines()
        if len(lines) == 0:
            print('file empty')
            return 0
        n_particles = int(lines[0].strip())
        if (len(lines)< n_particles + 2):
            print('file empty')
            return 0
        i = 0
        line_counter = len(lines)
        while i < len(lines):
            n_particles = int(lines[i].strip())
            positions = []
            i += 1 # Skip the comment line and move to the next timestep
            line_counter -= 2
            if line_counter < n_particles:
                return timesteps_positions
            for j in range(n_particles):
                i += 1
                parts = lines[i].strip().split()
                try:
                    positions.append([float(parts[1]), float(parts[2]), 0.0])
                except:
                    return timesteps_positions
                line_counter -= 1
            timesteps_positions.append(np.array(positions))
            i += 1  # Skip the comment line and move to the next timestep
    
    return timesteps_positions


def load_txt_data(filename, maximum_energy_per_particle = 10):
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


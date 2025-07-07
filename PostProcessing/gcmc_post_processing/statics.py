import pandas as pd
import re
from gcmc_post_processing.data_loader import load_txt_data  # Assuming this function exists
import numpy as np


def extract_mu_and_temperature(filename):
    """
    Extracts mu and temperature from a filename by looking for floating-point numbers,
    including handling negative values.

    Parameters:
    - filename: The filename to extract from.

    Returns:
    - mu: The first floating-point number found in the filename.
    - temperature: The second floating-point number found in the filename.
    """
    # Regular expression to find all floating-point numbers including negative numbers
    numbers = re.findall(r"-?\d+\.\d+|-?\d+", filename)
    
    if len(numbers) < 2:
        raise ValueError(f"Filename {filename} does not contain enough numeric values to extract mu and temperature.")
    
    # Convert the found numbers to floats and assume the first is mu and the second is temperature
    mu = float(numbers[0])
    temperature = float(numbers[1])
    
    return mu, temperature
def read_simulation_input_Old(input_file):
    """
    Reads the simulation input file and extracts relevant parameters.
    
    Parameters:
    - input_file: Path to the input.txt file.
    
    Returns:
    - params: Dictionary containing simulation parameters.
    """
    params = {}
    
    with open(input_file, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) < 2:
                continue  # Skip lines that don't have key-value pairs
            try:
                if "mu" in line:
                    params['mu'] = float(parts[-1])
                elif "f" in line:
                    params['f'] = float(parts[-1])
                elif "boxLengthX" in line:
                    params['boxLengthX'] = float(parts[-1])
                elif "boxLengthY" in line:
                    params['boxLengthY'] = float(parts[-1])
                elif "temperature" in line:
                    params['T'] = float(parts[-1])
                elif "kappa" in line:
                    params['kappa'] = float(parts[-1])
            except ValueError:
                # print(f"Skipping line due to conversion error: {line}")
                pass
    
    return params
def read_simulation_input(input_file):
    """
    Reads the simulation input file and extracts relevant parameters.
    
    Parameters:
    - input_file: Path to the input.txt file.
    
    Returns:
    - params: Dictionary containing simulation parameters.
    """
    params = {}
    
    with open(input_file, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) < 2:
                continue  # Skip lines that don't have key-value pairs
            try:
                if "mu" in line:
                    params['mu'] = float(parts[-1])
                elif "f" in line:
                    params['f'] = float(parts[-1])
                elif "Lx" in line:
                    params['boxLengthX'] = float(parts[-1])
                elif "Ly" in line:
                    params['boxLengthY'] = float(parts[-1])
                elif "T" in line:
                    params['T'] = float(parts[-1])
                elif "kappa" in line:
                    params['kappa'] = float(parts[-1])
            except ValueError:
                # print(f"Skipping line due to conversion error: {line}")
                pass
    
    return params


def process_simulation_data(data_files, input_files, lag):
    """
    Process simulation data from a list of files, extracting mu and temperature from filenames.

    Parameters:
    - data_files: List of file paths to process.
    - box_area: Area of the simulation box.

    Returns:
    - detailed_df: DataFrame containing detailed data from all simulations.
    - avg_df: DataFrame containing averaged data from all simulations.
    """
    detailed_records = []
    avg_records = []

    for filename, input_file in zip(data_files, input_files):
        # Extract mu and temperature from the filename
        try:
            params = read_simulation_input(input_file)
        except:
            params = read_simulation_input_Old(input_file)
        print(filename)
        f = params['f']
        mu = params['mu']
        l = -(params['kappa'] - 6.56 )/7.71
        if (l==6.56/7.71):
            l = 0
        box_area = params['boxLengthX'] * params['boxLengthY']
        temperature = params['T']
        # Load simulation data
        _, num_particles, pressures, energies = load_txt_data(filename, 1000)
        print(pressures)
        num_particles = num_particles[lag:]
        pressures = pressures[lag:]
        energies = energies[lag:]
        if len(pressures) > 0:
            print('1')
            sim_avgN = np.mean(num_particles)
            avg_pressure = np.mean(pressures)
            avg_energy = np.mean(energies)
            stddevN = np.std(num_particles) / np.sqrt(len(num_particles))
            stddevP = np.std(pressures) / np.sqrt(len(pressures))
            stddevE = np.std(energies) / np.sqrt(len(energies))

            # Record detailed data for each timestep
            for n, p, e in zip(num_particles, pressures, energies):
                detailed_records.append({
                    'mu': mu,
                    'temperature': temperature,
                    'num_particles': n,
                    'density': n / box_area,
                    'pressure': p,
                    'energy': e,
                    'f': f,
                    'l': l,
                    'bx_area' : box_area
                })

            # Record averaged data for the entire file
            avg_records.append({
                'mu': mu,
                'temperature': temperature,
                'l': l,
                'sim_avgN': sim_avgN,
                'avg_pressure': avg_pressure,
                'avg_energy': avg_energy,
                'avg_density': sim_avgN/box_area,
                'stddevN': stddevN,
                'stddevP': stddevP,
                'stddevE': stddevE,
                'stddevrho': stddevN/box_area,
                'f': f,
                'bx_area' : box_area
            })
    
    detailed_df = pd.DataFrame(detailed_records)
    avg_df = pd.DataFrame(avg_records)
    # avg_df.sort_values(by=['mu'], inplace=True)

    return detailed_df, avg_df


def bin_data_by_density(df, density_bins, tolerance=0.01):
    binned_data = []
    for density_bin in density_bins:
        bin_indices = np.abs(df['density'] - density_bin) < tolerance
        if bin_indices.any():
            avg_pressure = df[bin_indices]['pressure'].mean()
            stddev_pressure = df[bin_indices]['pressure'].std() / np.sqrt(bin_indices.sum())
            binned_data.append({
                'density_bin': density_bin,
                'avg_pressure': avg_pressure,
                'stddev_pressure': stddev_pressure,
                'count': bin_indices.sum()
            })
    return pd.DataFrame(binned_data)



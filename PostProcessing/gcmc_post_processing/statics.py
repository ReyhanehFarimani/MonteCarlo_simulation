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

def process_simulation_data(data_files, box_area=1.0):
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

    for filename in data_files:
        # Extract mu and temperature from the filename
        mu, temperature = extract_mu_and_temperature(filename)
        print(mu, temperature)
        # Load simulation data
        _, num_particles, pressures, energies = load_txt_data(filename)
        
        if len(pressures) > 0:
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
                    'energy': e
                })

            # Record averaged data for the entire file
            avg_records.append({
                'mu': mu,
                'temperature': temperature,
                'sim_avgN': sim_avgN,
                'avg_pressure': avg_pressure,
                'avg_energy': avg_energy,
                'avg_density': sim_avgN/box_area,
                'stddevN': stddevN,
                'stddevP': stddevP,
                'stddevE': stddevE,
                'stddevrho': stddevN/box_area
            })
    
    detailed_df = pd.DataFrame(detailed_records)
    avg_df = pd.DataFrame(avg_records)
    avg_df.sort_values(by=['mu'], inplace=True)

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
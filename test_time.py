import numpy as np
import matplotlib.pyplot as plt
import subprocess
import time

# Base configuration for the input file template
base_config = {
    "boxLengthX": 100.0,
    "boxLengthY": 10.0,
    "numParticles": 600,
    "randomPlacement": 1,
    "temperature": 1,
    "numSteps": 1000,
    "timeStep": 0.1,
    "outputFrequency": 100,
    "equilibrationTime": 1,
    "r2cut": 2.5,
    "potentialType": "WCA",
    "simulationType": "GCMC",
    "mu": 500.0,
    "seed": 0,
    "useCellList": 1,
    "cellListUpdateFrequency": 50,
    "positionFile": "particle_positions.xyz",
    "dataFile": "simulation_data.dat"
}

# Define simulation configurations
configurations = [
    {"cpus": 1, "exec_command": "./Monte_carlo_serial"},
    {"cpus": 4, "exec_command": "mpirun -n 4 ./Monte_carlo_MPI"},
    {"cpus": 8, "exec_command": "mpirun -n 8 ./Monte_carlo_MPI"}
]

# Define box size scaling for each case to keep density constant
scaling_factors = [1, 2, 3, 4, 5, 6, ]  # Increase box size to match density across cases

# Store results for plotting
timing_results = []

# Create input files, run simulations, and measure execution time
for config, scale in zip(configurations, scaling_factors):
    # Update box size and particle count for each case
    adjusted_config = base_config.copy()
    adjusted_config["boxLengthX"] = base_config["boxLengthX"] * scale
    adjusted_config["boxLengthY"] = base_config["boxLengthY"] * scale
    adjusted_config["numParticles"] = int(base_config["numParticles"] * (scale ** 2))

    # Write input file
    with open("input.txt", "w") as f:
        for key, value in adjusted_config.items():
            f.write(f"{key} {value}\n")
    
    # Run the simulation and measure the time taken
    exec_command = config["exec_command"]
    print(f"Running: {exec_command} with scale {scale}x")

    start_time = time.time()
    try:
        # Execute the command and wait for it to finish
        subprocess.run(exec_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running {exec_command}: {e}")
        continue  # Skip this configuration if there was an error
    end_time = time.time()

    # Calculate execution time
    execution_time = end_time - start_time
    timing_results.append({"cpus": config["cpus"], "scale": scale, "execution_time": execution_time})

# Extract data for plotting
cpus = [result["cpus"] for result in timing_results]
execution_times = [result["execution_time"] for result in timing_results]

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(cpus, execution_times, marker='o', linestyle='-', color='b')
plt.title("Execution Time vs CPU Count for Constant Density Simulation")
plt.xlabel("Number of CPUs")
plt.ylabel("Execution Time (s)")
plt.xticks(cpus)
plt.grid(True)
plt.show()

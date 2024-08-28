from gcmc_post_processing.data_loader import load_txt_data, load_xyz
import matplotlib.pyplot as plt

def test_load_txt_data():
    # Test loading data from a .txt file
    txt_file = 'simulation_data.txt'  # Adjust the filename as needed
    timesteps, num_particles, pressures, energies = load_txt_data(txt_file)
    
    print("TXT Data:")
    print(f"Timesteps: {timesteps[:5]}")
    print(f"Number of Particles: {num_particles[:5]}")
    print(f"Pressures: {pressures[:5]}")
    print(f"Energies: {energies[:5]}")

def test_load_xyz_and_plot():
    # Test loading positions from an .xyz file and plot the second timestep
    xyz_file = 'particle_positions.xyz' 
    positions = load_xyz(xyz_file)
    
    # Ensure there are at least two timesteps
    if len(positions) < 2:
        print("Not enough timesteps in the file to plot the second one.")
        return

    # Extract the positions for the second timestep
    second_timestep_positions = positions[1]
    
    # Plot the positions using a scatter plot (assume 2D scatter plot using x and y coordinates)
    plt.scatter(second_timestep_positions[:, 0], second_timestep_positions[:, 1], c='blue', marker='o')
    plt.title("Particle Positions at Second Timestep")
    plt.xlabel("X Position")
    plt.ylabel("Y Position")
    plt.show()

if __name__ == "__main__":
    test_load_txt_data()
    test_load_xyz_and_plot()

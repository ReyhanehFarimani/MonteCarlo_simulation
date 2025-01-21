import numpy as np

def generate_fake_crystal_structure(lattice_constant=1.0, grid_size=(10, 10), output_file='fake_crystal_structure.xyz'):
    """
    Generates a fake 2D crystal structure with a simple square lattice.
    
    Parameters
    ----------
    lattice_constant : float
        The lattice constant (distance between adjacent particles).
    grid_size : tuple of int
        The number of particles along each dimension (x, y).
    output_file : str
        The filename for the output .xyz file.
    
    Returns
    -------
    None
    """
    x_particles, y_particles = grid_size
    num_particles = int(x_particles/ lattice_constant) * int(y_particles/ lattice_constant) 
    
    positions = []
    
    for i in range(int(x_particles/ lattice_constant)):
        for j in range(int(y_particles/ lattice_constant)):
            x = i * lattice_constant
            y = j * lattice_constant
            positions.append([x, y, 0.0])  # Adding 0.0 as the z-coordinate for a 2D structure
    
    positions = np.array(positions)
    
    # Write the positions to an .xyz file
    with open(output_file, 'w') as file:
        file.write(f"{num_particles}\n")
        file.write("Fake crystal structure\n")
        for pos in positions:
            file.write(f"C {pos[0]:.5f} {pos[1]:.5f} {pos[2]:.5f}\n")
    
    print(f"Fake crystal structure written to {output_file}")




def generate_fake_triangular_lattice(lattice_constant=1.0, grid_size=(10, 10), output_file='fake_triangular_lattice.xyz'):
    """
    Generates a fake 2D crystal structure with a triangular lattice.

    Parameters
    ----------
    lattice_constant : float
        The lattice constant (distance between adjacent particles).
    grid_size : tuple of int
        The number of particles along each dimension (x, y).
    output_file : str
        The filename for the output .xyz file.

    Returns
    -------
    None
    """
    x_particles, y_particles = grid_size
    num_particles = 0
    
    positions = []

    for i in range(int(x_particles / lattice_constant)):
        for j in range(int(y_particles / lattice_constant)):
            x = i * lattice_constant
            y = j * (lattice_constant * (3**0.5 / 2))  # Vertical spacing for triangular lattice
            if j % 2 == 1:  # Offset every other row
                x += lattice_constant / 2
            positions.append([x, y, 0.0])  # Adding 0.0 as the z-coordinate for a 2D structure
            num_particles += 1

    positions = np.array(positions)

    # Write the positions to an .xyz file
    with open(output_file, 'w') as file:
        file.write(f"{num_particles}\n")
        file.write("Fake triangular lattice\n")
        for pos in positions:
            file.write(f"C {pos[0]:.5f} {pos[1]:.5f} {pos[2]:.5f}\n")

    print(f"Fake triangular lattice written to {output_file}")
    

def generate_fake_triangular_lattice_defected(lattice_constant=1.0, grid_size=(10, 10), output_file='fake_triangular_lattice.xyz', defect = 0.01):
    """
    Generates a fake 2D crystal structure with a triangular lattice.

    Parameters
    ----------
    lattice_constant : float
        The lattice constant (distance between adjacent particles).
    grid_size : tuple of int
        The number of particles along each dimension (x, y).
    output_file : str
        The filename for the output .xyz file.

    Returns
    -------
    None
    """
    x_particles, y_particles = grid_size
    num_particles = 0
    
    positions = []

    for i in range(int(x_particles / lattice_constant)):
        for j in range(int(y_particles / lattice_constant)):
            x_R = (np.random.rand() - 0.5) * defect * lattice_constant
            y_R = (np.random.rand() - 0.5) * defect * lattice_constant
            x = i * lattice_constant
            y = j * (lattice_constant * (3**0.5 / 2))  # Vertical spacing for triangular lattice
            if j % 2 == 1:  # Offset every other row
                x += lattice_constant / 2
            positions.append([x + x_R, y + y_R, 0.0])  # Adding 0.0 as the z-coordinate for a 2D structure
            num_particles += 1

    positions = np.array(positions)

    # Write the positions to an .xyz file
    with open(output_file, 'w') as file:
        file.write(f"{num_particles}\n")
        file.write("Fake triangular lattice\n")
        for pos in positions:
            file.write(f"C {pos[0]:.5f} {pos[1]:.5f} {pos[2]:.5f}\n")

    print(f"Fake triangular lattice written to {output_file}")
    


def generate_random_particles(num_particles=100, box_size=(10.0, 10.0), output_file='random_particles.xyz'):
    """
    Generates a set of randomly distributed 2D particles within a specified box size.
    
    Parameters
    ----------
    num_particles : int
        The number of particles to generate.
    box_size : tuple of float
        The size of the simulation box in the x and y dimensions.
    output_file : str
        The filename for the output .xyz file.
    
    Returns
    -------
    None
    """
    x_size, y_size = box_size
    
    # Generate random positions for the particles
    positions = np.random.rand(num_particles, 2)
    positions[:, 0] *= x_size  # Scale x positions to the box size
    positions[:, 1] *= y_size  # Scale y positions to the box size
    
    # Add a zero z-coordinate to make it compatible with 3D .xyz format
    positions = np.hstack((positions, np.zeros((num_particles, 1))))
    
    # Write the positions to an .xyz file
    with open(output_file, 'w') as file:
        file.write(f"{num_particles}\n")
        file.write("Random particle distribution\n")
        for pos in positions:
            file.write(f"C {pos[0]:.5f} {pos[1]:.5f} {pos[2]:.5f}\n")
    
    print(f"Random particle distribution written to {output_file}")
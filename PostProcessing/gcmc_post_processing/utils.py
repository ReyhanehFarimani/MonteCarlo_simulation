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


def generate_fake_honeycomb_structure(a=1.0, grid_size=(5 * 3 ** 0.5, 10), output_file='fake_honeycomb_structure.xyz'):
    """
    Generates a fake 2D honeycomb structure with a given lattice constant.
    
    Parameters
    ----------
    a : float
        The lattice constant (distance between adjacent particles).
    grid_size : tuple of int
        The number of hexagons along each dimension (x, y).
    output_file : str
        The filename for the output .xyz file.
    
    Returns
    -------
    None
    """
    x_hexagons, y_hexagons = grid_size
    num_particles = int(x_hexagons/ a) * int(y_hexagons/ a) * 2
    
    positions = []

    # The unit vectors for the honeycomb lattice
    a1 = np.array([a, 0])
    a2 = np.array([a/2, np.sqrt(3)*a/2])

    for i in range(int(x_hexagons)):
        for j in range(int(y_hexagons)):
            # Position of the first atom in the unit cell
            r1 = i * a1 + j * a2
            # Position of the second atom in the unit cell
            r2 = r1 + np.array([a/2, np.sqrt(3)*a/2])
            positions.append([r1[0], r1[1], 0.0])
            positions.append([r2[0], r2[1], 0.0])

    positions = np.array(positions)

    # Write the positions to an .xyz file
    with open(output_file, 'w') as file:
        file.write(f"{num_particles}\n")
        file.write("Fake honeycomb structure\n")
        for pos in positions:
            file.write(f"C {pos[0]:.5f} {pos[1]:.5f} {pos[2]:.5f}\n")
    
    print(f"Fake honeycomb structure written to {output_file}")



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
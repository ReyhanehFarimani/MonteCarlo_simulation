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
    num_particles = x_particles * y_particles
    
    positions = []
    
    for i in range(x_particles):
        for j in range(y_particles):
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
    num_particles = x_hexagons * y_hexagons * 2  # Two particles per hexagon

    positions = []

    # The unit vectors for the honeycomb lattice
    a1 = np.array([a, 0])
    a2 = np.array([a/2, np.sqrt(3)*a/2])

    for i in range(x_hexagons):
        for j in range(y_hexagons):
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




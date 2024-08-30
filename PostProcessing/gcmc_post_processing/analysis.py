import numpy as np
import freud

def compute_rdf(particles, box_length, dr, rcutoff=0.9):
    """
    Computes the 2D radial distribution function g(r) of a set of particle coordinates,
    considering periodic boundary conditions using freud.

    Parameters
    ----------
    particles : (N, 2) np.array
        Array of particle positions with shape (N, 2), where N is the number of particles.
    box_length : float or np.ndarray
        Lengths of the simulation box in each dimension. If float, assumed to be the same in both dimensions.
    dr : float
        The width of the radial bins. Determines the spacing between successive radii over which g(r) is computed.
    rcutoff : float, optional
        The cutoff fraction of the maximum radius. Default is 0.9, meaning g(r) will be computed up to 90% of the half-box length.

    Returns
    -------
    rdf : np.ndarray
        The computed g(r) values.
    r_bins : np.ndarray
        The radii corresponding to the g(r) values.
    """
    if isinstance(box_length, float) or isinstance(box_length, int):
        box_size = [box_length, box_length]
    else:
        box_size = box_length

    # Define the box
    box = freud.box.Box(Lx=box_size[0], Ly=box_size[1], Lz=0.0)

    # Calculate the maximum radius
    r_max = (min(box_size) / 2) * rcutoff

    # Initialize the RDF computation
    rdf = freud.density.RDF(bins=int(r_max / dr), r_max=r_max)

    # Compute the RDF
    rdf.compute(system=(box, particles))

    # Return the RDF values and the corresponding radii
    return rdf.rdf, rdf.bin_centers

def average_rdf_over_trajectory(particle_trajectories, box_length, dr, rcutoff=0.9):
    """
    Averages the radial distribution function g(r) over all time steps of a trajectory.

    Parameters
    ----------
    particle_trajectories : list of (N, 2) np.array
        List of particle position arrays with shape (N, 2) for each timestep.
    box_length : float or np.ndarray
        Lengths of the simulation box in each dimension. If float, assumed to be the same in both dimensions.
    dr : float
        The width of the radial bins.
    rcutoff : float, optional
        The cutoff fraction of the maximum radius.
    
    Returns
    -------
    avg_rdf : np.ndarray
        The averaged g(r) values.
    r_bins : np.ndarray
        The radii corresponding to the g(r) values.
    """
    rdf_accumulator = None
    num_timesteps = len(particle_trajectories)

    for timestep in particle_trajectories:
        rdf, r_bins = compute_rdf(timestep, box_length, dr, rcutoff)
        
        if rdf_accumulator is None:
            rdf_accumulator = np.zeros_like(rdf)
        
        rdf_accumulator += rdf

    avg_rdf = rdf_accumulator / num_timesteps

    return avg_rdf, r_bins


def voronoi_neighborhood_list(points, box_size):
    """
    Generates the Voronoi neighborhood list for a set of 2D points with periodic boundary conditions using freud.

    Parameters
    ----------
    points : np.ndarray
        Array of points with shape (N, 2), where N is the number of points.
    box_size : tuple of float
        The size of the box in the x and y dimensions.

    Returns
    -------
    neighbors : dict
        A dictionary where each key is a point index, and the value is a list of neighboring point indices.
    """
    # Define the box and Voronoi analysis
    box = freud.box.Box(Lx=box_size[0], Ly=box_size[1])
    voronoi = freud.locality.Voronoi(box)
    
    # Compute the Voronoi diagram
    voronoi.compute(system=(box, points))
    
    # Extract the Voronoi neighbors
    neighbors = {i: set() for i in range(len(points))}
    for bond in voronoi.nlist:
        neighbors[bond[0]].add(bond[1])
        neighbors[bond[1]].add(bond[0])
    
    # Convert sets to sorted lists for consistency
    for i in neighbors:
        neighbors[i] = sorted(neighbors[i])
    
    return neighbors


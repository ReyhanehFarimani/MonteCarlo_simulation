
'''
Reyhaneh 28 Aug 2024
'''
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def gr_only_centered_particles(particles, box, dr, rho=None, rcutoff=0.8, eps=1e-15):
    """
    Computes the 2D radial distribution function g(r) of a set of particle coordinates.
    It only considers the particles in the center of the simulation box to avoid
    issues related to periodic boundary conditions (PBC).

    Parameters
    ----------
    particles : (N, 2) np.array
        Array of particle positions with shape (N, 2), where N is the number of particles.
    box : (2,) np.array or list
        The dimensions of the simulation box.
    dr : float
        The width of the radial bins. Determines the spacing between successive radii over which g(r) is computed.
    rho : float, optional
        The number density of the system. If None, the density will be calculated from the particle positions and box dimensions.
    rcutoff : float, optional
        The cutoff fraction of the maximum radius. Default is 0.8, meaning g(r) will be computed up to 80% of the half-box length.
    eps : float, optional
        A small value used to avoid edge cases in distance comparisons (default is 1e-15).

    Returns
    -------
    g_r : (n_radii,) np.array
        Radial distribution function values g(r) corresponding to each radius.
    radii : (n_radii,) np.array
        Array of radii at which g(r) is computed.
    """

    # Ensure particles are in the correct format (N, 2)
    if not isinstance(particles, np.ndarray):
        particles = np.array(particles)

    # Translate particles such that the minimum coordinates start at the origin
    particles -= np.array([0, 0])

    # Determine the dimensions of the box (width and height)
    dims = np.array(box)
    
    # Calculate the maximum radius to consider, limited by the rcutoff
    r_max = (np.min(dims) / 2) * rcutoff
    
    # Generate an array of radii where g(r) will be computed
    radii = np.arange(dr, r_max, dr)

    # Determine the number of particles
    N, d = particles.shape
    
    # If density (rho) is not provided, calculate it from the number of particles and the box area
    if rho is None:
        rho = N / np.prod(dims)  # Density in 2D (particles per unit area)

    # Use a KDTree for efficient neighbor searching
    tree = cKDTree(particles)

    # Initialize an array to hold the g(r) values
    g_r = np.zeros(shape=(len(radii)))

    # Loop over each radius in the radii array
    for r_idx, r in enumerate(radii):
        # Find all particles that are at least r + dr away from the edges of the box
        valid_idxs = np.bitwise_and.reduce([
            (particles[:, i] - (r + dr) >= 0) & (particles[:, i] + (r + dr) <= dims[i])
            for i in range(d)
        ])
        valid_particles = particles[valid_idxs]

        # Compute the number of neighbors within the shell for each valid particle
        for particle in valid_particles:
            # Number of particles within the shell between r and r + dr
            n = tree.query_ball_point(particle, r + dr - eps, return_length=True) - \
                tree.query_ball_point(particle, r, return_length=True)
            g_r[r_idx] += n

        # Normalize g(r) for this radius
        n_valid = len(valid_particles)
        if n_valid > 0:  # Ensure there are valid particles before dividing
            shell_area = np.pi * ((r + dr) ** 2 - r ** 2)  # Area of the shell in 2D
            g_r[r_idx] /= n_valid * shell_area * rho  # Normalize by the number of particles and the shell area

    return g_r, radii


def compute_gr_only_centered_particles(particles_timestep, box, dr, rho=None, rcutoff=0.8, eps=1e-15):
    """
    Computes the 2D radial distribution function g(r) over multiple timesteps.
    It only considers the particles in the center of the simulation box to avoid 
    issues related to periodic boundary conditions (PBC).

    Parameters
    ----------
    particles_timestep : list of (N, 2) np.array
        List of arrays containing particle positions for each timestep.
    box : (2,) np.array or list
        The dimensions of the simulation box.
    dr : float
        The width of the radial bins. Determines the spacing between successive radii over which g(r) is computed.
    rho : float, optional
        The number density of the system. If None, the density will be calculated from the particle positions and box dimensions.
    rcutoff : float, optional
        The cutoff fraction of the maximum radius. Default is 0.8, meaning g(r) will be computed up to 80% of the half-box length.
    eps : float, optional
        A small value used to avoid edge cases in distance comparisons (default is 1e-15).

    Returns
    -------
    radii : (n_radii,) np.array
        Array of radii at which g(r) is computed.
    g_r : (n_radii,) np.array
        Averaged radial distribution function values g(r) over all timesteps.
    """

    # Initialize a list to store g(r) for each timestep
    g_r_all = []

    # Loop through each timestep
    for ts in particles_timestep:
        g_r_tmp, radii = gr_only_centered_particles(ts, box, dr, rho, rcutoff, eps)
        g_r_all.append(g_r_tmp)
    
    # Average g(r) over all timesteps
    g_r = np.mean(g_r_all, axis=0)

    return radii, g_r

    
    
    
def minimum_image_distance(particles, box_length):
    """
    Computes the minimum image distance between all pairs of particles in a periodic box.

    Parameters
    ----------
    particles : np.ndarray
        Array of particle positions with shape (N, 2), where N is the number of particles.
    box_length : np.ndarray or float
        Lengths of the simulation box in each dimension. If float, assumed to be the same in both dimensions.

    Returns
    -------
    distances : np.ndarray
        A 2D array of shape (N, N) containing the minimum image distances between all particle pairs.
    """
    # Compute the differences between all pairs of particles
    delta = particles[:, np.newaxis, :] - particles[np.newaxis, :, :]
    
    # Apply minimum image convention
    delta = np.where(delta > 0.5 * box_length, delta - box_length, delta)
    delta = np.where(delta < -0.5 * box_length, delta + box_length, delta)
    
    # Compute the Euclidean distances
    distances = np.sqrt(np.sum(delta**2, axis=-1))
    
    return distances

def gr(particles, box_length, dr, rho=None, rcutoff=0.9, eps=1e-15):
    """
    Computes the 2D radial distribution function g(r) of a set of particle coordinates,
    considering periodic boundary conditions using the minimum image convention.

    Parameters
    ----------
    particles : (N, 2) np.array
        Array of particle positions with shape (N, 2), where N is the number of particles.
    box_length : np.ndarray or float
        Lengths of the simulation box in each dimension. If float, assumed to be the same in both dimensions.
    dr : float
        The width of the radial bins. Determines the spacing between successive radii over which g(r) is computed.
    rho : float, optional
        The number density of the system. If None, the density will be calculated from the particle positions and box dimensions.
    rcutoff : float, optional
        The cutoff fraction of the maximum radius. Default is 0.9, meaning g(r) will be computed up to 90% of the half-box length.
    eps : float, optional
        A small value used to avoid edge cases in distance comparisons (default is 1e-15).

    Returns
    -------
    g_r : (n_radii,) np.array
        Radial distribution function values g(r) corresponding to each radius.
    radii : (n_radii,) np.array
        Array of radii at which g(r) is computed.
    """

    # Ensure particles are in the correct format (N, 2)
    if not isinstance(particles, np.ndarray):
        particles = np.array(particles)

    # Ensure box_length is a NumPy array
    if np.isscalar(box_length):
        box_length = np.array([box_length, box_length])

    # Calculate the maximum radius to consider, limited by the rcutoff
    r_max = (np.min(box_length) / 2) * rcutoff
    
    # Generate an array of radii where g(r) will be computed
    radii = np.arange(dr, r_max, dr)

    # Determine the number of particles and dimensionality (should be 2)
    N, d = particles.shape
    
    assert(d==2)
    
    # If density (rho) is not provided, calculate it from the number of particles and the box area
    if rho is None:
        rho = N / np.prod(box_length)  # Density in 2D (particles per unit area)

    # Calculate all pairwise distances using the minimum image convention
    distances = minimum_image_distance(particles, box_length)
    # Only consider the upper triangle of the distance matrix to avoid double counting
    i_upper = np.triu_indices(N, k=1)
    distances = distances[i_upper]
    
    # Initialize an array to hold the g(r) values
    g_r = np.zeros(shape=(len(radii)))

    # Loop over each radius in the radii array and count distances that fall into each radial bin
    for r_idx, r in enumerate(radii):
        # Count the number of particle pairs with distances within the current shell [r, r+dr)
        in_shell = np.logical_and(distances >= r, distances < r + dr - eps)
        g_r[r_idx] = np.sum(in_shell)

    # Normalize g(r) for each radius
    shell_areas = np.pi * ((radii + dr)**2 - radii**2)  # Area of the shell in 2D
    g_r /= (N * (N - 1) / 2) * shell_areas * rho  # Normalize by the number of pairs and the shell area

    return g_r, radii

def compute_gr(particles_timestep, box_length, dr, rho=None, rcutoff=0.9, eps=1e-15):
    """
    Computes the 2D radial distribution function g(r) of a set of particle coordinates,
    considering periodic boundary conditions using the minimum image convention.

    Parameters
    ----------
    particles_timestep : list of (N, 2) np.array
        List of array of particle positions with shape (N, 2), where N is the number of particles.
    box_length : np.ndarray or float
        Lengths of the simulation box in each dimension. If float, assumed to be the same in both dimensions.
    dr : float
        The width of the radial bins. Determines the spacing between successive radii over which g(r) is computed.
    rho : float, optional
        The number density of the system. If None, the density will be calculated from the particle positions and box dimensions.
    rcutoff : float, optional
        The cutoff fraction of the maximum radius. Default is 0.9, meaning g(r) will be computed up to 90% of the half-box length.
    eps : float, optional
        A small value used to avoid edge cases in distance comparisons (default is 1e-15).

    Returns
    -------
    g_r : (n_radii,) np.array
        Radial distribution function values g(r) corresponding to each radius.
    radii : (n_radii,) np.array
        Array of radii at which g(r) is computed.
    """
    g_r = []
    for ts in particles_timestep:
        g_r_tmp, radii = gr(ts, box_length, dr, rho, rcutoff, eps )
        g_r.append(g_r_tmp)
    
    # Average g(r) over all timesteps
    g_r = np.mean(g_r, axis=0)

    return radii, g_r
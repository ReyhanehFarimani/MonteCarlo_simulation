import numpy as np

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
    # Ensure box_length is a numpy array
    box_length = np.asarray(box_length)
    if box_length.ndim == 0:
        box_length = np.array([box_length, box_length])
    
    # Compute the differences between all pairs of particles
    delta = particles[:, np.newaxis, :] - particles[np.newaxis, :, :]
    
    # Apply minimum image convention
    delta = delta - np.round(delta / box_length) * box_length
    
    # Compute the Euclidean distances
    distances = np.sqrt(np.sum(delta**2, axis=-1))
    
    return distances

def gr(particles, box_length, dr, rho=None, rcutoff=0.9):
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

    Returns
    -------
    g_r : (n_radii,) np.array
        Radial distribution function values g(r) corresponding to each radius.
    radii : (n_radii,) np.array
        Array of radii at which g(r) is computed.
    """

    # Ensure particles are in the correct format (N, 2)
    particles = np.asarray(particles)

    # Ensure box_length is a NumPy array
    box_length = np.asarray(box_length)
    if box_length.ndim == 0:
        box_length = np.array([box_length, box_length])

    # Calculate the maximum radius to consider, limited by the rcutoff
    r_max = (np.min(box_length) / 2) * rcutoff
    
    # Generate an array of radii where g(r) will be computed
    radii = np.arange(0, r_max + dr, dr)

    # Determine the number of particles and dimensionality (should be 2)
    N, d = particles.shape
    
    assert d == 2
    
    # If density (rho) is not provided, calculate it from the number of particles and the box area
    if rho is None:
        rho = N / np.prod(box_length)  # Density in 2D (particles per unit area)

    # Calculate all pairwise distances using the minimum image convention
    distances = minimum_image_distance(particles, box_length)
    
    # Only consider the upper triangle of the distance matrix to avoid double counting
    i_upper = np.triu_indices(N, k=1)
    distances = distances[i_upper]

    # Calculate the histogram of distances
    g_r, bin_edges = np.histogram(distances, bins=radii)
    # Convert g_r to float to avoid casting issues
    g_r = g_r.astype(np.float64)
    
    # Calculate the bin centers (midpoints of the bins)
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2.0
    print(bin_centers)
    # Compute the area of each annular shell
    shell_areas = np.pi * ((bin_edges[1:]**2) - (bin_edges[:-1]**2))
    shell_areas = np.pi * bin_centers * dr
    # Normalize g(r) by the expected number of particles in a random distribution
    normalization = N * rho * shell_areas

    # Handle cases where normalization is very small to avoid division by zero or huge values
    g_r = np.divide(g_r, normalization, out=np.zeros_like(g_r), where=normalization > 0)
    
    return g_r, bin_centers



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
        g_r_tmp, radii = gr(ts, box_length, dr = dr, rho=rho, rcutoff=rcutoff )
        g_r.append(g_r_tmp)
    
    # Average g(r) over all timesteps
    if len(g_r)>1:
        g_r = np.mean(g_r, axis=0)
    else:
        [g_r] = g_r
    return radii, g_r




def voronoi_neighborhood_list(points):
    from scipy.spatial import Voronoi   
    """
    Generates the Voronoi neighborhood list for a set of 2D points.

    Parameters
    ----------
    points : np.ndarray
        Array of points with shape (N, 2), where N is the number of points.

    Returns
    -------
    neighbors : dict
        A dictionary where each key is a point index, and the value is a list of neighboring point indices.
    """
    vor = Voronoi(points)
    neighbors = {i: set() for i in range(len(points))}
    
    for simplex in vor.ridge_points:
        i, j = simplex
        neighbors[i].add(j)
        neighbors[j].add(i)
    
    # Convert sets to sorted lists for consistency
    for i in neighbors:
        neighbors[i] = sorted(neighbors[i])
    
    return neighbors
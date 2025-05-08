import numpy as np
import freud
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

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
    box_size = np.array([box_size[0], box_size[1]])
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


def compute_psi(all_points, box_size, order_number=6):
    """
    Compute the average and standard deviation of the hexatic or other order parameter
    (e.g., sigma_6 or sigma_4) over a number of samples.

    Parameters:
    ----------
    points: list of np.ndarray
        Array of points with shape (N, 2), where N is the number of points.
    box_size : tuple of float
        The size of the box in the x and y dimensions.

    order_number : int
        Order of the orientational order parameter (e.g., 6 for hexatic).

    Returns:
    -------
    mean_sigma : float
        Mean value of the order parameter.
    std_sigma : float
        Standard deviation of the order parameter.
    """
    sigma_abs = []

    for positions in all_points:
        # Compute Voronoi neighbors and hexatic order parameter
        voro = freud.locality.Voronoi()
        voro.compute(
            system=({"Lx": box_size[0], "Ly": box_size[1], "dimensions": 2}, positions)
        )

        op = freud.order.Hexatic(k=order_number)
        op.compute(
            system=({"Lx": box_size[0], "Ly": box_size[1], "dimensions": 2}, positions),
            neighbors=voro.nlist,
        )

        sigma_abs.append(np.mean(np.abs(op.particle_order)))

    return np.mean(sigma_abs), np.std(sigma_abs)/(len(sigma_abs) - 1)**0.5


def compute_g6(all_positions, box_size, r_max=10.0, nbins=100):
    """
    Compute the g6 function for bond-orientational order using Freud.

    Parameters:
    - system_snap: A freud box or particle snapshot containing particle positions.
    - r_max: Maximum distance to consider for calculating g6.
    - nbins: Number of bins for distance calculation.

    Returns:
    - r_bins: Radial distances.
    - g6_vals: g6 values for each radial distance.
    """
    # if np.min(box)/2 < r_max:
    #     raise ValueError
    
    r_bins = np.linspace(0, r_max, nbins)
    g6_vals = np.zeros_like(r_bins)
    nrmalize_vals = np.zeros_like(r_bins)
    box = freud.box.Box(Ly = box_size[1], Lx= box_size[0])
    CF = []
    for points in all_positions:
        # Compute Voronoi neighbors and hexatic order parameter
        voro = freud.locality.Voronoi()
        voro.compute((box, points))
        op = freud.order.Hexatic(k=6)
        op.compute(
            system=({"Lx": box_size[0], "Ly": box_size[1], "dimensions": 2}, points),
            neighbors=voro.nlist,
        )
        values = op.particle_order
        
        cf = freud.density.CorrelationFunction(bins=nbins, r_max=r_max)
        
        cf.compute(
                    system=(box, points), values=values, query_points=points, query_values=values
                )
        g6_vals = np.real(cf.correlation)
        CF.append(g6_vals)
    return r_bins, CF



def compute_gG(all_positions, box_size, density = 0, r_max = 10.0, n_bins = 100):
    """
    Compute g_G(r) for the given reciprocal lattice vector G using freud.density.CorrelationFunction.

    Parameters
    ----------
    all_positions : list of np.ndarray
        List of particle positions for each time step.
    box_size : tuple of float
        The size of the simulation box (Lx, Ly).
    r_max : float, optional
        Maximum distance for computing g_G(r).
    n_bins : int, optional
        Number of bins for the radial distance.

    Returns
    -------
    r_bins : np.ndarray
        The radial distance bins.
    gG_r : np.ndarray
        The computed g_G(r) values (complex).
    """
    # Compute the average RDF over all timesteps
    if density ==0:
        avg_rdf, r_bins = average_rdf_over_trajectory(all_positions, box_size, dr=0.1, rcutoff=0.9)
        # plt.plot(r_bins, avg_rdf)
        # plt.show()
        # Find the first peak of the RDF using scipy's find_peaks
        peaks, _ = find_peaks(avg_rdf)

        if len(peaks) == 0:
            raise ValueError("No peaks found in the RDF. Check the input data or parameters.")

        # Extract the lattice constant as the radius corresponding to the first peak
        lattice_constant = r_bins[peaks[0]]
    
    else:
        lattice_constant = np.sqrt(2 / density / 3**0.5)
    
    print(f"lattice_constan is {lattice_constant}.")
    
    
    Lx = box_size[0]
    Ly = box_size[1]
    # Compute reciprocal lattice vectors based on lattice_constant
    ratio = (3 ** 0.5) / 2
    if np.abs(Lx / Ly - ratio) < 1e-4:
        b1 = np.array([0, 4 * np.pi / (3 ** 0.5) / lattice_constant, 0])
        b2 = np.array([1, -1 / (3 ** 0.5), 0]) * (2 * np.pi / lattice_constant)
    elif np.abs(Ly / Lx - ratio) < 1e-4:
        b1 = np.array([4 * np.pi / (3 ** 0.5) / lattice_constant, 0, 0])
        b2 = np.array([-1 / (3 ** 0.5), 1, 0]) * (2 * np.pi / lattice_constant)
    else:
        raise ValueError("Box does not match the triangular lattice symmetry.")

    # Define the wave vector G (e.g., b1 or b2, or a combination)
    G_vector = b1 + b2 

    # Initialize CorrelationFunction
    cf = freud.density.CorrelationFunction(bins=n_bins, r_max=r_max)

    # Initialize the box
    box = freud.box.Box(Lx=Lx, Ly=Ly)

    # Accumulate g_G(r) over all time steps
    gG_values = []
    for positions in all_positions:
        # Compute the dot product of G_vector with particle positions
        values = np.exp(1j * np.dot(positions, G_vector))

        # Compute g_G(r) using freud's CorrelationFunction
        cf.compute(
            system=(box, positions), values=values, query_points=positions, query_values=values
        )

        gG_values.append(cf.correlation)

    # Average g_G(r) over all time steps
    gG_avg = np.mean(gG_values, axis=0)

    return cf.bin_centers, gG_avg

def sub_system_translational(all_positions, box_size, Lb):
    
    
# Compute the average RDF over all timesteps
    avg_rdf, r_bins = average_rdf_over_trajectory(all_positions, box_size, dr=0.01, rcutoff=0.9)
    # plt.plot(r_bins, avg_rdf)
    # plt.show()
    # Find the first peak of the RDF using scipy's find_peaks
    peaks, _ = find_peaks(avg_rdf)

    if len(peaks) == 0:
        raise ValueError("No peaks found in the RDF. Check the input data or parameters.")

    # Extract the lattice constant as the radius corresponding to the first peak
    lattice_constant = r_bins[peaks[0]]
    print(f"lattice_constan is {lattice_constant}.")
    
    Lx = box_size[0]
    Ly = box_size[1]
    # Compute reciprocal lattice vectors based on lattice_constant
    ratio = (3 ** 0.5) / 2
    if np.abs(Lx / Ly - ratio) < 1e-4:
        b1 = np.array([0, 4 * np.pi / (3 ** 0.5) / lattice_constant, 0])
        b2 = np.array([1, -1 / (3 ** 0.5), 0]) * (2 * np.pi / lattice_constant)
    elif np.abs(Ly / Lx - ratio) < 1e-4:
        b1 = np.array([4 * np.pi / (3 ** 0.5) / lattice_constant, 0, 0])
        b2 = np.array([1 / (3 ** 0.5), -1, 0]) * (2 * np.pi / lattice_constant)
    else:
        raise ValueError("Box does not match the triangular lattice symmetry.")

    # Define the wave vector G (e.g., b1 or b2, or a combination)
    G_vector = b1 + b2 
    Psi_G = []
    if Lb == 1:
        for points in all_positions:
            values = np.mean(np.exp(1j * np.dot(points, G_vector)))
            Psi_G.append(values)
        return Psi_G
    sub_systems_0 = np.arange(Lb/2, 1.0 - Lb, Lb/2)

    sub_systems_1 = sub_systems_0 + Lb

    for points in all_positions:
        for sub_system_0, sub_system_1 in zip(sub_systems_0, sub_systems_1):
            x_bound = [sub_system_0 * Lx, sub_system_1 * Lx]
            y_bound = [sub_system_0 * Ly, sub_system_1 * Ly]
            bounded_point = points[points[:, 0]<x_bound[1]]
            bounded_point = bounded_point[bounded_point[:,0]>x_bound[0]]
            bounded_point = bounded_point[bounded_point[:,1]<y_bound[1]]
            bounded_point = bounded_point[bounded_point[:,1]<y_bound[0]]
            
            ### computing \psi_G
            if bounded_point != np.array([]):
                values = np.mean(np.exp(1j * np.dot(bounded_point, G_vector)))
                if ~np.isnan(values).all():
                    Psi_G.append(values)
    
    
    return Psi_G


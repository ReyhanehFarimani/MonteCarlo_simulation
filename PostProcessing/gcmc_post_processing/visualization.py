import matplotlib.pyplot as plt

def plot_gr(radii, g_r, title="Radial Distribution Function g(r)", xlabel="r", ylabel="g(r)", output_file=None):
    """
    Plots the radial distribution function g(r).

    Parameters
    ----------
    radii : np.array
        Array of radii at which g(r) is computed.
    g_r : np.array
        Radial distribution function values g(r).
    title : str, optional
        Title of the plot.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    output_file : str, optional
        If provided, the plot will be saved to this file.

    Returns
    -------
    None
    """
    plt.figure(figsize=(8, 6))
    plt.plot(radii, g_r, linestyle='-', color='blue')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    
    if output_file:
        plt.savefig(output_file)
        print(f"Plot saved to {output_file}")
    
    plt.show()
    
    
    

def plot_voronoi(points):
    """
    Plots the Voronoi diagram for a set of 2D points.

    Parameters
    ----------
    points : np.ndarray
        Array of points with shape (N, 2), where N is the number of points.

    Returns
    -------
    None
    """
    import matplotlib.pyplot as plt
    from scipy.spatial import voronoi_plot_2d, Voronoi
    vor = Voronoi(points)
    
    # Plot Voronoi diagram
    fig, ax = plt.subplots()
    voronoi_plot_2d(vor, ax=ax, show_vertices=False, line_colors='orange', line_width=2)
    
    # Plot points
    ax.plot(points[:, 0], points[:, 1], 'o', color='black')
    
    ax.set_xlim([points[:, 0].min() - 1, points[:, 0].max() + 1])
    ax.set_ylim([points[:, 1].min() - 1, points[:, 1].max() + 1])
    
    ax.set_aspect('equal')
    ax.set_title('Voronoi Diagram')
    plt.show()





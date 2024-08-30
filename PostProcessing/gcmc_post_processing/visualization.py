import matplotlib.pyplot as plt
import pandas as pd

# Step 1: Set Global Visualization Settings
plt.rcParams.update({
    'font.family': 'serif',
    # 'font.serif': ['Times New Roman'],
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.figsize': (10, 6),
    'figure.dpi': 100,
    'savefig.dpi': 300,
    'savefig.format': 'png',
    'grid.linestyle': '--',
    'grid.color': 'gray',
    'grid.alpha': 0.2,
})

def plot_data(ax, df, x_column, y_column, yerr_column=None, xlabel="", ylabel="", title="", plot_type="errorbar", **kwargs):
    """
    General function to plot data on a provided ax.

    Parameters:
    - ax: Matplotlib axis object to plot on.
    - df: DataFrame containing the data to plot.
    - x_column: Name of the column in df to use for the x-axis.
    - y_column: Name of the column in df to use for the y-axis.
    - yerr_column: (Optional) Name of the column in df to use for y-axis error bars.
    - xlabel: (Optional) Label for the x-axis.
    - ylabel: (Optional) Label for the y-axis.
    - title: (Optional) Title for the plot.
    - plot_type: Type of plot ("errorbar", "scatter", "line").
    - **kwargs: Additional keyword arguments to pass to the plotting function.
    """
    
    if plot_type == "errorbar":
        ax.errorbar(df[x_column], df[y_column], yerr=df[yerr_column] if yerr_column else None, fmt=kwargs.get('fmt', 'o-'), label=kwargs.get('label', 'Data'))
    elif plot_type == "scatter":
        ax.scatter(df[x_column], df[y_column], label=kwargs.get('label', 'Data'))
    elif plot_type == "line":
        ax.plot(df[x_column], df[y_column], label=kwargs.get('label', 'Data'))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True)
    ax.legend()
    
def plot_rdf(rdf, r_bins):
    """
    Plots the radial distribution function g(r).

    Parameters
    ----------
    rdf : np.ndarray
        The computed g(r) values.
    r_bins : np.ndarray
        The radii corresponding to the g(r) values.

    Returns
    -------
    None
    """
    plt.figure()
    plt.plot(r_bins, rdf, '-o')
    plt.xlabel('r')
    plt.ylabel('g(r)')
    plt.title('Radial Distribution Function (RDF)')
    plt.grid(True)
    plt.show()

def plot_voronoi_with_ids(points, box_size, neighbors):
    """
    Plots the Voronoi diagram with particle IDs and shows the Voronoi neighbors using freud.

    Parameters
    ----------
    points : np.ndarray
        Array of points with shape (N, 2), where N is the number of points.
    box_size : tuple of float
        The size of the box in the x and y dimensions.
    neighbors : dict
        A dictionary where each key is a point index, and the value is a list of neighboring point indices.

    Returns
    -------
    None
    """
    # Define the box and Voronoi analysis
    box = freud.box.Box(Lx=box_size[0], Ly=box_size[1], Lz = 0.0)
    voronoi = freud.locality.Voronoi(box)
    
    # Compute the Voronoi diagram
    voronoi.compute(system=(box, points))
    
    # Plot Voronoi diagram
    fig, ax = plt.subplots()
    for poly in voronoi.polytopes:
        polygon = np.asarray(poly)
        polygon = np.vstack([polygon, polygon[0]])  # Close the polygon
        ax.plot(polygon[:, 0], polygon[:, 1], color='orange', lw=2)
    
    # Plot points and their IDs
    for i, point in enumerate(points):
        ax.plot(point[0], point[1], 'o', color='black')
        ax.text(point[0] + 0.1, point[1] + 0.1, str(i), color='blue', fontsize=12)
        
        # Annotate Voronoi neighbors
        neighbor_ids = neighbors[i]
        neighbor_text = f"{', '.join(map(str, neighbor_ids))}"
        ax.text(point[0] + 0.1, point[1] - 0.3, neighbor_text, color='red', fontsize=8)

    ax.set_xlim([0, box_size[0]])
    ax.set_ylim([0, box_size[1]])
    
    ax.set_aspect('equal')
    ax.set_title('Voronoi Diagram with Periodic Boundary Conditions (freud)')
    plt.show()

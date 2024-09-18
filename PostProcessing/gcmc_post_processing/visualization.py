import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import warnings
import freud

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


def user_defined_voronoi_plot(box, polytopes, ax=None, color_by_sides=True, cmap=None, color_array = None):
    """Helper function to draw 2D Voronoi diagram.

    Args:
        box (:class:`freud.box.Box`):
            Simulation box.
        polytopes (:class:`numpy.ndarray`):
            Array containing Voronoi polytope vertices.
        ax (:class:`matplotlib.axes.Axes`): Axes object to plot.
            If :code:`None`, make a new axes and figure object.
            (Default value = :code:`None`).
        color_by_sides (bool):
            If :code:`True`, color cells by the number of sides.
            If :code:`False`, random colors are used for each cell.
            (Default value = :code:`True`).
        cmap (str):
            Colormap name to use (Default value = :code:`None`).

    Returns:
        :class:`matplotlib.axes.Axes`: Axes object with the diagram.
    """
    from matplotlib import cm
    from matplotlib.collections import PatchCollection
    from matplotlib.colorbar import Colorbar
    from matplotlib.patches import Polygon
    from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

    if ax is None:
        fig = plt.figure()
        ax = fig.subplots()

    # Draw Voronoi polytopes
    patches = [Polygon(poly[:, :2]) for poly in polytopes]
    patch_collection = PatchCollection(patches, edgecolors="black", alpha=0.4)

    if color_by_sides:
        colors = np.array([len(poly) for poly in polytopes])
        num_colors = np.ptp(colors) + 1
    else:
        colors = np.random.RandomState().permutation(np.arange(len(patches)))
        num_colors = np.unique(colors).size
    if color_array is not None:
        colors = color_array
        num_colors = np.unique(colors).size

    # Ensure we have enough colors to uniquely identify the cells
    if cmap is None:
        if color_by_sides and num_colors <= 10:
            cmap = "tab10"
        else:
            if num_colors > 20:
                warnings.warn(
                    "More than 20 unique colors were requested. "
                    "Consider providing a colormap to the cmap "
                    "argument.",
                    UserWarning,
                )
            cmap = "tab20"
    cmap = cm.get_cmap(cmap, num_colors)
    bounds = np.arange(-1, 2)

    patch_collection.set_array(np.array(colors))
    patch_collection.set_cmap(cmap)
    patch_collection.set_clim(bounds[0], bounds[-1] )
    ax.add_collection(patch_collection)

    # Draw box
    corners = [[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0]]
    # Need to copy the last point so that the box is closed.
    corners.append(corners[0])
    corners = box.make_absolute(corners)[:, :2]
    ax.plot(corners[:, 0], corners[:, 1], color="k")

    # Set title, limits, aspect
    ax.set_title("Voronoi Diagram")
    ax.set_xlim((np.min(corners[:, 0]), np.max(corners[:, 0])))
    ax.set_ylim((np.min(corners[:, 1]), np.max(corners[:, 1])))
    ax.set_aspect("equal", "datalim")

    # Add colorbar for number of sides
    if color_by_sides:
        ax_divider = make_axes_locatable(ax)
        cax = ax_divider.append_axes("right", size="7%", pad="10%")
        cb = Colorbar(cax, patch_collection)
        cb.set_label(r"$\vec{\Psi}(\rho, f).\vec{\psi}(\vec{r})$")
        cb.set_ticks(bounds)
    return ax


def plot_voro(ts, box, filename, label):
    fig, ax = plt.subplots(1, 1, figsize=(20, 20))
    points = ts - box/2
    points[:, 2] = 0
    box = freud.box.Box(Ly = box[1], Lx= box[0])
    voro = freud.locality.Voronoi()
    op = freud.order.Hexatic(k=6)
    cells = voro.compute((box, points)).polytopes
    op.compute( (box , points),neighbors=voro.nlist)
    sigma_6 = op.particle_order
    S = sigma_6.mean()
    colors = np.array([(s.real * S.real + s.imag * S.imag)/(np.absolute(s)*np.absolute(S)) for s in sigma_6])
    user_defined_voronoi_plot(box = box, polytopes= cells, ax = ax, color_array=colors, cmap = "terrain")
    ax.quiver(points[:, 0], points[:, 1], np.real(sigma_6),np.imag(sigma_6),color = "k",  label = label, scale = 100, alpha = 0.6)
    ax.legend(loc = "upper center")
    plt.savefig(filename + ".pdf", dpi = 800)
    plt.savefig(filename + ".eps", dpi = 800)
    plt.show()
from gcmc_post_processing.analysis import voronoi_neighborhood_list
from gcmc_post_processing.visualization import plot_voronoi_with_ids
import numpy as np

# Example usage:
points = np.random.rand(10, 3) * np.array([10, 10, 0.0])  # Generate 10 random 2D points within a 10x10 box
box_size = (10, 10)  # Define the size of the box

# Compute Voronoi neighborhood list using freud
neighbors = voronoi_neighborhood_list(points, box_size)

# Plot Voronoi diagram with particle IDs and neighbors using freud
plot_voronoi_with_ids(points, box_size, neighbors)

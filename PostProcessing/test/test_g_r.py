import numpy as np
from gcmc_post_processing.data_loader import load_xyz
from gcmc_post_processing.analysis import average_rdf_over_trajectory, compute_rdf
from gcmc_post_processing.utils import generate_fake_honeycomb_structure, generate_fake_crystal_structure, generate_random_particles
from gcmc_post_processing.visualization import plot_gr
import matplotlib.pyplot as plt
import freud

if __name__ == "__main__":
    box = np.array([100, 100])
    generate_random_particles(box_size=box, output_file="random.xyz", num_particles=10000)
    particles = load_xyz("random.xyz")
    # print(particles[0].shape)
    
    gr, radii = compute_rdf(particles[0], box, dr = 0.01, rcutoff=0.9, )
    plot_gr(radii, gr)
    generate_fake_crystal_structure(grid_size=box, output_file="crystal.xyz")
    particles = load_xyz("crystal.xyz")
    gr, radii = average_rdf_over_trajectory(particles, box, dr = 0.1, rcutoff=0.9, )
    plot_gr(radii, gr)

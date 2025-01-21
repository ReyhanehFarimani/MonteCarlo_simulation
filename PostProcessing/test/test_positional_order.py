import numpy as np
from gcmc_post_processing import load_xyz, compute_gG, average_rdf_over_trajectory, generate_fake_triangular_lattice, generate_fake_triangular_lattice_defected
import matplotlib.pyplot as plt

def test_the_lattice_constant():
    box = np.array([100 * 3, 200 * np.sqrt(3)])
    generate_fake_triangular_lattice(lattice_constant= 2, output_file="triangular_2.xyz", grid_size=box)
    p = load_xyz("triangular_2.xyz")
    a = compute_gG(p, box)
    generate_fake_triangular_lattice(lattice_constant= 1, output_file="triangular_1.xyz", grid_size=box)
    p = load_xyz("triangular_1.xyz")
    a = compute_gG(p, box)
    
def test_the_G_g():    
    box = np.array([100 * 3, 200 * np.sqrt(3)])
    generate_fake_triangular_lattice(lattice_constant= 2, output_file="triangular_2.xyz", grid_size=box)
    p = load_xyz("triangular_2.xyz")
    r1, g = compute_gG(p, box, r_max=100, n_bins=10000)
    plt.plot(r1, np.real(g), label = "0 defect")
    generate_fake_triangular_lattice_defected(lattice_constant= 2, output_file="triangular_2.xyz", grid_size=box, defect=0.05)
    p = load_xyz("triangular_2.xyz")
    r1, g = compute_gG(p, box, r_max=100, n_bins=10000)
    plt.plot(r1, np.real(g), label = "0.05 defect")
    generate_fake_triangular_lattice_defected(lattice_constant= 2, output_file="triangular_2.xyz", grid_size=box, defect=0.1)
    p = load_xyz("triangular_2.xyz")
    r1, g = compute_gG(p, box, r_max=100, n_bins=10000)
    plt.plot(r1, np.real(g), label = "0.1 defect")
    generate_fake_triangular_lattice_defected(lattice_constant= 2, output_file="triangular_2.xyz", grid_size=box, defect=0.3)
    p = load_xyz("triangular_2.xyz")
    r1, g = compute_gG(p, box, r_max=100, n_bins=10000)
    plt.plot(r1, np.real(g), label = "0.3 defect")
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(0.001, 10.0)
    plt.legend()
    plt.show()
if __name__ == "__main__":
    # test_the_lattice_constant()
    test_the_G_g()
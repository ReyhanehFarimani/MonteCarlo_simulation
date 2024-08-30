# Test Scenarios

This document outlines the test scenarios used to validate the Grand Canonical Monte Carlo (GCMC) simulation code.

For running the simulations by yourself please follow this procedure:
```
cd sys_tests/
make clean
make all
mpirun -np 'number-of-cores' ./sys_tests

```
## 1. Ideal Gas

### Description

In the grand canonical ensemble, the chemical potential \( \mu \) controls the exchange of particles between the system and a reservoir. For an ideal gas, the fugacity \(z\) is defined as:

$$
z = e^{\beta \mu}
$$

where \( \beta = \frac{1}{k_B T} \). Fugacity can be interpreted as the "effective pressure" exerted by the particles in the system, simplifying the calculation of the number of particles.

The average number of particles \( \langle N \rangle \) in the system can be derived from the fugacity:

$$
\langle N \rangle = \frac{V z}{\lambda^2}
$$

where:
- \( V \) is the volume of the system,
- \( \lambda \) is the thermal de Broglie wavelength, which depends on the temperature and mass of the particles.

This equation indicates that the number of particles in an ideal gas increases linearly with the volume and fugacity, and inversely with the square of the thermal wavelength.

The pressure \( P \) in the system can also be related to the average number of particles \( \langle N \rangle \) and the temperature:

$$
P = \frac{k_B T \langle N \rangle}{V}
$$

This is consistent with the ideal gas law, showing that the pressure is directly proportional to the number of particles and temperature, and inversely proportional to the volume.

### Summary

This ideal gas test serves to validate the Grand Canonical Monte Carlo (GCMC) simulation by comparing the simulation results for the average number of particles \( \langle N \rangle \) and the pressure \( P \) with the theoretical predictions. The fugacity is a key factor in determining the number of particles in the system, and the pressure is expected to follow the ideal gas law.

### Results

[Results Plot](test_results/N_P_plot_ideal.html)

The discrepancies observed between the theoretical predictions and the simulation results were due to the fact that the program is designed in a way that the maximum density of particle to be set 3, this cause a dispersety in the number of particles particulary at high \(\mu\).


---

## 2. Lennard-Jones Fluid (NVT simulation, and GCMC simulation)

### Description

This test simulates a fluid using the Lennard-Jones potential, a common model for intermolecular interactions.

### Expected Outcome

DOI: <10.1103/PhysRevE.102.062101>

### Results:

---

## 3. Phase diagram of Logarithmic interacting particles

### Description

This test simulates a system of athermal stars.

### Expected Outcome

DOI: <https://doi.org/10.1039/C8SM02100G>

### Results:


---


### Additional Considerations

- **Parameter Sensitivity**: Ensure that the simulation parameters, such as temperature, chemical potential, and volume, are varied systematically to test the robustness of the simulation.
- **Performance Metrics**: Evaluate the performance of the simulation in terms of computational efficiency and convergence rates, especially for more complex scenarios like the Lennard-Jones fluid.

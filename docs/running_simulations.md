## Basic Usage

To use the program one should make an 'input.txt' in a directory, with the following information and structure:
```
# Simulation box dimensions
boxLengthX 25.0
boxLengthY 25.0

# Number of particles in the simulation
numParticles 100

# Particle initialization: 1 for random placement, 0 for grid placement
randomPlacement 1

# Temperature of the simulation
temperature 1

# Number of simulation steps
numSteps 1000000000

# Time step for the simulation
timeStep 0.1

# Output frequency: how often to log the results
outputFrequency 10000

# Equilibration time: number of steps before starting to log
equilibrationTime 1000

# Cutoff squared distance for potential calculations
r2cut 2.0

# Type of potential: LennardJones, WCA, Yukawa, AthermalStar, Ideal
potentialType AthermalStar
# polymer functionality in case of Athermal star polymer
f 10
# Type of simulation: NVT, GCMC
simulationType GCMC
# chemical potential in case of GCMC simulation
mu 0.1

seed 1234

# Use cell list: 1 for yes, 0 for no
useCellList 1

# Frequency for updating the cell list
cellListUpdateFrequency 1000

# File names for logging
positionFile particle_positions.xyz
dataFile simulation_data.dat

```
### Running a GCMC Simulation

To run a GCMC simulation, use the following command:
```
cd /path_to_input
./path_to_installation/simulation 
```

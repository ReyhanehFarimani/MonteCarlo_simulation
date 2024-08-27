# GCMC Simulation

Welcome to the Grand Canonical Monte Carlo (GCMC) simulation project. This repository contains code for performing two-dimensional Monte Carlo simulations of interacting particles. The project supports both NVT (canonical) and GCMC (grand canonical) ensemble simulations.

## Table of Contents

- [About the Project](#about-the-project)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [Contact](#contact)

## About the Project

This project is designed to simulate interacting particles in a two-dimensional space using Monte Carlo methods. The main focus is on the Grand Canonical Monte Carlo (GCMC) simulation, which allows for the exchange of particles with a reservoir at a fixed chemical potential. The project is written in C++ with Python integration for post-processing and analysis.

## Features

- **NVT and GCMC Simulations**: Perform simulations in both the canonical and grand canonical ensembles.
- **Customizability**: Easily extend the code with your own Monte Carlo moves or interaction potentials.
- **High Performance**: The code is optimized for computational efficiency, leveraging C++.
- **Python Integration**: Python scripts are provided for post-processing simulation results and creating interactive plots.

## Installation

### Prerequisites

Before you begin, ensure you have the following installed:

- A C++ compiler (e.g., `g++`, `clang`)
- MPI (Message Passing Interface) for parallel computation
- Python 3.x with the necessary libraries (`numpy`, `matplotlib`, etc.)

### Installation Steps

   ```
   git clone https://github.com/ReyhanehFarimani/MonteCarlo_simulation.git
   cd MonteCarlo_simulation
   make all
   ```
### Usage:
  Please read the documentry provided at <https://reyhanehfarimani.github.io/MonteCarlo_simulation/>

### Contribuation:

   Contributions are welcome!

### Contact

For more information, please contact:

Reyhaneh Afghahi Farimani

Email: <reyhaneh.afghahi.farimani@univie.ac.at>

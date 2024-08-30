# GCMC Simulation Documentation

Welcome to the documentation for the Grand Canonical Monte Carlo (GCMC) simulation project. This resource provides comprehensive guidance on using the simulation software, including the API documentation, installation instructions, and steps to run the program.

## Overview

This project is designed for two-dimensional Monte Carlo simulations of interacting particles. Currently, the program supports both the NVT (canonical) and GCMC (grand canonical) ensembles. The codebase is structured to be highly modular, allowing users to easily extend the functionality by adding custom Monte Carlo moves or interaction potentials.

The simulation is implemented in C++ with additional Python scripting for post-processing and extended functionality. To run this software, you will need to have a C++ compiler, MPI (Message Passing Interface) for parallelization -Still in progress-, and a Python environment set up.

## Features

- **NVT and GCMC Simulations**: Perform simulations in both canonical and grand canonical ensembles.
- **Modular Design**: Easily extend the codebase with custom Monte Carlo moves or interaction potentials.
- **High Performance**: Optimized for efficient computation using C++.
- **Python Integration**: Utilize Python for post-processing, analysis, and additional scripting.

## Installation

Follow the [Installation Guide](installation.md) to set up the necessary environment and compile the code.

## Running Simulations

Learn how to configure and execute simulations by visiting the [Running Simulations](running_simulations.md) section.

## Test Scenarios

Explore various test cases used to validate the simulation code in the [Test Scenarios](test_scenarios.md) section.

## Examples

For practical examples of using the software, refer to the [Examples](example.md) section.

## API Documentation

Detailed API documentation is available [here](doxygen/html/index.html), generated using Doxygen.

## Credits

This program was developed by Reyhaneh Afghahi Farimani. For inquiries or further information, please contact me at:

<reyhaneh.afghahi.farimani@univie.ac.at>

## Table of Contents

- [Installation](installation.md)
- [Running Simulations](running_simulations.md)
- [Test Scenarios](test_scenarios.md)
- [Examples](example.md)
- [API Documentation](doxygen/html/index.html)

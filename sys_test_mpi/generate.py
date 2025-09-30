#!/usr/bin/env python3
"""
Generate input folders/files for 256^2 Lennard-Jones systems at specified densities,
with box ratio Lx:Ly = sqrt(3):2 and rcut = 2.5, T = 10.0.

Creates one folder per density (e.g. rho_1.190/) containing input.txt.
"""

import os
import math
import pathlib

# ---- Config ----
densities = [1.19, 1.20, 1.21, 1.22, 1.205, 1.215]
N = 256 * 256               # number of particles
T = 10.0
rcut = 2.5

# MC defaults (adjust if needed)
nSteps = 200000
eSteps = 50
outputFreq = 100
cellUpdateFreq = 50
delta = 0.10

# ---- Helpers ----
def dims_for_density(rho: float):
    """
    Given density rho (N/V) and the fixed aspect ratio Lx:Ly = sqrt(3):2,
    compute box lengths Lx, Ly.
    Area V = N / rho
    Lx = sqrt(3)*a, Ly = 2*a, so Lx*Ly = (2*sqrt(3)) * a^2 = V  =>  a = sqrt(V/(2*sqrt(3)))
    """
    V = N / rho
    a = math.sqrt(V / (2.0 * math.sqrt(3.0)))
    Lx = math.sqrt(3.0) * a
    Ly = 2.0 * a
    return V, Lx, Ly

def write_input(path: pathlib.Path, params: dict):
    with open(path, "w") as f:
        f.write(f"Lx {params['Lx']:.10f}\n")
        f.write(f"Ly {params['Ly']:.10f}\n")
        f.write(f"N {params['N']}\n")
        f.write(f"rcut {params['rcut']}\n")
        f.write(f"T {params['T']}\n")
        f.write(f"nSteps {params['nSteps']}\n")
        f.write(f"eSteps {params['eSteps']}\n")
        f.write(f"outputFreq {params['outputFreq']}\n")
        f.write(f"cellUpdateFreq {params['cellUpdateFreq']}\n")
        f.write(f"potential LennardJones\n")
        f.write(f"out_xyz {params['out_xyz']}\n")
        f.write(f"out_data {params['out_data']}\n")
        f.write(f"delta {params['delta']}\n")
        f.write(f"mu 0.0\n")
        f.write(f"f 0.0\n")
        f.write(f"alpha 0.0\n")
        f.write(f"A_0 0.0\n")
        f.write(f"kappa 1.0\n")
        f.write(f"seed {params['seed']}\n")
        f.write(f"position_file \n")
        f.write(f"ensemble NVT\n")
        f.write(f"P 0.0\n")
        f.write(f"delta_V 0.0\n")

def main():
    base = pathlib.Path(".").resolve()
    created = 0
    for i, rho in enumerate(densities):
        V, Lx, Ly = dims_for_density(rho)
        name = f"rho_{rho:.3f}"
        d = base / name
        d.mkdir(parents=True, exist_ok=True)

        seed = 10000 + i
        params = {
            "Lx": Lx,
            "Ly": Ly,
            "N": N,
            "rcut": rcut,
            "T": T,
            "nSteps": nSteps,
            "eSteps": eSteps,
            "outputFreq": outputFreq,
            "cellUpdateFreq": cellUpdateFreq,
            "out_xyz": f"traj_{rho:.3f}.xyz",
            "out_data": f"thermo_{rho:.3f}.csv",
            "delta": delta,
            "seed": seed,
        }
        write_input(d / "input.txt", params)

        with open(d / "README.txt", "w") as rf:
            rf.write(
                f"rho = {rho:.3f}\n"
                f"N = {N}\n"
                f"V = {V:.10f}\n"
                f"Lx = {Lx:.10f}\n"
                f"Ly = {Ly:.10f}\n"
                f"ratio Lx:Ly ≈ sqrt(3):2\n"
            )
        created += 1

    print(f"Done.")

if __name__ == "__main__":
    main()

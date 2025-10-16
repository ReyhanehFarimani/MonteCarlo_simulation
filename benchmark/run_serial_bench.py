#!/usr/bin/env python3
"""
Serial Monte Carlo benchmark harness (space-separated input format).

- Runs ../Monte_carlo_serial (default) directly.
- Generates inputs for 3 densities; varies N.
- Chooses serial nSteps to match MPI total attempts (default: equal).
- Writes all inputs, logs, CSVs, plots under ./results_serial/.
"""

import argparse, csv, math, os, subprocess, sys, time
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Benchmark plan (edit as needed) ----
DENSITIES = [0.10, 0.30, 0.60]           # rho = N/(Lx*Ly)
N_BY_RHO = {
    0.10: [2_000, 8_000, 32_000],
    0.30: [2_000, 8_000, 32_000],
    0.60: [2_000, 8_000, 32_000],
}
TRIALS  = 3
TEMP    = 1.0
POTENTIAL = "LennardJones"
RCUT    = 2.5
DELTA   = 0.1 * RCUT
# Serial-only knobs (kept simple & IO-light for timing):
F_OUTPUT_STEP = 100      # effectively disables periodic output in serial code paths
F_UPDATE_CELL = 80         # update cell list every k steps (tune if you like)
SEED    = 1234
ENSEMBLE = "NVT"

def make_input_text(Lx, Ly, N, rcut, T, nSteps, eSteps, outputFreq, cellUpdateFreq,
                    potential, out_xyz, out_data, delta, seed, ensemble):
    # SPACE-SEPARATED format (no '='), per your requirement.
    return f"""Lx {Lx:.8f}
Ly {Ly:.8f}
N {N}
rcut {rcut}
T {T}
nSteps {nSteps}
eSteps {eSteps}
outputFreq {outputFreq}
cellUpdateFreq {cellUpdateFreq}
potential {potential}
out_xyz {out_xyz}
out_data {out_data}
delta {delta}
seed {seed}
ensemble {ensemble}
"""

def run_once(exe, input_path):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = env.get("OMP_NUM_THREADS", "1")
    cmd = [exe, str(input_path)]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, env=env)
    t1 = time.perf_counter()
    return (t1 - t0), proc.returncode, proc.stdout

def ensure_dir(p: Path): p.mkdir(parents=True, exist_ok=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default="../Monte_carlo_serial",
                    help="Path to serial executable (default ../Monte_carlo_serial)")
    ap.add_argument("--outdir", default="results_serial", help="Output dir")
    ap.add_argument("--densities", nargs="*", type=float, default=DENSITIES)
    ap.add_argument("--trials", type=int, default=TRIALS)
    ap.add_argument("--rcut", type=float, default=RCUT)
    ap.add_argument("--delta", type=float, default=DELTA)
    ap.add_argument("--seed", type=int, default=SEED)
    # Matching with MPI:
    ap.add_argument("--mpi-nsteps", type=int, default=500 * 4,
                    help="MPI sweeps used in your MPI runs (default 50)")
    ap.add_argument("--scale-nsteps", type=float, default=1.0,
                    help="Optional scaling if you want serial work != MPI work (default 1.0)")
    ap.add_argument("--f-update-cell", type=int, default=F_UPDATE_CELL,
                    help="Serial cell-list update period (steps)")
    ap.add_argument("--f-output-step", type=int, default=F_OUTPUT_STEP,
                    help="Serial output period (steps)")
    args = ap.parse_args()

    exe = Path(args.exe)
    if not exe.exists():
        print(f"ERROR: Executable not found: {exe.resolve()}", file=sys.stderr)
        sys.exit(1)

    outdir = Path(args.outdir); ensure_dir(outdir)
    inputs_dir = outdir / "inputs"; ensure_dir(inputs_dir)
    logs_dir = outdir / "logs"; ensure_dir(logs_dir)

    # Build (rho, N, Lx, Ly)
    plan = []
    for rho in args.densities:
        Ns = N_BY_RHO.get(rho, [])
        if not Ns:
            print(f"[warn] No N list for rho={rho}, skipping.", file=sys.stderr)
            continue
        for N in Ns:
            L = math.sqrt(N / rho)
            plan.append((rho, N, L, L))

    # Matching rule:
    #   MPI: attempts ≈ N * mpi_nsteps
    #   Serial: attempts = N * serial_nsteps
    # => serial_nsteps = mpi_nsteps * scale_nsteps
    serial_nsteps = max(1, int(round(args.mpi_nsteps * args.scale_nsteps)))

    csv_path = outdir / "serial_bench.csv"
    with open(csv_path, "w", newline="") as fcsv:
        wr = csv.writer(fcsv)
        wr.writerow(["density","N","Lx","Ly","trial","seconds","retcode","log_path","nSteps_serial"])

        for (rho, N, Lx, Ly) in plan:
            # Per-(rho,N) input
            tag_base = f"rho{rho:.2f}_N{N}"
            input_text = make_input_text(
                Lx=Lx, Ly=Ly, N=N, rcut=args.rcut, T=TEMP,
                nSteps=serial_nsteps, eSteps=0,
                outputFreq=args.f_output_step, cellUpdateFreq=args.f_update_cell,
                potential=POTENTIAL,
                out_xyz=f"{tag_base}_serial.xyz",
                out_data=f"{tag_base}_serial.dat",
                delta=args.delta, seed=args.seed, ensemble=ENSEMBLE
            )
            in_path = inputs_dir / f"{tag_base}.inp"
            in_path.write_text(input_text)

            for trial in range(1, args.trials+1):
                tag = f"{tag_base}_t{trial}"
                print(f"[run-serial] {tag} (nSteps={serial_nsteps})")
                sec, ret, out = run_once(str(exe), in_path)
                logp = logs_dir / f"{tag}.log"
                logp.write_text(out)
                wr.writerow([rho, N, f"{Lx:.8f}", f"{Ly:.8f}", trial, f"{sec:.6f}", ret, str(logp), serial_nsteps])
                fcsv.flush()

    # ---- Plotting (runtime vs N) ----
    df = pd.read_csv(csv_path)
    agg = df.groupby(["density","N"])["seconds"].median().reset_index()

    for rho, sub in agg.groupby("density"):
        fig, ax = plt.subplots(figsize=(6,4), dpi=130)
        ax.plot(sub["N"], sub["seconds"], marker="o", color="tab:orange", label="serial")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(f"Serial Runtime vs N  (rho={rho:.2f}, nSteps={serial_nsteps})")
        ax.set_xlabel("N (particles)")
        ax.set_ylabel("Median runtime (s)")
        ax.grid(True, which="both", ls="--", alpha=0.4)
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / f"serial_runtime_vs_N_rho{rho:.2f}.png")

    print("\nDone.")
    print(f"CSV:   {csv_path}")
    print(f"Plots: {outdir}/*.png")
    print(f"Logs:  {logs_dir}/*.log")

if __name__ == "__main__":
    main()

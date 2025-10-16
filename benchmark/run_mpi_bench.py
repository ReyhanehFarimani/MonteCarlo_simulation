#!/usr/bin/env python3
"""
MPI Monte Carlo benchmark harness (no Makefile).

- Runs ../Monte_carlo_mpi via mpirun from inside the benchmarking/ folder.
- Generates inputs for 3 densities; varies N; sweeps np (ranks).
- Writes all inputs, logs, CSVs, plots under ./results/.
"""

import argparse, csv, math, os, subprocess, sys, time
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Benchmark plan (adjust freely) ----
DENSITIES = [0.10, 0.30, 0.60]           # rho = N / (Lx*Ly) in 2D
N_BY_RHO = {
    0.10: [2_000, 8_000, 32_000],
    0.30: [2_000, 8_000, 32_000],
    0.60: [2_000, 8_000, 32_000],
}
PROCS   = [ 4, 8]                   # MPI ranks
TRIALS  = 3
SEED    = 1234
TEMP    = 1.0
POTENTIAL = "LennardJones"               # must match your selectPotentialType()
RCUT    = 2.5
DELTA   = 0.1 * RCUT
NSTEPS  = 500                           # keep constant for comparability
ESTEPS  = 0
OUTPUT_FREQ = 40                          # keep I/O off for clean timing
CELL_UPDATE_FREQ = 50
ENSEMBLE = "NVT"

def make_input_text(Lx, Ly, N, rcut, T, nSteps, eSteps, outputFreq, cellUpdateFreq,
                    potential, out_xyz, out_data, delta, seed, ensemble):
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

def run_once(exe, np, input_path, extra_mpirun):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = env.get("OMP_NUM_THREADS", "1")
    cmd = ["mpirun", "-np", str(np)] + extra_mpirun + [exe, str(input_path)]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, env=env)
    t1 = time.perf_counter()
    return (t1 - t0), proc.returncode, proc.stdout

def ensure_dir(p: Path): p.mkdir(parents=True, exist_ok=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default="../Monte_carlo_mpi",
                    help="Path to MPI executable (relative to benchmarking/). Default: ../Monte_carlo_mpi")
    ap.add_argument("--outdir", default="results", help="Results directory (created if missing).")
    ap.add_argument("--densities", nargs="*", type=float, default=DENSITIES)
    ap.add_argument("--procs", nargs="*", type=int, default=PROCS)
    ap.add_argument("--trials", type=int, default=TRIALS)
    ap.add_argument("--nsteps", type=int, default=NSTEPS)
    ap.add_argument("--rcut", type=float, default=RCUT)
    ap.add_argument("--delta", type=float, default=DELTA)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--extra-mpirun", nargs="*", default=[],
                    help="Extra mpirun args (e.g. --bind-to core --map-by socket)")
    args = ap.parse_args()

    exe = Path(args.exe)
    if not exe.exists():
        print(f"ERROR: Executable not found: {exe.resolve()}", file=sys.stderr)
        sys.exit(1)

    outdir = Path(args.outdir); ensure_dir(outdir)
    inputs_dir = outdir / "inputs"; ensure_dir(inputs_dir)
    logs_dir = outdir / "logs"; ensure_dir(logs_dir)

    # Build (rho, N, Lx, Ly) plan with constant density
    plan = []
    for rho in args.densities:
        Ns = N_BY_RHO.get(rho, [])
        if not Ns:
            print(f"[warn] No N list for rho={rho}, skipping.", file=sys.stderr)
            continue
        for N in Ns:
            L = math.sqrt(N / rho)      # Lx=Ly to preserve rho
            plan.append((rho, N, L, L))

    csv_path = outdir / "mpi_bench.csv"
    with open(csv_path, "w", newline="") as fcsv:
        wr = csv.writer(fcsv)
        wr.writerow(["density","N","Lx","Ly","np","trial","seconds","retcode","log_path"])

        for (rho, N, Lx, Ly) in plan:
            # Per-(rho,N) input
            tag_base = f"rho{rho:.2f}_N{N}"
            input_text = make_input_text(
                Lx=Lx, Ly=Ly, N=N, rcut=args.rcut, T=TEMP,
                nSteps=args.nsteps, eSteps=ESTEPS,
                outputFreq=OUTPUT_FREQ, cellUpdateFreq=CELL_UPDATE_FREQ,
                potential=POTENTIAL,
                # out files remain in benchmarking/ results/ via working dir
                out_xyz=f"{tag_base}.xyz", out_data=f"{tag_base}.dat",
                delta=args.delta, seed=args.seed, ensemble=ENSEMBLE
            )
            in_path = inputs_dir / f"{tag_base}.inp"
            in_path.write_text(input_text)

            for np in args.procs:
                for trial in range(1, args.trials+1):
                    tag = f"{tag_base}_np{np}_t{trial}"
                    print(f"[run] {tag}")
                    sec, ret, out = run_once(str(exe), np, in_path, args.extra_mpirun)
                    logp = logs_dir / f"{tag}.log"
                    logp.write_text(out)
                    wr.writerow([rho, N, f"{Lx:.8f}", f"{Ly:.8f}", np, trial, f"{sec:.6f}", ret, str(logp)])
                    fcsv.flush()

    # ---- Plotting ----
    df = pd.read_csv(csv_path)

    # Median over trials
    gcols = ["density","N","Lx","Ly","np"]
    agg = df.groupby(gcols)["seconds"].median().reset_index()

    # 1) Runtime vs N for each density (log-log)
    for rho, sub in agg.groupby("density"):
        fig, ax = plt.subplots(figsize=(6,4), dpi=130)
        for np, subnp in sub.groupby("np"):
            ax.plot(subnp["N"], subnp["seconds"], marker="o", label=f"np={np}")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(f"MPI Runtime vs N  (rho={rho:.2f})")
        ax.set_xlabel("N (particles)")
        ax.set_ylabel("Median runtime (s)")
        ax.grid(True, which="both", ls="--", alpha=0.4)
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / f"runtime_vs_N_rho{rho:.2f}.png")

    # 2) Speedup vs np at largest N for each density
    for rho, sub in agg.groupby("density"):
        subN = sub[sub["N"] == sub["N"].max()]
        if subN.empty: continue
        if 1 not in set(subN["np"]):   # need np=1 for baseline
            continue
        T1 = float(subN[subN["np"]==1]["seconds"].iloc[0])
        subN = subN.sort_values("np").copy()
        subN["speedup"] = T1 / subN["seconds"]

        fig, ax = plt.subplots(figsize=(6,4), dpi=130)
        ax.plot(subN["np"], subN["speedup"], marker="o", label="measured")
        ax.plot(subN["np"], subN["np"], ls="--", color="gray", label="ideal")
        ax.set_title(f"Speedup vs np (rho={rho:.2f}, N={int(subN['N'].iloc[0])})")
        ax.set_xlabel("MPI ranks (np)")
        ax.set_ylabel("Speedup S = T1/Tp")
        ax.grid(True, ls="--", alpha=0.4)
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / f"speedup_rho{rho:.2f}_N{int(subN['N'].iloc[0])}.png")

    # 3) Efficiency (optional)
    for rho, sub in agg.groupby("density"):
        subN = sub[sub["N"] == sub["N"].max()]
        if subN.empty or 1 not in set(subN["np"]): continue
        T1 = float(subN[subN["np"]==1]["seconds"].iloc[0])
        subN = subN.sort_values("np").copy()
        subN["speedup"] = T1 / subN["seconds"]
        subN["eff"] = 100.0 * subN["speedup"] / subN["np"]

        fig, ax = plt.subplots(figsize=(6,4), dpi=130)
        ax.plot(subN["np"], subN["eff"], marker="s", color="tab:green")
        ax.set_title(f"Efficiency (rho={rho:.2f}, N={int(subN['N'].iloc[0])})")
        ax.set_xlabel("MPI ranks (np)")
        ax.set_ylabel("Efficiency (%)")
        ax.set_ylim(0, 110)
        ax.grid(True, ls="--", alpha=0.4)
        fig.tight_layout()
        fig.savefig(outdir / f"efficiency_rho{rho:.2f}_N{int(subN['N'].iloc[0])}.png")

    print("\nDone.")
    print(f"CSV: {csv_path}")
    print(f"Plots: {outdir}/*.png")
    print(f"Logs:  {logs_dir}/*.log")

if __name__ == "__main__":
    main()

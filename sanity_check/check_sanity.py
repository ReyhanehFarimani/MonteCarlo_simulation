# sanity_check/check_sanity.py
import csv
import math
import subprocess
import sys
from pathlib import Path

PARALLEL_EXE = "./run_parallel"   # adjust if needed
SERIAL_EXE   = "./run_serial"

def run(cmd):
    print(">>", " ".join(cmd))
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(proc.stdout)
    if proc.returncode != 0:
        raise SystemExit(f"Command failed: {' '.join(cmd)}")

def read_csv(path):
    rows = []
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append({
                "frame":    int(row["frame"]),
                "N":        int(row.get("N_global", row.get("N", "0"))),
                "rho":      float(row["rho"]),
                "Energy":   float(row["Energy"]),
                "Virial":   float(row["Virial"]),
                "Pressure": float(row["Pressure"]),
            })
    rows.sort(key=lambda x: x["frame"])
    return rows

def approx_equal(a, b, atol, rtol):
    if not (math.isfinite(a) and math.isfinite(b)):
        return a == b
    return abs(a - b) <= max(atol, rtol * max(abs(a), abs(b)))

def main():
    # 1) Run the parallel program (modify -np if desired)
    if not Path(PARALLEL_EXE).exists():
        print(f"Missing {PARALLEL_EXE}. Build it first.", file=sys.stderr)
        sys.exit(2)
    run(["mpirun", "-np", "10", PARALLEL_EXE])

    # 2) Run the serial program
    if not Path(SERIAL_EXE).exists():
        print(f"Missing {SERIAL_EXE}. Build it first.", file=sys.stderr)
        sys.exit(2)
    run([SERIAL_EXE])

    # 3) Load CSVs
    par = read_csv("sanity_parallel_data.csv")
    ser = read_csv("sanity_serial_data.csv")

    if len(par) != len(ser):
        raise SystemExit(f"Frame count mismatch: parallel={len(par)} serial={len(ser)}")

    # Field-wise tolerances
    # Keep rho quite strict; allow Energy/Virial/Pressure some breathing room.
    TOLS = {
        "rho":      dict(atol=1e-12, rtol=1e-12),
        "Energy":   dict(atol=1e-6 ,  rtol=1e-6),
        "Virial":   dict(atol=1e-6,  rtol=1e-6),
        "Pressure": dict(atol=1e-6,  rtol=1e-6),
    }

    # 4) Compare per frame
    fails = 0
    for p, s in zip(par, ser):
        ok = True
        msgs = []
        if p["frame"] != s["frame"]:
            ok = False; msgs.append(f"frame id {p['frame']} vs {s['frame']}")
        if p["N"]    != s["N"]:
            ok = False; msgs.append(f"N {p['N']} vs {s['N']}")

        for k in ("rho", "Energy", "Virial", "Pressure"):
            tol = TOLS[k]
            if not approx_equal(p[k], s[k], tol["atol"], tol["rtol"]):
                ok = False; msgs.append(f"{k}: {p[k]} vs {s[k]} (atol={tol['atol']}, rtol={tol['rtol']})")

        if not ok:
            print(f"[FAIL] frame {p['frame']}: " + "; ".join(msgs))
            fails += 1

    if fails == 0:
        print("✅ Sanity check passed: serial == parallel within tolerances for all frames.")
        sys.exit(0)
    else:
        print(f"❌ {fails} frame(s) mismatched.")
        sys.exit(1)

if __name__ == "__main__":
    main()

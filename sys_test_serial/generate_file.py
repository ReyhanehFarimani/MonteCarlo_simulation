from pathlib import Path

def generate_input_file(
    Lx=50.0,
    Ly=50.0,
    N=50,
    rcut=2.5,
    T=0.45,
    nSteps=1000000,
    eSteps=60000,
    outputFreq=10000,
    cellUpdateFreq=50,
    mu=0.024,
    seed=42,
    potential="LennardJones",
    ensemble="GCMC",
    position_file=None,
    output_dir="."
):
    # Construct filenames based on parameters
    mu_tag = f"{mu:4f}".replace('.', '_')
    T_tag = f"{T:.2f}".replace('.', '_')
    out_xyz = f"traj_z_{mu_tag}_T_{T_tag}.dump"
    out_data = f"data_z_{mu_tag}_T_{T_tag}.txt"
    filename = f"input_LJ_T_{T_tag}_z_{mu_tag}.txt"

    lines = []
    lines.append("# Required constants")
    lines += [f"{k} {v}" for k, v in {
        "Lx": Lx, "Ly": Ly, "N": N, "rcut": rcut, "T": T,
        "nSteps": nSteps, "eSteps": eSteps,
        "outputFreq": outputFreq, "cellUpdateFreq": cellUpdateFreq
    }.items()]

    lines.append("\n# Optional constants")
    lines.append(f"mu {mu}")
    lines.append(f"seed {seed}")

    lines.append("\n# Filenames (must be valid paths)")
    lines.append(f"potential {potential}")
    lines.append(f"out_xyz {out_xyz}")
    lines.append(f"out_data {out_data}")

    lines.append("\n# Optional: for reading particle positions or choosing ensemble")
    if position_file:
        lines.append(f"position_file {position_file}")
    lines.append(f"ensemble {ensemble}")

    output_path = Path(output_dir) / filename
    output_path.write_text("\n".join(lines))
    return output_path.name


for mu in [0.017, 0.018, 0.0182, 0.0184, 0.0186, 0.0188, 0.019, 0.0192, 0.0194, 0.0196]:
    generate_input_file(mu = mu, T = 0.42)
    

for mu in [0.0241, 0.0242, 0.0243, 0.0244, 0.0245, 0.0246, 0.0247, 0.0248, 0.0249]:
    generate_input_file(mu = mu, T = 0.45)


for mu in [0.0297, 0.03, 0.0305, 0.031, 0.0315, 0.032, 0.0325 , 0.033, 0.0335, 0.034, 0.0345, 0.035, 0.0355, 0.036, 0.0365]:
    generate_input_file(mu = mu, T = 0.47)

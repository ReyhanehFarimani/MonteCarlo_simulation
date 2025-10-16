#!/bin/bash
set -e
EXE=../build/MC_mpi    # adjust to your compiled binary
INPUT=../inputs/test_input.inp
OUTDIR=results
mkdir -p $OUTDIR

for p in 1 2 4 8; do
  for trial in 1 2 3; do
    mpirun -np $p $EXE $INPUT > $OUTDIR/run_${p}r_${trial}.log
    grep "timestep" $OUTDIR/run_${p}r_${trial}.log \
      | awk '{print $2}' >> $OUTDIR/time_${p}.dat
  done
done

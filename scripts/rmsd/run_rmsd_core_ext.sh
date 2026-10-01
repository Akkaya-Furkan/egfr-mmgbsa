#!/bin/bash
# WT extended: production (0-10 ns) + ext1 (10-20 ns), same core-CA fit as run_rmsd_core.sh
set -uo pipefail
CPPTRAJ=/home/furkanakkaya/.conda/envs/md_env/bin/cpptraj
P02=/mnt/HHD/portfolio-projects/P02-MMGBSA
cd "$(dirname "$0")"
for s in 1 2 3; do
  { echo "parm $P02/system_WT.prmtop"
    echo "trajin $P02/production_WT_seed${s}.dcd"
    echo "trajin $P02/production_WT_seed${s}_ext1.dcd"
    echo "autoimage"
    echo "rms core :8-296@CA first out rmsdcoreext_WT_seed${s}.dat"
    echo "rms lig :MOL&!@/H first nofit out rmsdcoreext_WT_seed${s}.dat"
    echo "rms tail :297-327@CA first nofit out rmsdcoreext_WT_seed${s}.dat"
    echo "run"; } > rmsdcoreext_WT_seed${s}.in
  $CPPTRAJ -i rmsdcoreext_WT_seed${s}.in > rmsdcoreext_WT_seed${s}.log 2>&1 &
done
wait
echo finished

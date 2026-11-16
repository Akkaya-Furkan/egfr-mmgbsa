#!/bin/bash
# P02 core-only RMSD: fit on CA of the common core (EGFR 701-989) so the WT/C797S terminal tails
# (694-700, 990-1020/1018) do not distort fit or RMSD. Ligand heavy-atom RMSD after core fit (nofit).
# Tail CA RMSD (nofit, after core fit) is reported separately for WT and C797S_triple.
set -uo pipefail
CPPTRAJ=${CPPTRAJ:-cpptraj}
P02=/mnt/HHD/portfolio-projects/P02-MMGBSA
cd "$(dirname "$0")"
declare -A CORE=( [WT]="8-296" [T790M]="1-289" [C797S_triple]="6-294" )
declare -A TAIL=( [WT]="297-327" [C797S_triple]="295-323" )
for g in WT T790M C797S_triple; do
  for s in 1 2 3; do
    { echo "parm $P02/system_${g}.prmtop"
      echo "trajin $P02/production_${g}_seed${s}.dcd"
      echo "autoimage"
      echo "rms core :${CORE[$g]}@CA first out rmsdcore_${g}_seed${s}.dat"
      echo "rms lig :MOL&!@/H first nofit out rmsdcore_${g}_seed${s}.dat"
      [ -n "${TAIL[$g]:-}" ] && echo "rms tail :${TAIL[$g]}@CA first nofit out rmsdcore_${g}_seed${s}.dat"
      echo "run"; } > rmsdcore_${g}_seed${s}.in
    $CPPTRAJ -i rmsdcore_${g}_seed${s}.in > rmsdcore_${g}_seed${s}.log 2>&1 &
  done
done
wait
ls rmsdcore_*.dat | wc -l

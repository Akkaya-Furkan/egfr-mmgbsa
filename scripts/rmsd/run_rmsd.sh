#!/bin/bash
# P02 RMSD platosu: protein CA (ilk kareye hizalı) ve ligand ağır atom RMSD (protein hizasında, fit yok).
# 9 trajektori, 500 kare, 20 ps aralık. MM-GBSA penceresi: kare 251-500 (son 5 ns).
set -euo pipefail
CPPTRAJ=${CPPTRAJ:-cpptraj}
P02=/mnt/HHD/portfolio-projects/P02-MMGBSA
cd "$(dirname "$0")"
for g in WT T790M C797S_triple; do
  for s in 1 2 3; do
    cat > "rmsd_${g}_seed${s}.in" <<IN
parm $P02/system_${g}.prmtop
trajin $P02/production_${g}_seed${s}.dcd
autoimage
rms prot !:MOL&@CA first out rmsd_${g}_seed${s}.dat
rms lig :MOL&!@/H first nofit out rmsd_${g}_seed${s}.dat
run
IN
    "$CPPTRAJ" -i "rmsd_${g}_seed${s}.in" > "rmsd_${g}_seed${s}.log" 2>&1 &
  done
done
wait
grep -L "Run: [0-9]\+\|TIME: Total execution" rmsd_*.log || true

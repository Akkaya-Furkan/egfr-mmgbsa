#!/bin/bash
# Dry complex/receptor/ligand topologies with ante-MMPBSA.py (AmberTools 24.8, env md_env).
# Ligand residue is selected by NAME (:MOL): its index differs per system (WT 328, T790M 290, C797S_triple 324).
set -e
P02=/mnt/HHD/portfolio-projects/P02-MMGBSA
# (a) first pass: full constructs, solvent/ions stripped only            -> $P02/mmgbsa/
mkdir -p $P02/mmgbsa && cd $P02/mmgbsa
for g in WT T790M C797S_triple; do
  ante-MMPBSA.py -p ../system_${g}.prmtop -c complex_${g}.prmtop -r receptor_${g}.prmtop -l ligand_${g}.prmtop \
     -s ":WAT,Na+,Cl-" -n ":MOL" --radii mbondi2
done
# (b) common receptor range EGFR 701-989 (terminal tails removed)       -> $P02/mmgbsa_common/
#     WT      tails: residues 1-7 (694-700) and 297-327 (990-1020)
#     triple  tails: residues 1-5 (696-700) and 295-323 (990-1018)
#     T790M   is 701-989 already
mkdir -p $P02/mmgbsa_common && cd $P02/mmgbsa_common
declare -A MASK=( [WT]=":WAT,Na+,Cl-,1-7,297-327" [C797S_triple]=":WAT,Na+,Cl-,1-5,295-323" [T790M]=":WAT,Na+,Cl-" )
for g in WT T790M C797S_triple; do
  ante-MMPBSA.py -p ../system_${g}.prmtop -c complex_${g}.prmtop -r receptor_${g}.prmtop -l ligand_${g}.prmtop \
     -s "${MASK[$g]}" -n ":MOL" --radii mbondi2
done
# check: receptor net charge WT +2, T790M +2, triple +4 (full constructs: -8, +2, -6)

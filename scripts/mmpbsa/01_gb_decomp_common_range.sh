#!/bin/bash
# MM-GBSA (igb=5, 0.15 M) over the last 5 ns (frames 251-500, 250 frames) and per-residue decomposition
# (idecomp=2, every 10th frame = 25 frames), common receptor range. One MMPBSA.py process per run (no mpirun on the host).
source /etc/profile.d/conda.sh && conda activate md_env
P02=/mnt/HHD/portfolio-projects/P02-MMGBSA; cd $P02/mmgbsa_common
declare -A MASK=( [WT]=":WAT,Na+,Cl-,1-7,297-327" [C797S_triple]=":WAT,Na+,Cl-,1-5,295-323" [T790M]=":WAT,Na+,Cl-" )
for g in WT T790M C797S_triple; do
cat > gb_${g}.in << IN
GB common range 701-989, last 5 ns
&general
  startframe=251, endframe=500, interval=1, verbose=1, keep_files=0,
  strip_mask="${MASK[$g]}",
/
&gb
  igb=5, saltcon=0.150,
/
IN
cat > decomp_${g}.in << IN
Decomp common range 701-989 (+ligand = residue 290), sparse sampling
&general
  startframe=251, endframe=500, interval=10, verbose=1, keep_files=0,
  strip_mask="${MASK[$g]}",
/
&gb
  igb=5, saltcon=0.150,
/
&decomp
  idecomp=2, dec_verbose=1,
  print_res="1-290",
/
IN
done
for g in WT T790M C797S_triple; do for s in 1 2 3; do
  d=run_GB_${g}_seed${s}; mkdir -p $d
  (cd $d && MMPBSA.py -O -i ../gb_${g}.in -o FINAL_GB_${g}_seed${s}.dat -sp ../../system_${g}.prmtop \
     -cp ../complex_${g}.prmtop -rp ../receptor_${g}.prmtop -lp ../ligand_${g}.prmtop \
     -y ../../production_${g}_seed${s}.dcd > run.log 2>&1) &
done; done; wait
# decomposition was run for T790M seeds 2-3 and the triple mutant seeds 1-3 here;
# T790M seed 1 uses the identical input from the first pass (T790M needs no truncation).
for gs in "T790M 2" "T790M 3" "C797S_triple 1" "C797S_triple 2" "C797S_triple 3"; do set -- $gs; g=$1; s=$2
  d=run_DEC_${g}_seed${s}; mkdir -p $d
  (cd $d && MMPBSA.py -O -i ../decomp_${g}.in -o FINAL_DECOMP_${g}_seed${s}.dat -do FINAL_DECOMP_full_${g}_seed${s}.dat \
     -sp ../../system_${g}.prmtop -cp ../complex_${g}.prmtop -rp ../receptor_${g}.prmtop -lp ../ligand_${g}.prmtop \
     -y ../../production_${g}_seed${s}.dcd > run.log 2>&1) &
done; wait

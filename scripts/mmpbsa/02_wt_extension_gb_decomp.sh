#!/bin/bash
# WT after the 10 ns extension (extend_wt.py segment 1): same protocol on the extension trajectory,
# frames 251-500 of production_WT_seed{n}_ext1.dcd (= 15-20 ns of the 20 ns WT run).
source /etc/profile.d/conda.sh && conda activate md_env
P02=/mnt/HHD/portfolio-projects/P02-MMGBSA; cd $P02/mmgbsa_common
MASK=":WAT,Na+,Cl-,1-7,297-327"
cat > gb_WText.in << IN
GB common range 701-989, WT extension window (15-20 ns = ext1 frames 251-500)
&general
  startframe=251, endframe=500, interval=1, verbose=1, keep_files=0,
  strip_mask="${MASK}",
/
&gb
  igb=5, saltcon=0.150,
/
IN
cat > decomp_WText.in << IN
Decomp common range, WT extension window, sparse sampling (25 frames)
&general
  startframe=251, endframe=500, interval=10, verbose=1, keep_files=0,
  strip_mask="${MASK}",
/
&gb
  igb=5, saltcon=0.150,
/
&decomp
  idecomp=2, dec_verbose=1,
  print_res="1-290",
/
IN
for s in 1 2 3; do
  d=run_GBext_WT_seed${s}; mkdir -p $d
  (cd $d && MMPBSA.py -O -i ../gb_WText.in -o FINAL_GB_WText_seed${s}.dat -sp ../../system_WT.prmtop \
     -cp ../complex_WT.prmtop -rp ../receptor_WT.prmtop -lp ../ligand_WT.prmtop \
     -y ../../production_WT_seed${s}_ext1.dcd > run.log 2>&1) &
  d=run_DECext_WT_seed${s}; mkdir -p $d
  (cd $d && MMPBSA.py -O -i ../decomp_WText.in -o FINAL_DECOMP_WText_seed${s}.dat -do FINAL_DECOMP_full_WText_seed${s}.dat \
     -sp ../../system_WT.prmtop -cp ../complex_WT.prmtop -rp ../receptor_WT.prmtop -lp ../ligand_WT.prmtop \
     -y ../../production_WT_seed${s}_ext1.dcd > run.log 2>&1) &
done; wait

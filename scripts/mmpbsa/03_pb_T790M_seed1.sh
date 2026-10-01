#!/bin/bash
# MM-PBSA cross-check on one complex (T790M seed 1, last 5 ns, 250 frames; ~90 min single core).
source /etc/profile.d/conda.sh && conda activate md_env
cd /mnt/HHD/portfolio-projects/P02-MMGBSA/mmgbsa; mkdir -p run_PB_T790M_seed1
cat > mmpbsa_pb.in << IN
PB free energy, last 5 ns, T790M seed1 comparison
&general
  startframe=251, endframe=500, interval=1, verbose=2, keep_files=0,
/
&pb
  inp=2, radiopt=0, istrng=0.150,
/
IN
cd run_PB_T790M_seed1 && MMPBSA.py -O -i ../mmpbsa_pb.in -o FINAL_PB_T790M_seed1.dat -sp ../../system_T790M.prmtop \
  -cp ../complex_T790M.prmtop -rp ../receptor_T790M.prmtop -lp ../ligand_T790M.prmtop -y ../../production_T790M_seed1.dcd

## Does end-point free energy rank EGFR drug resistance better than docking? (osimertinib, WT / T790M / L858R-T790M-C797S)

Took the three most discriminating EGFR–osimertinib complexes from a docking study and re-scored them with
MM-GBSA on explicit-solvent MD (OpenMM on an AMD GPU, AmberTools; 3 independent-seed replicates per genotype,
WT extended to 20 ns), then compared the ranking with docking and decomposed it per residue.

**What came out.** Both methods rank T790M first. MM-GBSA separates T790M from the triple-mutant construct by
15 kcal/mol (95 % CI 6–24, n = 3 seeds), two thirds of it van der Waals; the gap to WT (12 kcal/mol) is not
statistically resolved because osimertinib drifts away from its docked pose in WT, giving ΔG between −27 and −56
kcal/mol across seeds. Per-residue decomposition puts the T790M gain on residue 718 (glycine-rich loop) and the hinge residues 792–796,
not on residue 790 itself.

**Things I checked that changed the answer.** (1) The three systems had different receptor lengths and net charges
(−8, +2, −6); individual electrostatic/GB terms were not comparable, so all energies were recomputed on a common
701–989 receptor (totals moved by < 1 kcal/mol). (2) Whole-protein Cα RMSD made WT look unequilibrated because of
floppy terminal tails; core-only RMSD and ligand contacts show the real issue is ligand pose drift. (3) GB and PB
disagree by ≈ 60 kcal/mol in absolute ΔG on the same frames, so only rankings are interpreted. (4) A first
single-seed residue decomposition over-stated residue 797; with three seeds it is not resolved.

**Stack:** OpenMM (HIP), AmberTools (tleap, antechamber, MMPBSA.py, cpptraj), ff14SB/GAFF2, Python (pandas, SciPy, matplotlib).

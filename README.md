# EGFR resistance — does MM-GBSA rank better than docking?

**Question.** Does end-point free energy (MM-GBSA) rank the EGFR resistance pattern for osimertinib better
than docking scores ([egfr-resistance-docking](https://github.com/Akkaya-Furkan/egfr-resistance-docking))?

**Status.** Complete (October 2026). Trajectories (7.5 GB) stay on the compute host and are not in this repository;
all derived data needed to regenerate every number and figure are.

## Result in brief

| Genotype (start structure) | Docking, mean of 3 seeds (kcal/mol) | ΔG_GBSA, mean ± SD of 3 MD seeds (kcal/mol) | MM-GBSA window |
|---|---|---|---|
| WT (4WKQ) | −7.87 | **−43.5 ± 15.1** (seeds: −26.8 / −47.8 / −56.1) | 15–20 ns (after extension) |
| T790M (6JX0) | −8.80 | **−55.8 ± 4.2** | 5–10 ns |
| L858R/T790M/C797S\* (6LUD) | −8.67 | **−40.8 ± 3.2** | 5–10 ns |

\* The 6LUD construct also carries E865A/E866A/K867A; differences from the other two are therefore not attributable to L858R and C797S alone.

* **Both methods put T790M first.** MM-GBSA resolves T790M from the triple-mutant construct: −15.0 kcal/mol
  (95 % CI −23.8 … −6.3; Welch p = 0.01; n = 3 seeds per genotype). About two thirds of that gap (−10.1) is van der Waals.
* **T790M vs WT is not resolved** (−12.2; 95 % CI −46.6 … +22.1; p = 0.29). WT ΔG spans −27 to −56 kcal/mol across seeds because osimertinib leaves its starting pose in WT (see below).
* **WT vs triple mutant is not resolved** (−2.8 ± 8.9 SE). Docking ranks the triple mutant above WT; MM-GBSA cannot say. The MM-GBSA order (T790M < WT < triple) matches the clinical pattern (osimertinib active on T790M, resistance with C797S) but only the T790M end of it is statistically supported.
* **Per-residue decomposition** (3 seeds vs 3 seeds, 25 frames each): the T790M gain concentrates on residues 718 (glycine-rich loop; −2.5 ± 0.5 SE), 796 (−1.8 ± 0.5), 793 (−1.3 ± 0.3), 792 (hinge) and 723, 726; residue 790 itself contributes −0.1 ± 0.3. Residue 797 (−1.66 ± 0.85) is not resolved.
* **Absolute values are model-dependent.** For T790M seed 1 the same frames give ΔG_GB = −60.7 and ΔG_PB = +0.8 ± 0.3 (SEM) kcal/mol with identical gas-phase energy (−102.2); only the ranking is interpretable.

Figures: `results/fig1_energy_components.png`, `results/fig2_residue_contribution.png`, `results/fig3_rmsd_core_contacts.png`.

## Methods

* **Systems.** Osimertinib (Vina poses from the docking project) in WT (4WKQ, 30 internal missing residues rebuilt), T790M (6JX0) and the L858R/T790M/C797S construct (6LUD). ff14SB protein, GAFF2/AM1-BCC ligand (one shared charge set), TIP3P, 0.15 M NaCl, PME 1.0 nm, H-bond constraints, 2 fs.
* **MD.** OpenMM 8.6 (HIP, AMD gfx1200, ≈ 210 ns/day). Minimisation → NVT/NPT → 3 × 10 ns per genotype (300 K, 1 bar, Langevin-middle, MC barostat); the 3 seeds re-draw velocities from **one** equilibrated state, so they are not independent starting structures. WT was extended by a further 10 ns per seed (positions and box from the last production frame, velocities re-drawn).
* **MM-GBSA.** AmberTools 24.8 `MMPBSA.py` 14.0, single-trajectory, igb = 5, 0.15 M salt, 250 frames from the last 5 ns of the chosen window. Receptor restricted to the **common range EGFR 701–989** for all genotypes (see next section). Residue decomposition: `idecomp = 2`, 25 frames. MM-PBSA on one complex. No entropy term.
* **Statistics.** Welch t-test and 95 % CI on seed means (n = 3 vs 3), no multiplicity correction (three pairwise comparisons; the T790M–triple p = 0.01 stays below 0.05 after Bonferroni).

## Construct length and net charge

The three systems were built with different receptor lengths: WT 694–1020, T790M 701–989, triple 696–1018. The C-terminal tail (990–1020) carries −10, so the receptor net charge was −8 (WT), +2 (T790M) and −6 (triple). Individual EEL and EGB terms were therefore not comparable (T790M EEL −17.7 vs WT −129.1 kcal/mol, cancelling in EGB), and only their sum is reported. All free energies here use residues 701–989 only (net charge +2, +2, +4; the +2 of the triple construct is L858R plus E865A/E866A/K867A). Effect on totals: WT 5–10 ns −41.7 → −41.1, triple −41.3 → −40.8, T790M unchanged (regression check, difference 0.0) — `results/P02_full_construct_sensitivity.csv`.

## Equilibration and pose stability

RMSD is computed on core Cα (701–989) so the floating tails do not distort fit or values (tail Cα RMSD is 8.5–11.9 Å in WT and 5.6–8.9 Å in the triple construct). Native contacts are ligand–protein heavy-atom pairs within 4.5 Å, as a fraction of those present in the first frame. Plateau thresholds are the largest |slope| in the last 5 ns of the T790M and triple runs (core 0.032, ligand 0.076 Å/ns).

| Window | Core Cα RMSD (Å), seeds 1/2/3 | Ligand RMSD (Å) | Native contacts kept | Plateau (core/ligand) |
|---|---|---|---|---|
| WT 5–10 ns | 2.2 / 1.5 / 2.5 | 4.3 / 1.9 / 3.5 | 0.20 / 0.39 / 0.20 | no/no · yes/yes · no/no |
| WT 15–20 ns | 3.0 / 1.6 / 2.5 | 5.6 / 2.1 / 3.6 | 0.11 / 0.42 / 0.19 | yes/yes · yes/yes · no/no |
| T790M 5–10 ns | 1.5 / 1.8 / 1.9 | 1.1 / 1.6 / 1.8 | 0.61 / 0.54 / 0.52 | all yes |
| Triple 5–10 ns | 1.5 / 1.2 / 1.2 | 1.8 / 1.6 / 1.7 | 0.48 / 0.52 / 0.48 | all yes |

In WT, osimertinib relaxes away from its docked pose in all three seeds (seed 1 to 5.6 Å, seed 3 not yet plateaued at +0.22 Å/ns), whereas it stays near the starting pose in T790M and the triple construct. The WT ΔG therefore reflects different poses (−27 to −56 kcal/mol) and is **not converged**.

## Limitations

1. **Starting-structure bias.** WT starts from a gefitinib-bound crystal, T790M and the triple construct from structures with the co-crystal ligand YY3; osimertinib was docked into all three. Pocket shape differs by construction, which can favour the mutants and contributes to the WT pose drift.
2. **n = 3 seeds from one equilibrated state** per genotype: SDs underestimate structural uncertainty; the T790M–WT comparison is under-powered.
3. **WT not converged** at 20 ns (seed 3 drifting, seed 1 settled in a displaced pose). A longer WT run from several starting poses is the remaining open item.
4. Osimertinib binds C797 covalently; everything here is non-covalent, so C797S resistance is not expected to appear in these energies.
5. No entropy; GB/PB disagree by ≈ 60 kcal/mol in absolute terms; only rankings are used.
6. Per-residue decomposition uses 25 frames per run and three seeds; an earlier single-seed version (one seed vs an unequilibrated WT) suggested residue 797 as the largest contributor (−3.4); with three seeds it is −1.66 ± 0.85 and not resolved, and residue 718 leads.

## Reproduce

```bash
python analysis/make_analysis.py   # regenerates results/*.csv and results/fig1-3 from data/
```
`scripts/md/` builds, minimises, equilibrates, runs and extends the MD (GPU host); `scripts/mmpbsa/` holds the exact `ante-MMPBSA.py` / `MMPBSA.py` commands and inputs (`inputs_first_pass/` are the full-construct inputs used before the common-range re-analysis); `scripts/rmsd/` the cpptraj RMSD/contact analyses (`run_rmsd.sh` is the all-Cα version, kept for the record). Raw MMPBSA outputs are in `data/mmpbsa/`, RMSD/contact series in `data/rmsd/`, docking scores in `data/docking/`.

## Layout

| Path | Content |
|---|---|
| `results/P02_energy_table.csv` | ΔG and components per genotype (mean ± SD) |
| `results/P02_ddG_components.csv` | Welch differences with 95 % CI for ΔG, vdW, EEL+EGB, ESURF |
| `results/P02_docking_vs_gbsa_ranking.csv` | docking vs MM-GBSA ranks |
| `results/P02_residue_delta.csv` | per-residue means, SDs and T790M − WT / T790M − triple differences |
| `results/P02_rmsd_core_stats.csv` | window means, slopes, plateau flags, native contacts |
| `results/P02_gb_runs_all.csv`, `P02_full_construct_sensitivity.csv` | per-run values; full-construct vs common-range |
| `UPWORK.md` | write-up for the portfolio |

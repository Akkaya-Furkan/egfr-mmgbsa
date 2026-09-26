# EGFR resistance — does MM-GBSA rank better than docking?

**Question.** Does end-point free energy (MM-GBSA) rank the EGFR resistance pattern better than docking scores?

## Data
The three most discriminating complexes from egfr-resistance-docking.

## Approach
1. GAFF2 ligands (acpype), amber99sb-ildn protein, TIP3P water, 0.15 M NaCl
2. Minimisation → NVT/NPT equilibration → 3 × 10 ns production (GROMACS)
3. MM-GBSA on the last 5 ns (gmx_MMPBSA); MM-PBSA cross-check on one complex
4. Per-residue decomposition: where does the T790M effect come from?
5. Docking rank vs ΔG_GBSA rank

## Deliverables
Energy table, component bar chart, per-residue contribution map, write-up.

## Limitations
- Three systems give a trend, not a statistically robust correlation.
- MM-GBSA usually omits entropy; values are relative rankings, not absolute ΔG.

## Status
Planned — part of a 30-project computational biology and scientific software portfolio.

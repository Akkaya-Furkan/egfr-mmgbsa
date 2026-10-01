"""P02 - build a full-hydrogen, acpype-atom-named ligand PDB for each docked complex.

Vina's docked pose only carries heavy atoms (+ a couple of polar H's); MD needs all
34 hydrogens with names matching the antechamber/acpype-derived GAFF2 parameters
(osimertinib_AC.lib/frcmod), so that tleap can build protein+ligand+water using one
shared, once-derived charge set (not re-fit per system, which would confound the
receptor-driven energy comparison the project is testing).

Approach: acpype's charge derivation preserved the original SDF's coordinates
exactly (checked: max deviation 0.0008 A), so heavy_atom_template_index -> acpype
atom name is a fixed identity map, established once from the original conformer,
and valid for any new conformation built on the same atom-indexed RDKit mol (i.e.
every docked pose, which is the same molecule, just re-oriented/re-torsioned).
"""
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, "/mnt/HHD/portfolio-projects/egfr-resistance-docking/analysis")
from dock_utils import build_pdbqt_to_template_map, coords_in_template_order, parse_pose_pdbqt

P01 = Path("/mnt/HHD/portfolio-projects/egfr-resistance-docking")
P02 = Path("/mnt/HHD/portfolio-projects/P02-MMGBSA")
ACPYPE = P02 / "osimertinib.acpype"

OSIMERTINIB_SMILES = "C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C"

SYSTEMS = [
    ("WT", 3),
    ("T790M", 1),
    ("C797S_triple", 1),
]


def main():
    # 1) original (protonated, +1) sdf, all-atom
    sdf_full = Chem.MolFromMolFile(str(P01 / "data" / "ligands" / "osimertinib.sdf"), removeHs=False)
    sdf_full_coords = sdf_full.GetConformer().GetPositions()
    charged_atoms = [a.GetIdx() for a in sdf_full.GetAtoms() if a.GetFormalCharge() != 0]
    print("charged atom(s) in original sdf:", charged_atoms,
          [sdf_full.GetAtomWithIdx(i).GetSymbol() for i in charged_atoms])

    # 2) acpype NEW.pdb: name + coord per atom (same coordinates as sdf, just renamed/reordered)
    acpype_names, acpype_coords = [], []
    for line in (ACPYPE / "osimertinib_NEW.pdb").read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            acpype_names.append(line[12:16].strip())
            acpype_coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    acpype_coords = np.array(acpype_coords)

    def nearest_map(src_coords, dst_coords):
        mapping, max_dev = {}, 0.0
        for i, xyz in enumerate(src_coords):
            d = np.linalg.norm(dst_coords - xyz, axis=1)
            j = int(np.argmin(d))
            mapping[i] = j
            max_dev = max(max_dev, float(d[j]))
        return mapping, max_dev

    sdf_full_to_acpype_idx, dev1 = nearest_map(sdf_full_coords, acpype_coords)
    print("sdf_full<->acpype NEW.pdb match max dev:", dev1)

    # 3) heavy-atom template <-> sdf_full correspondence (template built via SMILES, heavy only)
    template = Chem.MolFromSmiles(OSIMERTINIB_SMILES)
    sdf_heavy = Chem.MolFromMolFile(str(P01 / "data" / "ligands" / "osimertinib.sdf"), removeHs=True, sanitize=True)
    match = sdf_heavy.GetSubstructMatch(template)  # position=template_idx, value=sdf_heavy_idx
    assert len(match) == template.GetNumAtoms(), "template match incomplete"
    sdf_heavy_coords = sdf_heavy.GetConformer().GetPositions()
    sdf_heavy_to_full, dev2 = nearest_map(sdf_heavy_coords, sdf_full_coords)
    print("sdf_heavy<->sdf_full match max dev:", dev2)

    template_to_acpype_name = {}
    for t_idx, heavy_idx in enumerate(match):
        full_idx = sdf_heavy_to_full[heavy_idx]
        acpype_idx = sdf_full_to_acpype_idx[full_idx]
        template_to_acpype_name[t_idx] = acpype_names[acpype_idx]
    print("heavy template->acpype name map (first 5):", list(template_to_acpype_name.items())[:5])

    # 4) mark the correct template N with +1 (matches original protonation) so AddHs adds 34 H's
    charge_full_idx = charged_atoms[0]
    charge_template_idx = None
    for t_idx, heavy_idx in enumerate(match):
        if sdf_heavy_to_full[heavy_idx] == charge_full_idx:
            charge_template_idx = t_idx
    assert charge_template_idx is not None
    template.GetAtomWithIdx(charge_template_idx).SetFormalCharge(1)
    Chem.SanitizeMol(template)
    print("protonated template atom idx:", charge_template_idx, "->", template_to_acpype_name[charge_template_idx])

    pos_to_template, dev3 = build_pdbqt_to_template_map(
        template, P01 / "data" / "ligands" / "osimertinib.sdf", P01 / "data" / "ligands" / "osimertinib.pdbqt")
    print("pdbqt<->template match max dev:", dev3)

    # 5) reference all-atom build (heavy template order + AddHs) to learn H-atom naming ONCE
    ref_heavy_coords = np.array([sdf_heavy_coords[i] for i in match])  # template-index order

    def build_all_atom(heavy_coords_template_order):
        mol = Chem.RWMol(template)
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, heavy_coords_template_order[i].tolist())
        mol.RemoveAllConformers()
        mol.AddConformer(conf, assignId=True)
        mol = mol.GetMol()
        Chem.SanitizeMol(mol)
        mol_h = Chem.AddHs(mol, addCoords=True)
        return mol_h

    ref_mol_h = build_all_atom(ref_heavy_coords)
    ref_h_coords = ref_mol_h.GetConformer().GetPositions()
    # map every atom (heavy+H) of ref_mol_h to an acpype name by nearest 3D match
    # (heavy atoms: should be ~exact; H atoms: RDKit's added positions approx match antechamber's)
    ref_to_acpype_idx, dev4 = nearest_map(ref_h_coords, acpype_coords)
    print("ref all-atom (RDKit AddHs) <-> acpype NEW.pdb match max dev (H's expected ~0.1-0.5 A):", dev4)
    rdkit_idx_to_acpype_name = {i: acpype_names[j] for i, j in ref_to_acpype_idx.items()}

    # sanity: heavy atoms in ref_mol_h are the first N atoms (template order preserved by AddHs)
    n_heavy = template.GetNumAtoms()
    mismatches = 0
    for t_idx in range(n_heavy):
        if rdkit_idx_to_acpype_name[t_idx] != template_to_acpype_name[t_idx]:
            mismatches += 1
    print(f"heavy-atom naming cross-check mismatches (should be 0): {mismatches}/{n_heavy}")

    boxes_dummy = None
    for genotype, seed in SYSTEMS:
        pose_path = P01 / "analysis" / "poses_remote" / f"{genotype}_osimertinib_seed{seed}.pdbqt"
        _, pose_heavy = parse_pose_pdbqt(pose_path, model=1)
        docked_heavy_template_order = coords_in_template_order(pose_heavy, pos_to_template, n_heavy)

        mol_h = build_all_atom(docked_heavy_template_order)
        assert mol_h.GetNumAtoms() == len(rdkit_idx_to_acpype_name), "atom count mismatch after AddHs"

        out_pdb = P02 / f"ligand_{genotype}.pdb"
        conf = mol_h.GetConformer()
        lines = []
        for i, atom in enumerate(mol_h.GetAtoms()):
            name = rdkit_idx_to_acpype_name[i]
            x, y, z = conf.GetAtomPosition(i)
            elem = atom.GetSymbol()
            lines.append(
                f"HETATM{i+1:>5} {name:<4} MOL A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2}"
            )
        lines.append("END")
        out_pdb.write_text("\n".join(lines) + "\n")
        print(f"{genotype}: wrote {out_pdb} ({mol_h.GetNumAtoms()} atoms)")


if __name__ == "__main__":
    main()

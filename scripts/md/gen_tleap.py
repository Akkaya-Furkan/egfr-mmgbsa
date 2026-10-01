"""Build each system with tleap in two passes: (1) solvate + neutralize only, read
back the water count and box volume from the tleap log to compute the ion pairs
needed for 0.15 M NaCl; (2) rebuild with that exact ion count. ff14SB substitutes
for the brief's "amber99sb-ildn" - AmberTools on this host only ships ff14SB (the
modern, widely-adopted successor; no ff99SB-ILDN leaprc present here). Same
underlying science (protein backbone/sidechain force field), different specific
parameterization; noted as a deviation from the brief's exact wording.
"""
import re
import subprocess
from pathlib import Path

P02 = Path("/mnt/HHD/portfolio-projects/P02-MMGBSA")
SYSTEMS = ["WT", "T790M", "C797S_triple"]

HEADER = """source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams {acpype_dir}/osimertinib_AC.frcmod
loadoff {acpype_dir}/osimertinib_AC.lib

protein = loadpdb {protein_pdb}
ligand = loadpdb {ligand_pdb}
complex = combine {{protein ligand}}
solvatebox complex TIP3PBOX 10.0 iso
"""

PASS1_TAIL = """charge complex
quit
"""

# addIonsRand's "0 for #ion1" neutralization form forbids a second ion type in the
# SAME call (tleap: "ion2 cannot be specified" when ion1 count is 0) - so
# neutralization (single ion, opposite sign to the net charge) and the 0.15 M
# salt (both ion types, both nonzero) have to be two separate addIonsRand calls.
PASS2_TAIL = """{neutralize_line}
addionsrand complex Na+ {n_pairs} Cl- {n_pairs}
saveamberparm complex {out_prefix}.prmtop {out_prefix}.inpcrd
savepdb complex {out_prefix}_solvated.pdb
quit
"""


def run_tleap(script_text, log_path):
    in_path = log_path.with_suffix(".in")
    in_path.write_text(script_text)
    result = subprocess.run(["tleap", "-f", str(in_path)], capture_output=True, text=True)
    log_path.write_text(result.stdout + "\n---STDERR---\n" + result.stderr)
    return result.stdout


for genotype in SYSTEMS:
    header = HEADER.format(
        acpype_dir=str(P02 / "osimertinib.acpype"),
        protein_pdb=str(P02 / f"protein_{genotype}.pdb"),
        ligand_pdb=str(P02 / f"ligand_{genotype}.pdb"),
    )

    log1 = run_tleap(header + PASS1_TAIL, P02 / f"tleap_{genotype}_pass1.log")
    wat_match = re.search(r"Added (\d+) residues", log1)
    n_wat_added = int(wat_match.group(1)) if wat_match else None
    charge_match = re.search(r"Total (?:perturbed |unperturbed )?charge:\s*(-?\d+\.\d+)", log1)
    net_charge = float(charge_match.group(1)) if charge_match else None
    # 0.15 M salt: n_ion_pairs ~= 0.15 / 55.5 * n_waters
    n_pairs = int(round(0.15 / 55.5 * (n_wat_added or 0)))
    print(f"{genotype}: waters added={n_wat_added}, net charge={net_charge}, NaCl pairs for 0.15 M={n_pairs}")

    if net_charge is None:
        raise SystemExit(f"{genotype}: could not parse net charge from tleap log:\n{log1[-2000:]}")
    n_neutralize = int(round(abs(net_charge)))
    if n_neutralize == 0:
        neutralize_line = "# system already neutral, no counter-ion needed"
    elif net_charge > 0:
        neutralize_line = f"addionsrand complex Cl- 0"
    else:
        neutralize_line = f"addionsrand complex Na+ 0"

    tail2 = PASS2_TAIL.format(n_pairs=n_pairs, neutralize_line=neutralize_line,
                              out_prefix=str(P02 / f"system_{genotype}"))
    log2 = run_tleap(header + tail2, P02 / f"tleap_{genotype}_pass2.log")
    ok = (P02 / f"system_{genotype}.prmtop").exists() and (P02 / f"system_{genotype}.inpcrd").exists()
    print(f"{genotype}: prmtop/inpcrd written = {ok}")

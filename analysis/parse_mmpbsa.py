"""Parsers for MMPBSA.py FINAL_*.dat (GB/PB summary) and idecomp per-residue CSV blocks."""
import re, io
import pandas as pd

def parse_mmpbsa_dat(path):
    text = open(path).read()
    m = re.search(r'\n(GENERALIZED BORN|POISSON BOLTZMANN):\n\n(.*?)(?=\n-{5,}\n-{5,}\n|\Z)', text, re.S)
    body = m.group(2) if m else text
    subpat = re.compile(r'(Complex|Receptor|Ligand|Differences \(Complex - Receptor - Ligand\)):\nEnergy Component.*?\n-+\n(.*?)(?=\n\n\n|\Z)', re.S)
    out = {}
    for sm in subpat.finditer(body):
        rows = {}
        for line in sm.group(2).strip().split("\n"):
            parts = line.split()
            if len(parts) >= 4:
                *name, avg, std, sem = parts
                try:
                    rows[" ".join(name)] = {"avg": float(avg), "std": float(std), "sem": float(sem)}
                except ValueError:
                    pass
        out[sm.group(1)] = rows
    return out

def delta_terms(path):
    d = parse_mmpbsa_dat(path)["Differences (Complex - Receptor - Ligand)"]
    return {"dG": d["DELTA TOTAL"]["avg"], "VDW": d["VDWAALS"]["avg"], "EEL": d["EEL"]["avg"],
            "EGB": d["EGB"]["avg"], "ESURF": d["ESURF"]["avg"], "Ggas": d["DELTA G gas"]["avg"],
            "Gsolv": d["DELTA G solv"]["avg"]}

DECOMP_COLS = ["Residue","Location","Int_avg","Int_std","Int_sem","VDW_avg","VDW_std","VDW_sem",
               "EEL_avg","EEL_std","EEL_sem","PB_avg","PB_std","PB_sem","NP_avg","NP_std","NP_sem",
               "TOTAL_avg","TOTAL_std","TOTAL_sem"]

def parse_decomp_total(path):
    """First (Total Energy Decomposition, DELTAS) block of an idecomp=2 output."""
    lines = open(path).read().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("Residue,Location"))
    data = []
    for l in lines[i+2:]:
        if not l.strip():
            break
        data.append(l)
    df = pd.read_csv(io.StringIO("\n".join(data)), header=None, names=DECOMP_COLS)
    df["local_idx"] = df["Residue"].str.split().str[-1].astype(int)
    df["resname"] = df["Residue"].str.split().str[0]
    return df

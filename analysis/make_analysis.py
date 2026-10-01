#!/usr/bin/env python3
"""P02 final analysis. Run from the repo root:  python analysis/make_analysis.py

Inputs : data/mmpbsa/{gb,decomp,pb}, data/rmsd/core, data/docking
Outputs: results/*.csv and results/fig{1,2,3}_*.png

All free energies use a common receptor range (EGFR 701-989) so that receptor net charge and
construct length do not differ between genotypes (see README). WT is evaluated on its extended
window (15-20 ns); T790M and the triple mutant on 5-10 ns. The WT 5-10 ns window is kept as a
sensitivity row.
"""
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.append(str(Path(__file__).resolve().parent))
from parse_mmpbsa import delta_terms, parse_decomp_total

ROOT = Path(__file__).resolve().parents[1]
DATA, RES = ROOT / "data", ROOT / "results"
RES.mkdir(exist_ok=True)
SEEDS = (1, 2, 3)
LABEL = {"WT": "WT", "T790M": "T790M", "C797S_triple": "L858R/T790M/C797S*"}
COLOR = {"WT": "#6c788a", "T790M": "#c0392b", "C797S_triple": "#6d8672"}
PRIMARY_WIN = {"WT": "w15-20ns", "T790M": "w5-10ns", "C797S_triple": "w5-10ns"}
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 6.5,
                     "ytick.labelsize": 6.5, "legend.fontsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 300, "pdf.fonttype": 42})

# ------------------------------------------------------------------ 1. GB per run
rows = []
for g, w in [("WT", "w5-10ns"), ("WT", "w15-20ns"), ("T790M", "w5-10ns"), ("C797S_triple", "w5-10ns")]:
    for s in SEEDS:
        t = delta_terms(DATA / "mmpbsa" / "gb" / f"GB_{g}_seed{s}_{w}.dat")
        t.update(genotype=g, seed=s, window=w)
        rows.append(t)
runs = pd.DataFrame(rows)
runs["polar_sum"] = runs["EEL"] + runs["EGB"]       # EEL and EGB are only meaningful as a sum across constructs
runs.round(3).to_csv(RES / "P02_gb_runs_all.csv", index=False)
prim = pd.concat([runs[(runs.genotype == g) & (runs.window == w)] for g, w in PRIMARY_WIN.items()])
wt_init = runs[(runs.genotype == "WT") & (runs.window == "w5-10ns")]

TERMS = ["dG", "VDW", "polar_sum", "ESURF", "EEL", "EGB", "Ggas", "Gsolv"]
def summarize(df, name=None):
    out = []
    for (g, w), sub in df.groupby(["genotype", "window"], sort=False):
        r = {"genotype": name or LABEL[g], "window": w, "n_seeds": len(sub)}
        for t in TERMS:
            r[f"{t}_mean"], r[f"{t}_sd"] = sub[t].mean(), sub[t].std(ddof=1)
        out.append(r)
    return pd.DataFrame(out)
energy = pd.concat([summarize(prim), summarize(wt_init, "WT (initial 5-10 ns window, sensitivity)")])
energy.round(2).to_csv(RES / "P02_energy_table.csv", index=False)

def welch(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = a.mean() - b.mean()
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    se = np.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    tc = stats.t.ppf(0.975, df)
    return d, se, d - tc * se, d + tc * se, 2 * stats.t.sf(abs(d / se), df)

def pick(g, w=None):
    w = w or PRIMARY_WIN[g]
    return runs[(runs.genotype == g) & (runs.window == w)]
comp_rows = []
cmp_sets = [("T790M - WT", pick("T790M"), pick("WT")), ("T790M - triple", pick("T790M"), pick("C797S_triple")),
            ("WT - triple", pick("WT"), pick("C797S_triple")),
            ("T790M - WT (initial 5-10 ns)", pick("T790M"), pick("WT", "w5-10ns"))]
for name, A, B in cmp_sets:
    for t in ["dG", "VDW", "polar_sum", "ESURF"]:
        d, se, lo, hi, p = welch(A[t], B[t])
        comp_rows.append(dict(comparison=name, term=t, diff=d, se=se, ci95_lo=lo, ci95_hi=hi, welch_p=p))
ddg = pd.DataFrame(comp_rows)
ddg.round(3).to_csv(RES / "P02_ddG_components.csv", index=False)

# ------------------------------------------------------------------ 2. docking vs GBSA
dk = pd.DataFrame(json.load(open(DATA / "docking" / "p01_main_matrix_results.json")))
dk = dk[dk.ligand == "osimertinib"].groupby("genotype").top_score_kcal_mol.mean()
rank = pd.DataFrame({"genotype": [LABEL[g] for g in PRIMARY_WIN],
                     "docking_kcal_mol": [dk[g] for g in PRIMARY_WIN],
                     "dG_GBSA_mean": [pick(g).dG.mean() for g in PRIMARY_WIN],
                     "dG_GBSA_sd": [pick(g).dG.std(ddof=1) for g in PRIMARY_WIN]})
rank["docking_rank"] = rank.docking_kcal_mol.rank().astype(int)
rank["GBSA_rank"] = rank.dG_GBSA_mean.rank().astype(int)
rank.round(2).to_csv(RES / "P02_docking_vs_gbsa_ranking.csv", index=False)


# ------------------------------------------------------------------ 2b. sensitivity: full constructs (first pass)
full = []
for g in PRIMARY_WIN:
    for s in SEEDS:
        t = delta_terms(DATA / "mmpbsa" / "gb_full_construct" / f"GBfull_{g}_seed{s}_w5-10ns.dat")
        t.update(genotype=g, seed=s)
        full.append(t)
full = pd.DataFrame(full)
fs = full.groupby("genotype", sort=False).agg(dG_mean=("dG", "mean"), dG_sd=("dG", lambda x: x.std(ddof=1)),
                                              EEL_mean=("EEL", "mean"), EGB_mean=("EGB", "mean"))
cmn = runs[runs.window == "w5-10ns"].groupby("genotype", sort=False).agg(dG_common_mean=("dG", "mean"))
sens = fs.join(cmn)
sens["shift_common_minus_full"] = sens.dG_common_mean - sens.dG_mean
sens.round(2).to_csv(RES / "P02_full_construct_sensitivity.csv")
print("full-construct vs common-range (all 5-10 ns):\n", sens.round(2).to_string())

# ------------------------------------------------------------------ 3. per-residue decomposition
def seed_matrix(g):
    cols, lig = {}, []
    for s in SEEDS:
        df = parse_decomp_total(DATA / "mmpbsa" / "decomp" / f"DECOMPfull_{g}_seed{s}_{PRIMARY_WIN[g]}.dat")
        lig.append(df[df.resname == "MOL"].TOTAL_avg.iloc[0])
        p = df[df.resname != "MOL"].copy()
        p["resnum"] = p.local_idx + 700                       # common range: local 1 == EGFR 701 ... 289 == 989
        cols[s] = p.set_index("resnum").TOTAL_avg
    return pd.DataFrame(cols), np.array(lig)
M = {g: seed_matrix(g) for g in PRIMARY_WIN}
assert all(M[g][0].index.equals(M["WT"][0].index) for g in M) and M["WT"][0].index[[0, -1]].tolist() == [701, 989]
res = pd.DataFrame({"resnum": M["WT"][0].index})
for g in M:
    res[f"{g}_mean"], res[f"{g}_sd"] = M[g][0].mean(axis=1).values, M[g][0].std(axis=1, ddof=1).values
res["d_T790M_WT"] = res.T790M_mean - res.WT_mean
res["se_T790M_WT"] = np.sqrt(res.T790M_sd ** 2 / 3 + res.WT_sd ** 2 / 3)
res["z_T790M_WT"] = res.d_T790M_WT / res.se_T790M_WT
res["d_T790M_triple"] = res.T790M_mean - res.C797S_triple_mean
res["se_T790M_triple"] = np.sqrt(res.T790M_sd ** 2 / 3 + res.C797S_triple_sd ** 2 / 3)
res.round(3).to_csv(RES / "P02_residue_delta.csv", index=False)
lig_self = pd.DataFrame({LABEL[g]: [M[g][1].mean(), M[g][1].std(ddof=1)] for g in M}, index=["mean", "sd"])
lig_self.round(2).to_csv(RES / "P02_ligand_self_term.csv")

# ------------------------------------------------------------------ 4. RMSD / contacts
def load_series(g, s):
    f = DATA / "rmsd" / "core" / (f"rmsdcoreext_WT_seed{s}.dat" if g == "WT" else f"rmsdcore_{g}_seed{s}.dat")
    d = np.loadtxt(f)
    c = np.loadtxt(DATA / "rmsd" / "core" / f"contacts_{g}_seed{s}.dat")
    n = min(len(d), len(c))
    t = (d[:n, 0] - 1) * 0.02
    return t, d[:n, 1], d[:n, 2], c[:n, 1] / c[0, 1]
stat_rows = []
WINS = {"WT": [("5-10 ns", 250, 500), ("15-20 ns", 750, 1000)], "T790M": [("5-10 ns", 250, 500)],
        "C797S_triple": [("5-10 ns", 250, 500)]}
SER = {(g, s): load_series(g, s) for g in WINS for s in SEEDS}
for g in WINS:
    for s in SEEDS:
        t, core, lig, nat = SER[(g, s)]
        for wn, a, b in WINS[g]:
            r = dict(genotype=LABEL[g], seed=s, window=wn)
            for nm, y in (("core_CA", core), ("ligand", lig), ("native_contact_frac", nat)):
                r[f"{nm}_mean"] = y[a:b].mean()
                r[f"{nm}_slope_per_ns"] = np.polyfit(t[a:b], y[a:b], 1)[0]
            stat_rows.append(r)
rst = pd.DataFrame(stat_rows)
ref = rst[(rst.genotype != "WT") & (rst.window == "5-10 ns")]
thr_core, thr_lig = ref.core_CA_slope_per_ns.abs().max(), ref.ligand_slope_per_ns.abs().max()
rst["plateau_core"] = rst.core_CA_slope_per_ns.abs() <= thr_core
rst["plateau_ligand"] = rst.ligand_slope_per_ns.abs() <= thr_lig
rst.round(3).to_csv(RES / "P02_rmsd_core_stats.csv", index=False)

# ------------------------------------------------------------------ 5. figures
SEEDS_MK = {1: "o", 2: "s", 3: "^"}
# Fig 1
fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.3), gridspec_kw={"width_ratios": [1.1, 1.2]})
ax = axes[0]
cats = [("WT", "w5-10ns", "WT\n5-10 ns"), ("WT", "w15-20ns", "WT\n15-20 ns"),
        ("T790M", "w5-10ns", "T790M\n5-10 ns"), ("C797S_triple", "w5-10ns", "L858R/T790M/\nC797S*\n5-10 ns")]
for i, (g, w, lab) in enumerate(cats):
    sub = pick(g, w)
    alpha = 0.45 if (g == "WT" and w == "w5-10ns") else 1.0
    ax.bar(i, sub.dG.mean(), width=0.6, color=COLOR[g], alpha=alpha * 0.85, zorder=2)
    ax.errorbar(i, sub.dG.mean(), yerr=sub.dG.std(ddof=1), fmt="none", ecolor="k", capsize=3, lw=0.9, zorder=3)
    for s, v in zip(sub.seed, sub.dG):
        ax.scatter(i + (s - 2) * 0.12, v, color="k", s=14, marker=SEEDS_MK[s], zorder=4)
ax.set_xticks(range(len(cats))); ax.set_xticklabels([c[2] for c in cats])
ax.set_ylabel(r"$\Delta G_{GBSA}$ (kcal/mol)")
ax.axhline(0, color="0.5", lw=0.7)
def _g(cmp_, term="dG"):
    return ddg[(ddg.comparison == cmp_) & (ddg.term == term)].iloc[0]
d_tw, d_tt = _g("T790M - WT"), _g("T790M - triple")
ax.set_title(f"T790M ranks first: {d_tt['diff']:.0f} kcal/mol vs triple mutant\n(p = {d_tt['welch_p']:.2f}); "
             f"{d_tw['diff']:.0f} vs WT is not resolved (p = {d_tw['welch_p']:.2f})", loc="left")
ax.margins(y=0.1)
ax.text(0.01, 0.02, "n = 3 seeds; markers = seed 1/2/3; bars = mean ± SD", transform=ax.transAxes, fontsize=6, color="0.35")
ax = axes[1]
terms = [("VDW", "van der Waals"), ("polar_sum", "Polar (EEL+EGB)"), ("ESURF", "Nonpolar solv.")]
w = 0.26
for k, g in enumerate(PRIMARY_WIN):
    sub = pick(g)
    for j, (t, _) in enumerate(terms):
        x = j + (k - 1) * w
        ax.bar(x, sub[t].mean(), width=w * 0.9, color=COLOR[g], zorder=2, label=LABEL[g] if j == 0 else None)
        ax.errorbar(x, sub[t].mean(), yerr=sub[t].std(ddof=1), fmt="none", ecolor="k", capsize=2, lw=0.8, zorder=3)
        for s, v in zip(sub.seed, sub[t]):
            ax.scatter(x, v, color="k", s=7, marker=SEEDS_MK[s], zorder=4)
ax.axhline(0, color="0.5", lw=0.7)
ax.set_xticks(range(len(terms))); ax.set_xticklabels([n for _, n in terms])
ax.set_ylabel("kcal/mol")
dv = _g("T790M - triple", "VDW")["diff"]
ax.set_title(f"{dv / d_tt['diff']:.0%} of the T790M-vs-triple gap\nis van der Waals", loc="left")
ax.legend(frameon=False, loc="lower right")
ax.margins(y=0.12)
for a, l in zip(axes, "ab"):
    a.text(-0.12, 1.05, l, transform=a.transAxes, fontsize=10, fontweight="bold")
fig.tight_layout()
fig.savefig(RES / "fig1_energy_components.png")
plt.close(fig)

# Fig 2
fig, axes = plt.subplots(2, 1, figsize=(7.4, 5.2), sharex=True, gridspec_kw={"height_ratios": [1, 1.25]})
ax = axes[0]
for g in ("WT", "T790M"):
    ax.plot(res.resnum, res[f"{g}_mean"], color=COLOR[g], lw=1.1, label=f"{LABEL[g]} (n=3 seeds)")
    ax.fill_between(res.resnum, res[f"{g}_mean"] - res[f"{g}_sd"], res[f"{g}_mean"] + res[f"{g}_sd"], color=COLOR[g], alpha=0.18, lw=0)
ax.axvline(790, color="0.5", lw=0.8, ls=(0, (3, 2)))
ax.text(790, ax.get_ylim()[1], " 790", fontsize=6.5, color="0.35", va="top")
ax.axhline(0, color="0.6", lw=0.5)
ax.set_ylabel("Per-residue contribution\n(kcal/mol)")
ax.set_title("Per-residue contribution to osimertinib binding (shaded = ±SD across seeds)", loc="left")
ax.legend(frameon=False, loc="upper right")
ax = axes[1]
sig = res.d_T790M_WT.abs() > 2 * res.se_T790M_WT
vmax = res.d_T790M_WT.abs().max()
cols = [plt.get_cmap("RdBu_r")(plt.Normalize(-vmax, vmax)(v)) if s else (0.82, 0.82, 0.82, 1) for v, s in zip(res.d_T790M_WT, sig)]
ax.bar(res.resnum, res.d_T790M_WT, width=1.0, color=cols, zorder=2)
ax.axhline(0, color="0.6", lw=0.5); ax.axvline(790, color="0.5", lw=0.8, ls=(0, (3, 2)))
order = res[sig].reindex(res[sig].d_T790M_WT.abs().sort_values(ascending=False).index)
picked = []
for _, r in order.iterrows():
    if all(abs(r.resnum - p) > 2 for p in picked):
        picked.append(r.resnum)
    if len(picked) == 5:
        break
for p in picked:
    v = float(res.loc[res.resnum == p, "d_T790M_WT"].iloc[0])
    ax.annotate(f"{int(p)}", xy=(p, v), xytext=(p, v - 0.55 if v < 0 else v + 0.55), fontsize=6.5, ha="center",
                va="top" if v < 0 else "bottom", arrowprops=dict(arrowstyle="-", lw=0.6, color="0.5"))
ax.set_ylim(res.d_T790M_WT.min() - 1.2, res.d_T790M_WT.max() + 1.0)
ax.set_ylabel("T790M − WT\n(kcal/mol; lower = T790M binds more)")
ax.set_xlabel("EGFR residue number (common range 701-989)")
ax.set_title("Where T790M gains binding energy (grey = within 2 SE of zero)", loc="left")
for a, l in zip(axes, "ab"):
    a.text(-0.07, 1.05, l, transform=a.transAxes, fontsize=10, fontweight="bold")
fig.tight_layout()
fig.savefig(RES / "fig2_residue_contribution.png")
plt.close(fig)

# Fig 3
fig, axes = plt.subplots(3, 3, figsize=(7.6, 6.2), sharex="col", sharey="row")
rowdefs = [("Core CA RMSD (Å)", 1), ("Ligand RMSD (Å)", 2), ("Native contacts kept", 3)]
shade = {"WT": [(5, 10, 0.07), (15, 20, 0.16)], "T790M": [(5, 10, 0.16)], "C797S_triple": [(5, 10, 0.16)]}
for c, g in enumerate(["WT", "T790M", "C797S_triple"]):
    for r_, (ylab, k) in enumerate(rowdefs):
        ax = axes[r_, c]
        for a, b, al in shade[g]:
            ax.axvspan(a, b, color="0.4", alpha=al, lw=0)
        for s in SEEDS:
            t, core, lig, nat = SER[(g, s)]
            y = (core, lig, nat)[k - 1]
            ax.plot(t, y, lw=0.7, color=COLOR[g], alpha=[0.45, 0.75, 1.0][s - 1], label=f"seed {s}")
        if c == 0:
            ax.set_ylabel(ylab)
        if r_ == 0:
            ax.set_title(LABEL[g], loc="left")
        if r_ == 2:
            ax.set_xlabel("Time (ns)")
axes[0, 0].legend(frameon=False, loc="lower right")
fig.tight_layout(rect=(0, 0.025, 1, 1))
fig.text(0.5, 0.004, "Grey bands: MM-GBSA windows. Rows share the y-axis; WT runs 20 ns (10 ns production + 10 ns extension).",
         ha="center", fontsize=6, color="0.35")
fig.savefig(RES / "fig3_rmsd_core_contacts.png")
plt.close(fig)

# ------------------------------------------------------------------ 6. console summary
pd.set_option("display.width", 200)
print(energy[["genotype", "window", "dG_mean", "dG_sd", "VDW_mean", "polar_sum_mean", "ESURF_mean"]].round(2).to_string(index=False))
print(ddg.round(2).to_string(index=False))
print(rank.round(2).to_string(index=False))
print("ligand self term:\n", lig_self.round(2))
print("plateau thresholds (max |slope| in T790M/triple 5-10 ns): core %.3f  ligand %.3f A/ns" % (thr_core, thr_lig))
print(rst[["genotype", "seed", "window", "core_CA_mean", "core_CA_slope_per_ns", "ligand_mean", "ligand_slope_per_ns",
           "native_contact_frac_mean", "plateau_core", "plateau_ligand"]].round(3).to_string(index=False))
print("top movers (|dev|>2SE):\n", res[sig].reindex(res[sig].d_T790M_WT.abs().sort_values(ascending=False).index)
      [["resnum", "WT_mean", "T790M_mean", "d_T790M_WT", "se_T790M_WT"]].head(12).round(2).to_string(index=False))
print("residues 790/797/858/865-867:\n", res[res.resnum.isin([790, 797, 858, 865, 866, 867])]
      [["resnum", "WT_mean", "T790M_mean", "C797S_triple_mean", "d_T790M_WT", "se_T790M_WT"]].round(2).to_string(index=False))
print("decomp sums (residues + ligand self) per genotype:",
      {g: round(float(M[g][0].sum().mean() + M[g][1].mean()), 1) for g in M})

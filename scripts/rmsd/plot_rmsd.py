"""P02 RMSD platosu: 9 trajektori, protein CA + ligand ağır atom. MM-GBSA penceresi kare 251-500."""
import csv, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).parent
GENOTYPES = [("WT", "WT"), ("T790M", "T790M"), ("C797S_triple", "L858R/T790M/C797S")]
DT_NS = 0.02
W0 = 250  # 0-indexed kare 250 = 251. kare = 5.00 ns

rows = []
fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex=True, sharey="row")
for col, (key, label) in enumerate(GENOTYPES):
    for seed in (1, 2, 3):
        data = np.loadtxt(HERE / f"rmsd_{key}_seed{seed}.dat")
        t = (data[:, 0] - 1) * DT_NS
        for row, (name, series) in enumerate((("protein CA", data[:, 1]), ("ligand", data[:, 2]))):
            ax = axes[row, col]
            ax.plot(t, series, lw=0.8, label=f"seed {seed}")
            win_t, win = t[W0:], series[W0:]
            slope = np.polyfit(win_t, win, 1)[0]
            half = len(win) // 2
            rows.append({"genotype": label, "seed": seed, "selection": name,
                         "window_mean_A": round(float(win.mean()), 3),
                         "window_std_A": round(float(win.std()), 3),
                         "slope_A_per_ns": round(float(slope), 3),
                         "mean_5_7.5ns_A": round(float(win[:half].mean()), 3),
                         "mean_7.5_10ns_A": round(float(win[half:].mean()), 3),
                         "max_A": round(float(series.max()), 3)})
    for row in (0, 1):
        axes[row, col].axvspan(5, 10, color="0.9", zorder=0)
    axes[0, col].set_title(label)
axes[0, 0].set_ylabel("Protein CA RMSD (Å)")
axes[1, 0].set_ylabel("Ligand ağır atom RMSD (Å)\n(protein hizasında)")
for ax in axes[1]:
    ax.set_xlabel("Zaman (ns)")
axes[0, 0].legend(fontsize=8, loc="lower right")
fig.suptitle("P02 — RMSD; gri alan = MM-GBSA penceresi (son 5 ns, kare 251–500)")
fig.tight_layout()
fig.savefig(HERE.parent / "fig3_rmsd_plateau.png", dpi=150)

with open(HERE.parent / "P02_rmsd_plateau.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
for r in rows:
    print(r)

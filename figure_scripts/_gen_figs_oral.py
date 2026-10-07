# -*- coding: utf-8 -*-
"""Figures for the rewritten Indian oral-cancer manuscript."""
import os, sys, csv, json
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding="utf-8")

FIG = r"./figures_cs"
EDGES = r"./checkpoints/edge_list.csv"
IND = r"{DATA_ROOT}/cancer_application/data/indian_oral"

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "xtick.direction": "out", "ytick.direction": "out",
    "figure.dpi": 150, "savefig.bbox": "tight",
})
ACC = "#B03A6E"      # emphasis (inflammatory)
BASE = "#9E9BC8"     # comparator cohort / baseline
NEU = "#8A8A8A"

DE = json.load(open(os.path.join(IND, "indian_genomewide_de.json")))
P = DE["p"]; D = DE["delta"]
FN = json.load(open(os.path.join(IND, "oral_final_numbers.json")))
REPL = json.load(open(os.path.join(IND, "indian_edge_replication.json")))
ROB = json.load(open(os.path.join(IND, "indian_robust.json")))

# network degrees
deg_out, deg_in = Counter(), Counter()
edges_disc = []
with open(EDGES, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        s = r["source_gene"].upper(); t = r["target_gene"].upper()
        src = (r.get("source") or "").strip()
        if src != "TRRUST" and not src.startswith("pathway"):
            deg_out[s] += 1; deg_in[t] += 1
            edges_disc.append((s, t))
print("discovered edges", len(edges_disc))

# ------------------------------------------------------------------ Fig 1
fig = plt.figure(figsize=(7.2, 7.4))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.55], hspace=0.75, wspace=0.55)

ax = fig.add_subplot(gs[0, 0])
top = deg_out.most_common(12)
names = [g for g, _ in top][::-1]; vals = [v for _, v in top][::-1]
cols = [ACC if g in ("NFKB1", "RELA", "STAT3", "STAT1", "JUN") else NEU for g in names]
ax.barh(names, vals, color=cols, height=0.72)
ax.set_xlabel("connections emanating", fontsize=8)
ax.tick_params(labelsize=7.5)
ax.set_title("a", loc="left", fontweight="bold", fontsize=10)

ax = fig.add_subplot(gs[0, 1])
top = deg_in.most_common(12)
names = [g for g, _ in top][::-1]; vals = [v for _, v in top][::-1]
cols = [ACC if g in ("CXCL8", "PTGS2", "TNF", "CCND1", "VEGFA") else NEU for g in names]
ax.barh(names, vals, color=cols, height=0.72)
ax.set_xlabel("connections terminating", fontsize=8)
ax.tick_params(labelsize=7.5)
ax.set_title("b", loc="left", fontweight="bold", fontsize=10)

# panel c: inflammatory core subgraph
ax = fig.add_subplot(gs[1:, :])
CORE_REG = ["NFKB1", "RELA", "STAT1", "STAT3", "JUN", "FOSL1", "ETS1", "SP1"]
CORE_TGT = ["CXCL8", "CXCL10", "MMP9", "PTGS2", "TNF", "CCND1", "CDKN1A", "BCL2", "VEGFA", "IL6"]
reg = [g for g in CORE_REG]
tgt = [g for g in CORE_TGT]
alledges = {}
with open(EDGES, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        s = r["source_gene"].upper(); t = r["target_gene"].upper()
        if s in reg and t in tgt:
            w = abs(float(r["weight"]))
            alledges[(s, t)] = max(alledges.get((s, t), 0.0), w)
sub = list(alledges.keys())
print("core subgraph edges:", len(sub), sorted(sub)[:14])
yp = {g: i for i, g in enumerate(reg)}
yt = {g: i for i, g in enumerate(tgt)}
n_r, n_t = len(reg), len(tgt)
for (s, t) in sub:
    y0 = (n_r - 1 - yp[s]) / max(n_r - 1, 1)
    y1 = (n_t - 1 - yt[t]) / max(n_t - 1, 1)
    lw = 0.45 + 1.1 * alledges[(s, t)]
    ax.annotate("", xy=(1.0, y1), xytext=(0.0, y0),
                arrowprops=dict(arrowstyle="-|>", color=ACC, lw=lw, alpha=0.5,
                                shrinkA=6, shrinkB=6))
for g, i in yp.items():
    y = (n_r - 1 - i) / max(n_r - 1, 1)
    ax.scatter([0], [y], s=70, color=ACC, zorder=3)
    ax.text(-0.05, y, g, ha="right", va="center", fontsize=8.5, fontweight="bold")
for g, i in yt.items():
    y = (n_t - 1 - i) / max(n_t - 1, 1)
    ax.scatter([1], [y], s=70, color=BASE, zorder=3)
    ax.text(1.05, y, g, ha="left", va="center", fontsize=8.5, fontweight="bold")
ax.set_xlim(-0.42, 1.42); ax.set_ylim(-0.12, 1.12)
ax.axis("off")
ax.text(0.0, 1.09, "regulators", ha="center", fontsize=8.5, color=ACC)
ax.text(1.0, 1.09, "effectors", ha="center", fontsize=8.5, color=BASE)
ax.set_title("c", loc="left", fontweight="bold", fontsize=10)
fig.savefig(os.path.join(FIG, "fig1_core.pdf"))
plt.close(fig)
print("fig1 done")

# ------------------------------------------------------------------ Fig 2
order = ["RELA", "STAT1", "STAT3", "JUN", "FOSL1", "ETS1", "E2F1",
         "CXCL8", "CXCL10", "MMP9", "PTGS2", "IL6", "CCL5", "IL1B",
         "CCND1", "CDKN1A", "VEGFA", "SP1", "NFKB1"]
order = [g for g in order if g in P["GSE85195"] and g in P["GSE23558"]]
d1 = [D["GSE85195"][g] for g in order]
d2 = [D["GSE23558"][g] for g in order]

fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.5, 1], hspace=0.5, wspace=0.28)

ax = fig.add_subplot(gs[0, :])
y = np.arange(len(order))
ax.barh(y - 0.2, d1, height=0.38, color=ACC, label="GSE85195 (tumour vs leukoplakia)")
ax.barh(y + 0.2, d2, height=0.38, color=BASE, label="GSE23558 (tumour vs normal)")
ax.set_yticks(y); ax.set_yticklabels(order, fontsize=8)
ax.invert_yaxis()
ax.axvline(0, color="k", lw=0.7)
ax.set_xlabel("difference of group means (log$_2$)", fontsize=8)
ax.legend(fontsize=7.5, frameon=False, loc="lower right")
ax.set_title("a", loc="left", fontweight="bold", fontsize=10)

ax = fig.add_subplot(gs[1, 0])
fish = FN["fisher"]
f9 = [r for r in fish if r[0] in ("STAT1", "FOSL1", "CXCL8", "MMP9", "E2F1",
                                  "CXCL10", "JUN", "RELA", "ETS1")]
f9 = sorted(f9, key=lambda r: r[5], reverse=True)
lbl = [r[0] for r in f9]
val = [-np.log10(r[5]) for r in f9]
ax.bar(lbl, val, color=ACC, width=0.7)
ax.set_ylabel("$-\\log_{10}$ combined $p$", fontsize=8)
ax.tick_params(labelsize=7.5, axis="x", rotation=60)
ax.set_title("b", loc="left", fontweight="bold", fontsize=10)

ax = fig.add_subplot(gs[1, 1])
ax.axis("off")
nine = ["STAT1", "FOSL1", "CXCL8", "MMP9", "E2F1", "CXCL10", "JUN", "RELA", "ETS1"]
txt = ("up-regulated in both cohorts\n"
       "(false-discovery rate 5%)\n\n"
       + ", ".join(nine[:5]) + ",\n" + ", ".join(nine[5:]))
ax.text(0.0, 1.0, txt, va="top", ha="left", fontsize=8.5, linespacing=1.7)
ax.set_title("c", loc="left", fontweight="bold", fontsize=10)
fig.savefig(os.path.join(FIG, "fig2_cohorts.pdf"))
plt.close(fig)
print("fig2 done")

# ------------------------------------------------------------------ Fig 3
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw={"wspace": 0.34})

ax = axes[0]
cohorts = ["GSE85195", "GSE23558"]
x = np.arange(len(cohorts)); w = 0.34
obs = [REPL[c]["mean_abs_rho_edges"] for c in cohorts]
nul = [REPL[c]["mean_abs_rho_null"] for c in cohorts]
ax.bar(x - w / 2, obs, w, color=ACC, label="network edges")
ax.bar(x + w / 2, nul, w, color=BASE, label="degree-matched null")
ax.set_xticks(x); ax.set_xticklabels(cohorts, fontsize=8)
ax.set_ylabel("mean |Spearman correlation|", fontsize=8)
for i, c in enumerate(cohorts):
    ax.text(i, max(obs[i], nul[i]) + 0.012, f"fold {REPL[c]['fold']:.2f}\n$p$={REPL[c]['p']:.2f}",
            ha="center", fontsize=7.5)
ax.set_ylim(0, max(max(obs), max(nul)) * 1.42)
ax.legend(fontsize=7.5, frameon=False)
ax.set_title("a", loc="left", fontweight="bold", fontsize=10)

ax = axes[1]
bars = [0.0, 1.0]
ax.bar([0], [0.0], 0.55, color=ACC)
ax.bar([1], [1.0], 0.55, color=BASE)
ax.set_xticks([0, 1])
ax.set_xticklabels(["estimator-\nproposed", "TRRUST-\nderived"], fontsize=8)
ax.set_ylabel("fraction with curated orientation", fontsize=8)
ax.set_ylim(0, 1.18)
ax.text(0, 0.03, "0 / 104", ha="center", fontsize=8.5, fontweight="bold")
ax.text(1, 1.02, "8,427 / 8,427", ha="center", fontsize=8.5)
ax.set_title("b", loc="left", fontweight="bold", fontsize=10)
fig.savefig(os.path.join(FIG, "fig3_transfer.pdf"))
plt.close(fig)
print("fig3 done")

# ------------------------------------------------------------------ Fig S2
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"wspace": 0.3})
for ax, c in zip(axes, cohorts):
    rows = ROB[c]["expr_matched"]
    xs = [r["decile"] for r in rows]
    a = [r["up_net"] for r in rows]
    b = [r["up_oth"] for r in rows]
    ax.plot(xs, b, "o-", color=BASE, ms=3.5, lw=1.2, label="other genes")
    ax.plot(xs, a, "o-", color=ACC, ms=3.5, lw=1.2, label="network genes")
    ax.set_xticks(xs)
    ax.set_xlabel("decile of mean expression", fontsize=8)
    ax.set_ylabel("fraction up-regulated", fontsize=8)
    ax.tick_params(labelsize=7.5)
    ax.set_title(c, fontsize=8.5)
axes[0].legend(fontsize=7.5, frameon=False)
fig.savefig(os.path.join(FIG, "figS2_robustness.pdf"))
plt.close(fig)
print("figS2 done")

for f in ["fig1_core.pdf", "fig2_cohorts.pdf", "fig3_transfer.pdf", "figS2_robustness.pdf"]:
    p = os.path.join(FIG, f)
    print(f, os.path.getsize(p), "bytes")

try:
    import fitz
    for f in ["fig1_core.pdf", "fig2_cohorts.pdf", "fig3_transfer.pdf", "figS2_robustness.pdf"]:
        d = fitz.open(os.path.join(FIG, f))
        d[0].get_pixmap(dpi=140).save(os.path.join(FIG, f.replace(".pdf", ".png")))
    print("png previews written")
except Exception as e:
    print("no fitz:", e)

"""
Figure 6 generator - extreme scale (up to 100M) + finance + Sachs DAG.
"""
import sys, os, json

CHK = r"C:\Users\高帅东\Desktop\ICLR2027_Submission_Clean\checkpoints"
OUT = r"C:\Users\高帅东\Desktop\ICLR2027_Submission_Clean\figures"
os.makedirs(OUT, exist_ok=True)

ext = json.load(open(os.path.join(CHK, "extreme_scale.json"), encoding="utf-8"))
fin = json.load(open(os.path.join(CHK, "finance_result.json"), encoding="utf-8"))
sd  = json.load(open(os.path.join(CHK, "sachs_dag.json"), encoding="utf-8"))
sf  = json.load(open(os.path.join(CHK, "sachs_fullrank.json"), encoding="utf-8"))

# --- Panel (a): ALL 6 scale test entries ---
all_st = ext["scale_tests"]
# Items with recovery_pct
rec_items = [t for t in all_st if "recovery_pct" in t]
# Items without (sparse mode, huge d)
sparse_items = [t for t in all_st if "recovery_pct" not in t]

d_rec  = [t["d"] for t in rec_items]
pct    = [t["recovery_pct"] for t in rec_items]

# All items: unify time (some have time_s, some have train_time_s)
def get_time(t):
    return t.get("time_s", t.get("train_time_s", 0))

all_d = [t["d"] for t in all_st]
all_t = [get_time(t) for t in all_st]
modes = ["full" if "recovery_pct" in t else "sparse" for t in all_st]

print("Recovery points:", list(zip(d_rec, pct)))
print("All d:", all_d)
print("All t:", all_t)
print("Modes:", modes)

# --- Plot ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

d_rec  = np.array(d_rec)
pct    = np.array(pct)
all_d  = np.array(all_d)
all_t  = np.array(all_t)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 5.5))

# ====== (a) Extreme Scale: dual-axis ======
ax1.plot(d_rec, pct, "o-", color="#27AE60", lw=3, ms=11,
         mfc="#27AE60", mec="white", mew=1.5, zorder=3, label="Recovery %")
# Only label first and last recovery points to avoid log-scale crowding
ax1.annotate("99.9%", (d_rec[0], pct[0] + 3), ha="center", fontsize=8,
             color="#1E8449", fontweight="bold")
ax1.annotate("89.1% at d=19K\n(degradation onset)", (150000, 86),
             ha="center", fontsize=7.5, color="#C0392B", fontweight="bold")
# Compact summary for middle points — shifted far right
ax1.annotate(">99% through d=10K", (150000, 95), fontsize=7.5, color="#1E8449",
             ha="center", fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                       edgecolor="#27AE60", alpha=0.85))

# Right axis: time
ax1t = ax1.twinx()
full_idx = [i for i, m in enumerate(modes) if m == "full"]
sparse_idx = [i for i, m in enumerate(modes) if m == "sparse"]

ax1t.plot(all_d[full_idx], all_t[full_idx], "D-", color="#2980B9", lw=2, ms=9,
          mfc="#2980B9", mec="white", mew=1, zorder=2, label="Time (full)")
ax1t.plot(all_d[sparse_idx], all_t[sparse_idx], "s--", color="#E74C3C", lw=2, ms=11,
          mfc="#E74C3C", mec="white", mew=1.5, zorder=2, label="Time (sparse)")

# Sparse annotations: placed ON the right-axis (time) coordinates, near red squares
# Use small xytext multiplier to keep text boxes anchored left of the squares
sparse_info = [(20, 3000, 15.6), (100, 3000, 738.1)]
for i, (d_m, edges, t_val) in zip(sparse_idx, sparse_info):
    ax1t.annotate("d=%dM sparse\n%d edges, %.0fs" % (d_m, edges, t_val),
                xy=(all_d[i], all_t[i]),
                xytext=(all_d[i] * 0.015, all_t[i] * 0.08),
                fontsize=7.5, color="#E74C3C", fontweight="bold", ha="left",
                arrowprops=dict(arrowstyle="->", color="#E74C3C", lw=1.5),
                bbox=dict(boxstyle="round,pad=0.35", facecolor="#FDEDEC",
                          edgecolor="#E74C3C", alpha=0.92))

ax1.set_xscale("log")
ax1t.set_yscale("log")
ax1.set_xlabel("Dimension (d, log scale)", fontsize=11, fontweight="bold")
ax1.set_ylabel("Edge Recovery (%)", fontsize=11, fontweight="bold", color="#27AE60")
ax1t.set_ylabel("Wall Time (s, log scale)", fontsize=11, fontweight="bold", color="#E74C3C")
ax1.set_title("Extreme Scale (d=2K to 100M)", fontsize=11, fontweight="bold", pad=8)
ax1.set_ylim(0, 105)
ax1.grid(True, alpha=0.2)
ax1.text(-0.08, 1.02, "(a)", transform=ax1.transAxes, fontsize=14,
         fontweight="bold", va="bottom")
# Legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax1t.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=7.5, loc="lower right", framealpha=0.9)

# ====== (b) Financial Panel ======
f1_vals = [fin["f1"], fin["dagma_f1"]]
bars_f = ax2.bar(["LowRankGNN", "DAGMA"], f1_vals,
                 color=["#27AE60", "#E74C3C"], edgecolor="black", lw=0.5, width=0.5)
for bar, val in zip(bars_f, f1_vals):
    ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
             "F1=%.3f" % val, ha="center", fontsize=11, fontweight="bold")
ax2.set_ylabel("F1 Score", fontsize=11, fontweight="bold")
ax2.set_title("Cross-Domain: Financial Panel (d=%d, n=%d)" % (fin["d"], fin["n"]),
              fontsize=11, fontweight="bold", pad=8)
ax2.grid(True, alpha=0.2, axis="y")
ax2.set_ylim(0, max(f1_vals) * 1.35)
ax2.text(-0.08, 1.02, "(b)", transform=ax2.transAxes, fontsize=14,
         fontweight="bold", va="bottom")
gap = f1_vals[0] / f1_vals[1]
ax2.annotate("%.1fx gap" % gap, xy=(0, f1_vals[0]),
             xytext=(0.5, f1_vals[0] * 1.2), fontsize=9, color="#1E8449",
             fontweight="bold", ha="center",
             arrowprops=dict(arrowstyle="->", color="#1E8449", lw=1.2))

# ====== (c) Sachs DAG Constraint ======
f1s   = [sd["f1"], sf["f1"]]
times = [sd["time_s"], sf["time_s"]]
bars_c = ax3.bar(["Low-Rank+DAG", "Full-Rank DAG"], f1s,
                 color=["#27AE60", "#95A5A6"], edgecolor="black", lw=0.5, width=0.5)
ax3_twin = ax3.twinx()
ax3_twin.plot([0, 1], times, "D-", color="#E74C3C", lw=2.5, ms=10,
              mfc="#E74C3C", mec="white", mew=1.5)
ax3_twin.set_ylim(0, max(times) * 1.35)
for bar, val in zip(bars_c, f1s):
    ax3.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
             "F1=%.3f" % val, ha="center", fontsize=10, fontweight="bold")
for i, t in enumerate(times):
    ax3_twin.text(i + 0.15, t + 0.3, "%.1fs" % t, fontsize=9,
                  color="#E74C3C", fontweight="bold")
ax3.set_ylabel("F1 Score", fontsize=11, fontweight="bold")
ax3_twin.set_ylabel("Time (s)", fontsize=11, fontweight="bold", color="#E74C3C")
ax3.set_title("Sachs: DAG Constraint Benefit (d=11, n=7466)",
              fontsize=11, fontweight="bold", pad=8)
ax3.grid(True, alpha=0.2, axis="y")
ax3.set_ylim(0, max(f1s) * 1.3)
ax3.text(-0.08, 1.02, "(c)", transform=ax3.transAxes, fontsize=14,
         fontweight="bold", va="bottom")

plt.tight_layout()
png_path = os.path.join(OUT, "fig6_failure.png")
pdf_path = os.path.join(OUT, "fig6_failure.pdf")
fig.savefig(png_path, dpi=250, bbox_inches="tight", facecolor="white")
fig.savefig(pdf_path, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("FIG6 DONE: PNG=%d KB  PDF=%d KB" % (
    os.path.getsize(png_path) // 1024, os.path.getsize(pdf_path) // 1024))

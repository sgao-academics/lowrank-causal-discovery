"""Figure 5: Ablation — rank violation, multi-scale benefit, cross-domain validation."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json, os

CHK = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'checkpoints')
OUT = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(OUT, exist_ok=True)

with open(os.path.join(CHK, 'rank_violation.json')) as f:
    rv_data = json.load(f)
with open(os.path.join(CHK, 'multiscale_ablation.json')) as f:
    ms_data = json.load(f)
with open(os.path.join(CHK, 'sachs_result.json')) as f:
    sachs = json.load(f)

density  = np.array([d['density'] for d in rv_data])
eff_rank = np.array([d['eff_rank'] for d in rv_data])
rv_f1    = np.array([d['f1'] for d in rv_data])
rv_prec  = np.array([d['prec'] for d in rv_data])
rv_rec   = np.array([d['rec'] for d in rv_data])

fig = plt.figure(figsize=(16, 5.5))

# ===== (a): Rank Violation =====
ax1 = fig.add_subplot(1, 3, 1)
ax1.plot(eff_rank, rv_f1, 'o-', color='#27AE60', lw=2.5, ms=10, mfc='#27AE60', mec='white', mew=1.5, label='F1')
ax1.plot(eff_rank, rv_prec, 's--', color='#2980B9', lw=2, ms=9, mfc='#2980B9', mec='white', mew=1, label='Precision')
ax1.plot(eff_rank, rv_rec, '^--', color='#E67E22', lw=2, ms=9, mfc='#E67E22', mec='white', mew=1, label='Recall')

ax1.axvline(x=52, color='#7F8C8D', linestyle=':', lw=1, alpha=0.5)
ax1.text(54, 0.88, 'Full rank', fontsize=8, color='#7F8C8D', fontweight='bold')

ax1.annotate('F1: 0.11-0.28, stable\nacross 2x rank change',
            xy=(70, 0.26), fontsize=8.5, color='#1E8449', fontweight='bold', ha='center',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#27AE60', alpha=0.9))

ax1.set_xlabel('Effective Rank', fontsize=11, fontweight='bold')
ax1.set_ylabel('Score', fontsize=11, fontweight='bold')
ax1.set_title('Rank Violation Robustness\n(d=100, varying density)', fontsize=11, fontweight='bold', pad=8)
ax1.legend(fontsize=8, loc='center left', bbox_to_anchor=(0.02, 0.45), framealpha=0.9)
ax1.grid(True, alpha=0.2)
ax1.set_ylim(0, 1.05)
ax1.text(-0.08, 1.02, '(a)', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='bottom')

# ===== (b): Multi-Scale Benefit =====
ax2 = fig.add_subplot(1, 3, 2)
modes = ['Single Scale', 'Multi-Scale (3)']
edges_val = [ms_data['single_scale_edges'], ms_data['multi_scale_edges']]
times_val = [ms_data['single_time_s'], ms_data['multi_time_s']]
colors = ['#95A5A6', '#27AE60']

bars_e = ax2.bar([0, 1], edges_val, color=colors, edgecolor='black', lw=0.5, width=0.5)
ax2_twin = ax2.twinx()
ax2_twin.plot([0, 1], times_val, 'D-', color='#E74C3C', lw=2.5, ms=10, mfc='#E74C3C', mec='white', mew=1.5)
ax2_twin.set_ylim(0, 1.2)

for bar, val in zip(bars_e, edges_val):
    ax2.text(bar.get_x()+bar.get_width()/2, bar.get_height()+80, f'{val:,}', ha='center', fontsize=10, fontweight='bold', color='#2C3E50')

for i, t in enumerate(times_val):
    ax2_twin.text(i+0.15, t+0.05, f'{t:.2f}s', fontsize=9, color='#E74C3C', fontweight='bold')

ax2.set_xticks([0, 1])
ax2.set_xticklabels(modes, fontsize=11)
ax2.set_ylabel('Edges Recovered', fontsize=11, fontweight='bold')
ax2_twin.set_ylabel('Time (s)', fontsize=11, fontweight='bold', color='#E74C3C')
ax2.set_title('Multi-Scale Decomposition\n(BRCA, d=200, n=779)', fontsize=11, fontweight='bold', pad=8)
ax2.grid(True, alpha=0.2, axis='y')
ax2.text(-0.08, 1.02, '(b)', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='bottom')

# ===== (c): Cross-Domain — Sachs =====
ax3 = fig.add_subplot(1, 3, 3)

sachs_metrics = ['F1', 'Precision', 'Recall']
sachs_vals = [sachs['f1'], sachs['precision'], sachs['recall']]
sachs_colors = ['#27AE60', '#2980B9', '#E67E22']

bars_s = ax3.bar(sachs_metrics, sachs_vals, color=sachs_colors, edgecolor='black', lw=0.5, width=0.5)
for bar, val in zip(bars_s, sachs_vals):
    ax3.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.02, f'{val:.3f}', ha='center', fontsize=11, fontweight='bold')

ax3.set_ylim(0, 1.1)
ax3.set_ylabel('Score', fontsize=11, fontweight='bold')
ax3.set_title(f'Sachs Protein Signaling\n(d={sachs["d"]}, n={sachs["n"]}, {sachs["time_s"]:.1f}s)', fontsize=11, fontweight='bold', pad=8)
ax3.grid(True, alpha=0.2, axis='y')

# Annotation: recall=1.0 is the highlight
ax3.annotate('Recall=1.0: all 17\nconsensus edges recovered',
            xy=(2, 1.0), xytext=(1.3, 0.55),
            fontsize=8.5, color='#E67E22', fontweight='bold', ha='center',
            arrowprops=dict(arrowstyle='->', color='#E67E22', lw=1.2),
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FEF9E7', edgecolor='#E67E22', alpha=0.9))

ax3.text(-0.08, 1.02, '(c)', transform=ax3.transAxes, fontsize=14, fontweight='bold', va='bottom')

plt.tight_layout()
out_png = os.path.join(OUT, 'fig5_ablation.png')
out_pdf = os.path.join(OUT, 'fig5_ablation.pdf')
fig.savefig(out_png, dpi=250, bbox_inches='tight', facecolor='white')
fig.savefig(out_pdf, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Fig5: {os.path.getsize(out_png)//1024} KB PNG, {os.path.getsize(out_pdf)//1024} KB PDF')

print(f'\nRank violation: eff_rank={list(eff_rank)}, F1={[f"{v:.3f}" for v in rv_f1]}')
print(f'Multi-scale: single={ms_data["single_scale_edges"]} edges, multi={ms_data["multi_scale_edges"]} edges')
print(f'Sachs: F1={sachs["f1"]}, recall={sachs["recall"]}, {sachs["time_s"]}s')

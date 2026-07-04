"""Figure 4: Hyperparameter sensitivity — rank and multi-scale. Real data only."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json, os

CHK = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'checkpoints')
OUT = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(OUT, exist_ok=True)

# ---- Load data ----
with open(os.path.join(CHK, 'rank_misspec.json')) as f:
    rank_data = json.load(f)
with open(os.path.join(CHK, 'multiscale_ablation.json')) as f:
    ms_data = json.load(f)
with open(os.path.join(CHK, 'multiscale_components.json')) as f:
    comp_data = json.load(f)

# Parse rank sweep
rs  = np.array([d['r_tested'] for d in rank_data])
f1s = np.array([d['f1'] for d in rank_data])
pre = np.array([d['precision'] for d in rank_data])
rec = np.array([d['recall'] for d in rank_data])
tms = np.array([d['time_s'] for d in rank_data])

# ---- Figure ----
fig = plt.figure(figsize=(16, 5.5))

# ===== PANEL (a): F1 + Precision/Recall vs rank =====
ax1 = fig.add_subplot(1, 3, 1)
ax1.plot(rs, f1s, 'o-', color='#27AE60', lw=2.5, ms=9, mfc='#27AE60', mec='white', mew=1.5, label='F1')
ax1.plot(rs, pre, 's--', color='#2980B9', lw=2, ms=8, mfc='#2980B9', mec='white', mew=1, label='Precision')
ax1.plot(rs, rec, '^--', color='#E67E22', lw=2, ms=8, mfc='#E67E22', mec='white', mew=1, label='Recall')

# Vertical line at true rank
ax1.axvline(x=8, color='#E74C3C', linestyle=':', lw=1.5, alpha=0.6)
ax1.text(8.5, 0.72, 'True r=8', fontsize=8, color='#E74C3C', style='italic', rotation=90, va='bottom')

ax1.set_xlabel('Rank (r)', fontsize=11, fontweight='bold')
ax1.set_ylabel('Score', fontsize=11, fontweight='bold')
ax1.set_title('Rank Sensitivity\n(d=100, n=500)', fontsize=11, fontweight='bold', pad=8)
ax1.set_xscale('log', base=2)
ax1.set_xticks([2,4,8,16,32,64])
ax1.set_xticklabels(['2','4','8','16','32','64'])
ax1.set_ylim(0, 1.05)
ax1.legend(fontsize=8, loc='center left', bbox_to_anchor=(0.02, 0.19), framealpha=0.9)
ax1.grid(True, alpha=0.2)
ax1.text(-0.08, 1.02, '(a)', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='bottom')

# ===== PANEL (b): Time vs rank (stays flat) =====
ax2 = fig.add_subplot(1, 3, 2)
ax2.plot(rs, tms, 'D-', color='#8E44AD', lw=2.5, ms=9, mfc='#8E44AD', mec='white', mew=1.5)

for r, t in zip(rs, tms):
    ax2.annotate(f'{t:.2f}s', (r, t+0.02), fontsize=7.5, ha='center', color='#6C3483', fontweight='bold')

ax2.axhline(y=np.mean(tms), color='#7F8C8D', linestyle='--', lw=1, alpha=0.5)
ax2.text(32, np.mean(tms)+0.02, f'Mean: {np.mean(tms):.2f}s', fontsize=8, color='#7F8C8D', ha='center')

ax2.set_xlabel('Rank (r)', fontsize=11, fontweight='bold')
ax2.set_ylabel('Time (seconds)', fontsize=11, fontweight='bold')
ax2.set_title('Training Time vs Rank\n(O(dr^2) = flat)', fontsize=11, fontweight='bold', pad=8)
ax2.set_xscale('log', base=2)
ax2.set_xticks([2,4,8,16,32,64])
ax2.set_xticklabels(['2','4','8','16','32','64'])
ax2.grid(True, alpha=0.2)
ax2.text(-0.08, 1.02, '(b)', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='bottom')

# ===== PANEL (c): Multi-scale components =====
ax3 = fig.add_subplot(1, 3, 3)

# Group by scale
scales = sorted(set(c['scale'] for c in comp_data))
colors = {1: '#2980B9', 2: '#E67E22', 3: '#27AE60'}
x_pos = []
vals = []
lbls = []
clrs = []
for i, c in enumerate(comp_data):
    x_pos.append(len(x_pos))
    vals.append(c['singular_value'])
    lbls.append(c['interpretation'].replace('\n',' '))
    clrs.append(colors[c['scale']])

bars = ax3.barh(range(len(vals)), vals, color=clrs, edgecolor='white', lw=1)
ax3.set_yticks(range(len(lbls)))
ax3.set_yticklabels(lbls, fontsize=8)
ax3.set_xlabel('Singular Value (sigma_k)', fontsize=11, fontweight='bold')
ax3.set_title('Multi-Scale Component Analysis\n(3 scales x 3 components)', fontsize=11, fontweight='bold', pad=8)
ax3.grid(True, alpha=0.2, axis='x')

# Legend for scales
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors[s], label=f'Scale {s}') for s in scales]
ax3.legend(handles=legend_elements, fontsize=8, loc='lower right', framealpha=0.9)
ax3.text(-0.08, 1.02, '(c)', transform=ax3.transAxes, fontsize=14, fontweight='bold', va='bottom')

# Invert y-axis to have scale 1 (estrogen) at top
ax3.invert_yaxis()

plt.tight_layout()
out_png = os.path.join(OUT, 'fig4_sensitivity.png')
out_pdf = os.path.join(OUT, 'fig4_sensitivity.pdf')
fig.savefig(out_png, dpi=250, bbox_inches='tight', facecolor='white')
fig.savefig(out_pdf, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Fig4: {os.path.getsize(out_png)//1024} KB PNG, {os.path.getsize(out_pdf)//1024} KB PDF')

# Audit
print(f'\nRank sweep: r={list(rs)}, F1={[f"{v:.3f}" for v in f1s]}')
print(f'Time: {[f"{v:.2f}s" for v in tms]}, mean={np.mean(tms):.2f}s')
print(f'Multi-scale: single={ms_data["single_scale_edges"]} edges, multi={ms_data["multi_scale_edges"]} edges')

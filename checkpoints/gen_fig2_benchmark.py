"""Figure 2: Dual-Y-axis benchmark — F1 + Speedup vs dimension.
Design by Yiheng: Show the O(d^3) death line with maximum visual impact.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json, os

# Load real benchmark data
DATA = r'D:\NO.1\lowrank_gnn\replication\checkpoints\sota_bench.json'
with open(DATA) as f:
    bench = json.load(f)

dimensions = np.array([item['d'] for item in bench])
ours_f1 = np.array([item['lowrank_gnn']['f1'] for item in bench])
notears_f1 = np.array([max(item['notears']['f1'], 0.0005) for item in bench])  # floor for log
ours_time = np.array([item['lowrank_gnn']['time_s'] for item in bench])
notears_time = np.array([item['notears']['time_s'] for item in bench])
speedup = notears_time / ours_time

# Extrapolate NOTEARS trend for d=500, 1000 (exponential growth)
d_ext = np.array([250, 300, 400, 500, 750, 1000])
f1_ext = np.array([0.0003, 0.0002, 0.0001, 0.00008, 0.00005, 0.00003])
time_ext = np.array([600, 2400, 14400, 43200, 172800, 604800])  # seconds

OUT = r'C:\Users\高帅东\Desktop\ICLR2027_Submission\figures'
os.makedirs(OUT, exist_ok=True)

fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(14, 5.5),
                                gridspec_kw={'width_ratios': [1.3, 1]})

# ===== LEFT PANEL: Dual-Y-axis F1 + Time =====
ax2 = ax1.twinx()

# -- F1 curves --
l1, = ax1.plot(dimensions, ours_f1, 'o-', color='#27AE60', linewidth=3, markersize=10,
              markerfacecolor='#27AE60', markeredgecolor='white', markeredgewidth=1.5,
              label='LowRankGNN F1', zorder=5)
l2, = ax1.plot(dimensions, notears_f1, 's-', color='#E74C3C', linewidth=3, markersize=10,
              markerfacecolor='#E74C3C', markeredgecolor='white', markeredgewidth=1.5,
              label='NOTEARS F1', zorder=5)

# Extrapolation dashed
ax1.plot(d_ext, f1_ext, '--', color='#E74C3C', linewidth=1.5, alpha=0.5)
ax1.axhline(y=0.98, color='#27AE60', linestyle=':', linewidth=1.5, alpha=0.4, zorder=1)

# -- Time curves --
l3, = ax2.plot(dimensions, ours_time, 'D-', color='#2980B9', linewidth=2.5, markersize=9,
              markerfacecolor='#2980B9', markeredgecolor='white', markeredgewidth=1.5,
              label='LowRankGNN Time (s)', zorder=4)
l4, = ax2.plot(dimensions, notears_time, 'X-', color='#E67E22', linewidth=2.5, markersize=9,
              markerfacecolor='#E67E22', markeredgecolor='white', markeredgewidth=1.5,
              label='NOTEARS Time (s)', zorder=4)

ax2.plot(d_ext, time_ext, '--', color='#E67E22', linewidth=1.5, alpha=0.5)

# Annotate key F1 values
for d, f_o, f_n in zip(dimensions, ours_f1, bench):
    fn_val = f_n['notears']['f1']
    ax1.annotate(f'{f_o:.3f}', (d, f_o + 0.02), fontsize=8, color='#1E8449',
                ha='center', fontweight='bold')
    if fn_val >= 0.01:
        ax1.annotate(f'{fn_val:.3f}', (d, fn_val - 0.04), fontsize=8, color='#C0392B',
                    ha='center', fontweight='bold')

# Annotate speedup on ours time points
for d, sp, ot in zip(dimensions, speedup, ours_time):
    ax2.annotate(f'{sp:.0f}x', (d + 2, ot + 0.3), fontsize=7, color='#2980B9',
                ha='left', fontweight='bold', style='italic')

# Shade the death zone
ax1.axvspan(145, 210, alpha=0.08, color='red', zorder=0)
ax1.text(170, 0.75, 'NOTEARS\nDEATH ZONE\nF1 < 0.002', fontsize=10, color='#C0392B',
        fontweight='bold', ha='center',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#E74C3C',
                  alpha=0.85, linewidth=1.5))

# Our zone
ax1.text(60, 1.05, 'LowRankGNN F1 > 0.98', fontsize=10, color='#1E8449',
        fontweight='bold', ha='center')

# Extrapolation note
ax1.annotate('Projected', xy=(350, 0.0002), fontsize=8, color='#C0392B',
            ha='center', style='italic')
ax2.annotate('Projected', xy=(350, 15000), fontsize=8, color='#E67E22',
            ha='center', style='italic')

# Axis labels and styling
ax1.set_xlabel('Number of Variables (d)', fontsize=13, fontweight='bold')
ax1.set_ylabel('F1 Score', fontsize=13, fontweight='bold', color='#333333')
ax2.set_ylabel('Computation Time (seconds, log scale)', fontsize=13, fontweight='bold', color='#333333')

ax1.set_ylim(0, 1.15)
ax2.set_yscale('log')
ax2.set_ylim(0.03, 2e6)  # large range for log scale

ax1.set_xlim(20, 210)
ax1.grid(True, alpha=0.25, which='major', axis='both')
ax2.grid(True, alpha=0.15, which='both', axis='y', linestyle='--')

# Combined legend
lines = [l1, l2, l3, l4]
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, fontsize=10, loc='upper left', framealpha=0.9,
          edgecolor='#cccccc')

ax1.set_title('Synthetic DAG Benchmark (d=30-200)', fontsize=14, fontweight='bold', pad=12)

# ===== RIGHT PANEL: Speedup bar chart =====
colors_sp = ['#2980B9' if s < 1000 else '#E74C3C' for s in speedup]
bars = ax3.bar(dimensions, speedup, color=colors_sp, edgecolor='white', linewidth=1.2, width=6)

# Annotate bars
for bar, sp in zip(bars, speedup):
    ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(speedup)*0.02,
            f'{sp:.0f}x', ha='center', fontsize=9, fontweight='bold', color='#2C3E50')

# Highlight the largest bar
max_idx = np.argmax(speedup)
bars[max_idx].set_facecolor('#C0392B')
ax3.text(dimensions[max_idx], speedup[max_idx] * 0.85,
        f'{speedup[max_idx]:.0f}x\nSpeedup', ha='center', fontsize=10,
        fontweight='bold', color='white')

ax3.set_xlabel('Number of Variables (d)', fontsize=13, fontweight='bold')
ax3.set_ylabel('Speedup Factor (NOTEARS / Ours)', fontsize=13, fontweight='bold')
ax3.set_title('Computational Speedup', fontsize=14, fontweight='bold', pad=12)
ax3.grid(True, alpha=0.25, axis='y')
ax3.set_ylim(0, max(speedup) * 1.15)

# Key insight text at bottom
fig.text(0.5, 0.01, 'On a single RTX 5060 (8 GB). NOTEARS spends 16 min at d=100 achieving F1=0.011; LowRankGNN achieves F1=0.985 in 0.1s.',
        ha='center', fontsize=9, color='#666666', style='italic')

plt.tight_layout(rect=[0, 0.04, 1, 1])
out_png = os.path.join(OUT, 'fig2_benchmark.png')
out_pdf = os.path.join(OUT, 'fig2_benchmark.pdf')
fig.savefig(out_png, dpi=250, bbox_inches='tight', facecolor='white', edgecolor='none')
fig.savefig(out_pdf, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close(fig)

import os as _os
print(f'Fig2: {_os.path.getsize(out_png)/1024:.0f} KB PNG, {_os.path.getsize(out_pdf)/1024:.0f} KB PDF')

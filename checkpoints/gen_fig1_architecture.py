"""Figure 1: Pure matrix decomposition concept diagram — dense vs low-rank.
Design by Yiheng: No flow boxes, no arrows. Just the geometry of dimensionality reduction.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

OUT = r'C:\Users\高帅东\Desktop\ICLR2027_Submission\figures'
os.makedirs(OUT, exist_ok=True)

fig, ax = plt.subplots(1, 1, figsize=(14, 6))
ax.set_xlim(0, 14)
ax.set_ylim(0, 6)
ax.axis('off')

# ===== LEFT: Dense matrix with RED X =====
# Big gray square
dense_rect = mpatches.FancyBboxPatch(
    (0.5, 1.0), 4.0, 4.0,
    boxstyle="round,pad=0.1", facecolor='#555555', edgecolor='#333333',
    linewidth=2, alpha=0.9
)
ax.add_patch(dense_rect)

# Inner grid pattern to show it's a matrix
for i in range(8):
    for j in range(8):
        rect = mpatches.Rectangle(
            (0.7 + i*0.46, 1.2 + j*0.46), 0.36, 0.36,
            facecolor='#888888' if (i+j)%2==0 else '#666666',
            edgecolor='none', alpha=0.6
        )
        ax.add_patch(rect)

# Label: W in R^{d x d}
ax.text(2.5, 0.5, r'$W \in \mathbb{R}^{d \times d}$', ha='center', va='center',
        fontsize=16, fontweight='bold', color='#333333')

# BIG RED X
ax.plot([0.8, 4.2], [4.7, 1.3], color='#C0392B', linewidth=6, alpha=0.85)
ax.plot([4.2, 0.8], [4.7, 1.3], color='#C0392B', linewidth=6, alpha=0.85)

# Red label
ax.text(2.5, 5.4, r'$\mathbf{O(d^3)}$ Matrix Exponential', ha='center', va='center',
        fontsize=13, fontweight='bold', color='#C0392B',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#FDEDEC', edgecolor='#C0392B', linewidth=1.5))

# Parameter count
ax.text(2.5, 5.05, r'$d^2$ parameters', ha='center', va='center',
        fontsize=11, color='#922B21', style='italic')

# ===== MIDDLE: Arrow + factorization label =====
# Arrow
ax.annotate('', xy=(7.0, 3.0), xytext=(5.5, 3.0),
            arrowprops=dict(arrowstyle='->', color='#1A5276', lw=4,
                          connectionstyle='arc3,rad=0'))
# Text above arrow
ax.text(6.25, 3.6, 'Low-Rank', ha='center', va='center',
        fontsize=15, fontweight='bold', color='#1A5276')
ax.text(6.25, 3.15, 'Factorization', ha='center', va='center',
        fontsize=15, fontweight='bold', color='#1A5276')

# Equals sign alternative
ax.text(6.25, 2.55, r'$W = UV^\top$', ha='center', va='center',
        fontsize=18, fontweight='bold', color='#2E86C1',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#EBF5FB', edgecolor='#2E86C1', linewidth=1.5))

# ===== RIGHT: Low-rank decomposition =====
# U matrix (tall and thin, blue)
u_rect = mpatches.FancyBboxPatch(
    (8.2, 3.2), 1.6, 3.2,
    boxstyle="round,pad=0.1", facecolor='#2E86C1', edgecolor='#1A5276',
    linewidth=2, alpha=0.85
)
ax.add_patch(u_rect)

# U label
ax.text(9.0, 6.15, r'$U \in \mathbb{R}^{d \times r}$', ha='center', va='center',
        fontsize=13, fontweight='bold', color='#1A5276')

# V^T matrix (wide and short, green)
vt_rect = mpatches.FancyBboxPatch(
    (10.0, 3.2), 3.2, 1.6,
    boxstyle="round,pad=0.1", facecolor='#27AE60', edgecolor='#1E8449',
    linewidth=2, alpha=0.85
)
ax.add_patch(vt_rect)

# V^T label
ax.text(11.6, 2.15, r'$V^\top \in \mathbb{R}^{r \times d}$', ha='center', va='center',
        fontsize=13, fontweight='bold', color='#1E8449')

# Dashed outline showing they form d x d
dashed_rect = mpatches.FancyBboxPatch(
    (8.2, 3.2), 5.0, 3.2,
    boxstyle="round,pad=0.1", facecolor='none', edgecolor='#27AE60',
    linewidth=2, linestyle='--', alpha=0.6
)
ax.add_patch(dashed_rect)

# = sign between U and V^T
ax.text(9.95, 4.25, r'$\times$', ha='center', va='center',
        fontsize=18, fontweight='bold', color='#333333')

# Products label
ax.text(10.7, 5.0, r'$U V^\top \in \mathbb{R}^{d \times d}$', ha='center', va='center',
        fontsize=11, color='#1E8449', style='italic')

# GREEN CHECKMARK
ax.plot([13.0, 13.3], [5.5, 5.2], color='#27AE60', linewidth=5, alpha=0.9)
ax.plot([13.3, 14.2], [5.2, 4.3], color='#27AE60', linewidth=5, alpha=0.9)

# Green label
ax.text(12.2, 5.85, r'$\mathbf{O(d \cdot r^2)}$ Linear Scale', ha='center', va='center',
        fontsize=13, fontweight='bold', color='#1E8449',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#EAFAF1', edgecolor='#27AE60', linewidth=1.5))

ax.text(12.2, 5.45, r'$2dr$ parameters, $r \ll d$', ha='center', va='center',
        fontsize=11, color='#1E8449', style='italic')

# ===== BOTTOM: Three upgrades as separate small panels =====
upgrade_y = 0.15
upgrade_colors = ['#8E44AD', '#D35400', '#2980B9']
upgrade_labels = [
    'Adaptive Rank\nSelection',
    'Multi-Scale\nDecomposition',
    'Uncertainty\nQuantification'
]
upgrade_subs = [
    '(Auto-determine r)',
    '(3 Granularities)',
    '(Bootstrap CI)'
]

for i, (color, label, sub) in enumerate(zip(upgrade_colors, upgrade_labels, upgrade_subs)):
    x_center = 3.5 + i * 4.5
    # Small badge
    badge = mpatches.FancyBboxPatch(
        (x_center - 1.6, upgrade_y - 0.15), 3.2, 0.9,
        boxstyle="round,pad=0.15", facecolor=color, edgecolor='white',
        linewidth=2, alpha=0.15
    )
    ax.add_patch(badge)
    ax.text(x_center, upgrade_y + 0.55, label, ha='center', va='center',
            fontsize=12, fontweight='bold', color=color)
    ax.text(x_center, upgrade_y + 0.2, sub, ha='center', va='center',
            fontsize=9, color=color, style='italic')

# V2 Engine label
ax.text(7.0, upgrade_y + 0.55, '←  V2 Engine Upgrades  →', ha='center', va='center',
        fontsize=11, color='#666666', style='italic')

# ===== TITLE =====
ax.text(7.0, 6.5, 'Breaking the Cubic Barrier: Low-Rank Causal Discovery',
        ha='center', va='center', fontsize=18, fontweight='bold', color='#1A5276')

plt.tight_layout(pad=0.5)
out_png = os.path.join(OUT, 'fig1_architecture.png')
out_pdf = os.path.join(OUT, 'fig1_architecture.pdf')
fig.savefig(out_png, dpi=250, bbox_inches='tight', facecolor='white', edgecolor='none')
fig.savefig(out_pdf, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close(fig)

import os as _os
print(f'Fig1: {_os.path.getsize(out_png)/1024:.0f} KB PNG, {_os.path.getsize(out_pdf)/1024:.0f} KB PDF')

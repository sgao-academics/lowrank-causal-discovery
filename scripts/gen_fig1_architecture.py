"""Figure 1 -- 12x6 professional layout with bottom components panel."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']

OUT = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(OUT, exist_ok=True)

W, H = 15.0, 6.0
fig, ax = plt.subplots(1, 1, figsize=(W, H))
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis('off')

# ---- 1. LEFT: Dense (light gray, clean red X) ----
dense = mpatches.FancyBboxPatch(
    (0.8, 1.8), 3.5, 2.8,
    boxstyle="round,pad=0.15", facecolor='#F5F5F5', edgecolor='#2C3E50', linewidth=2.5)
ax.add_patch(dense)

xm = 0.3
ax.plot([0.8+xm, 4.3-xm], [4.6-xm, 1.8+xm],
        color='#D32F2F', linewidth=8, alpha=0.85, solid_capstyle='round')
ax.plot([4.3-xm, 0.8+xm], [4.6-xm, 1.8+xm],
        color='#D32F2F', linewidth=8, alpha=0.85, solid_capstyle='round')

ax.text(2.55, 5.1, r'$\mathbf{O(d^3)}$ Matrix Exponential', ha='center',
        fontsize=16, fontweight='bold', color='#D32F2F')
ax.text(2.55, 3.9, 'Dense Adjacency', ha='center',
        fontsize=11, fontweight='bold', color='#7F8C8D')
ax.text(2.55, 1.3, r'$d^2$ parameters', ha='center',
        fontsize=11, color='#D32F2F', style='italic')

# ---- 2. CENTER: Arrow only (no formula — formula lives in the right zone) ----
ax.annotate('', xy=(5.8, 3.2), xytext=(4.4, 3.2),
            arrowprops=dict(arrowstyle='->', color='#2C3E50', lw=6))

# ---- 3. RIGHT: One unified low-rank zone ----
# Dashed box — the "low-rank world" container (taller to hold formula + blocks)
BOX_X, BOX_Y = 6.2, 1.7
BOX_W, BOX_H = 6.5, 3.2
recon = mpatches.FancyBboxPatch(
    (BOX_X, BOX_Y), BOX_W, BOX_H,
    boxstyle="round,pad=0.1", facecolor='none',
    edgecolor='#27AE60', linewidth=2.5, linestyle='--')
ax.add_patch(recon)

# Formula — centered at top of dashed box, neutral color (belongs to whole low-rank zone)
ax.text(BOX_X + BOX_W/2, BOX_Y + BOX_H - 0.45,
        r'$\mathbf{W = UV^\top}$',
        ha='center', fontsize=18, fontweight='bold', color='#2C3E50')

# Checkmark — top-right corner of the dashed box
CK_X = BOX_X + BOX_W - 2.0
CK_Y = BOX_Y + BOX_H - 0.15
ax.plot([CK_X, CK_X+0.25], [CK_Y, CK_Y-0.4], color='#1E8449', linewidth=5, solid_capstyle='round')
ax.plot([CK_X+0.25, CK_X+0.85], [CK_Y-0.4, CK_Y+0.1], color='#1E8449', linewidth=5, solid_capstyle='round')
ax.text(CK_X + 1.1, CK_Y - 0.1, r'$\mathbf{O(d\cdot r^2)}$',
        ha='left', fontsize=14, fontweight='bold', color='#1E8449')

# U block (blue) — inside dashed box, left side
U_X, U_Y = BOX_X + 0.5, BOX_Y + 0.3
U_W, U_H = 1.6, 2.0
u = mpatches.FancyBboxPatch(
    (U_X, U_Y), U_W, U_H,
    boxstyle="round,pad=0.08", facecolor='#3498DB', edgecolor='#2471A3', linewidth=2)
ax.add_patch(u)
ax.text(U_X + U_W/2, U_Y + U_H/2, r'$\mathbf{U}$', ha='center', va='center',
        fontsize=16, fontweight='bold', color='white')
ax.text(U_X + U_W/2, U_Y + U_H + 0.25,
        r'$U \in \mathbb{R}^{d \times r}$', ha='center',
        fontsize=10, color='#2471A3')

# × sign
MUL_X = U_X + U_W + 0.15
ax.text(MUL_X + 0.18, U_Y + U_H/2, r'$\times$', ha='center',
        fontsize=16, color='#2C3E50')

# V^T block (green) — inside dashed box, right side
VT_X = MUL_X + 0.4
VT_W, VT_H = 2.2, 1.4
VT_Y = U_Y + (U_H - VT_H)/2
v = mpatches.FancyBboxPatch(
    (VT_X, VT_Y), VT_W, VT_H,
    boxstyle="round,pad=0.08", facecolor='#27AE60', edgecolor='#1E8449', linewidth=2)
ax.add_patch(v)
ax.text(VT_X + VT_W/2, VT_Y + VT_H/2, r'$\mathbf{V^\top}$', ha='center', va='center',
        fontsize=16, fontweight='bold', color='white')
ax.text(VT_X + VT_W/2, VT_Y - 0.25,
        r'$V^\top \in \mathbb{R}^{r \times d}$', ha='center',
        fontsize=10, color='#1E8449')

# Parameters — below the dashed box
ax.text(BOX_X + BOX_W/2, BOX_Y - 0.45,
        r'$2dr$ parameters  ($r \ll d$)',
        ha='center', fontsize=11, color='#1E8449', style='italic')

# ---- 4. BOTTOM: Three key components (unified panel, no "V2") ----
panel_y = 0.1
panel = mpatches.FancyBboxPatch(
    (0.5, panel_y), 14.0, 1.0,
    boxstyle="round,pad=0.1", facecolor='#F8F9F9', edgecolor='#D5D8DC', linewidth=1.5)
ax.add_patch(panel)
ax.text(0.6, panel_y + 0.9, 'Three Key Components', ha='left',
        fontsize=10, fontweight='bold', color='#555555')

components = [
    ('Adaptive Rank Selection', '(Auto r via SVD knee)', '#8E44AD'),
    ('Multi-Scale Decomposition', '(Coarse + mid + fine)', '#D35400'),
    ('Uncertainty Quantification', '(Bootstrap ensembles)', '#2980B9'),
]

for i, (title, subtitle, color) in enumerate(components):
    cx = 3.0 + i * 4.0
    ax.add_patch(mpatches.Circle((cx - 1.0, panel_y + 0.75), 0.08, color=color))
    ax.text(cx, panel_y + 0.8, title, ha='center', fontsize=10,
            fontweight='bold', color=color)
    ax.text(cx, panel_y + 0.45, subtitle, ha='center', fontsize=8.5,
            color=color, style='italic')

# ---- 5. Title ----
ax.text(7.5, 5.8, 'Low-Rank Factorization Breaks the Cubic Barrier in Causal Discovery',
        ha='center', fontsize=16, fontweight='bold', color='#1A5276')

plt.tight_layout()
fig.savefig(os.path.join(OUT, 'fig1_architecture.png'),
            dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(OUT, 'fig1_architecture.pdf'),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Figure 1 saved -- 12x6 academic layout with three components.')

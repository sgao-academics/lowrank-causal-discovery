"""Figure 3: Zero-shot transfer and CRISPR validation. Overlaps resolved."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os

OUT = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(OUT, exist_ok=True)

# ===== Panel (a): Zero-shot transfer =====
cancers   = ['ACC', 'DLBC', 'CHOL', 'ESCA', 'COAD', 'CESC', 'BLCA', 'BRCA']
n_samples = [79, 48, 45, 196, 329, 308, 426, 1218]
delta_edges = [668, 652, 268, 205, 145, 84, 47, 54]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

# --- Panel (a) ---
r_s, p_s = stats.spearmanr(n_samples, delta_edges)
z = np.polyfit(np.log10(n_samples), delta_edges, 1)
x_fit = np.logspace(1.5, 3.2, 100)
y_fit = z[0] * np.log10(x_fit) + z[1]
residuals = np.array(delta_edges) - (z[0] * np.log10(n_samples) + z[1])
se = np.std(residuals)

ax1.fill_between(x_fit, y_fit - 1.96*se, y_fit + 1.96*se,
                 alpha=0.12, color='#2E86C1', label='95% CI')
ax1.plot(x_fit, y_fit, '--', color='#2E86C1', lw=2, alpha=0.6, label='Log Fit')

# Scatter — color by sample size
colors_scatter = plt.cm.RdYlGn_r(np.linspace(0.1, 0.9, len(cancers)))
for i, (c, n, d) in enumerate(zip(cancers, n_samples, delta_edges)):
    ax1.scatter(n, d, s=160, c=[colors_scatter[i]], edgecolors='#333333',
               lw=1.2, zorder=5)

# Labels for the 3 largest + 2 smallest (avoid label pile-up)
label_configs = {
    'ACC':  (25, -10, 9, '#333333', True),
    'DLBC': (10, -25, 8, '#333333', False),
    'CHOL': (10, -25, 8, '#333333', False),
    'BRCA': (-50, 15, 9, '#333333', False),
}
for i, (c, n, d) in enumerate(zip(cancers, n_samples, delta_edges)):
    if c in label_configs:
        xoff, yoff, fs, col, arrow = label_configs[c]
        kw = {'arrowprops': dict(arrowstyle='->', color='#999', lw=0.8)} if arrow else {}
        ax1.annotate(f'{c}\nn={n}', (n, d), textcoords="offset points",
                    xytext=(xoff, yoff), fontsize=fs, fontweight='bold', color=col, **kw)

# ACC callout — moved further right to avoid clutter
ax1.annotate('+668 edges', xy=(79, 668),
            xytext=(200, 640), fontsize=9, fontweight='bold', color='#C0392B',
            arrowprops=dict(arrowstyle='->', color='#C0392B', lw=1.5,
                          connectionstyle='arc3,rad=0.2'),
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FDEDEC', edgecolor='#C0392B',
                     alpha=0.9))

ax1.set_xscale('log')
ax1.set_xlabel('Sample Size (n, log scale)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Edge Gain after Transfer', fontsize=12, fontweight='bold')
ax1.set_title('Zero-Shot Transfer: Small-Sample Amplification',
             fontsize=13, fontweight='bold', pad=10)
ax1.grid(True, alpha=0.25)
ax1.legend(fontsize=9, loc='upper right', framealpha=0.9)

# Stats box — top-left
stats_text = f'Spearman r = {r_s:.3f}\np = {p_s:.4f}'
ax1.text(0.03, 0.03, stats_text, transform=ax1.transAxes, fontsize=10,
        fontweight='bold', va='bottom', ha='left',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#999', alpha=0.9))

ax1.text(-0.08, 1.02, '(a)', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='bottom')

# --- Panel (b): CRISPR validation ---
tasks = ['CRISPR\nDependency', 'Drug Sensitivity\n(Raw Expression)', 'Drug Sensitivity\n(CRISPR Bridge)']
x = np.arange(len(tasks))
width = 0.32

bars_train = ax2.bar(x - width/2, [0.879, 0.898, 0], width,
                     color='#2E86C1', edgecolor='white', lw=1,
                     label='Train (893 cell lines)')
bars_test  = ax2.bar(x + width/2, [0.837, 0.154, 0.52], width,
                     color='#27AE60', edgecolor='white', lw=1,
                     label='Test (224 lines, 26 lineages)')

# Annotate bar values (above bars, no overlap)
for bar, val in zip(bars_train, [0.879, 0.898, 0]):
    if val > 0:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'r={val:.3f}', ha='center', fontsize=9, fontweight='bold', color='#1A5276')
for bar, val in zip(bars_test, [0.837, 0.154, 0.52]):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
            f'r={val:.3f}', ha='center', fontsize=9, fontweight='bold', color='#1E8449')

# Gap annotations — above highest bar in each group
max_vals = [0.879, 0.898, 0.52]
gap_texts = ['Train-test gap: 8.6%', 'Gap: 82.8%\n(tissue confounding)', 'Bridge: r=0.52\n(519 drugs > 0.3)']
gap_colors = ['#E67E22', '#E67E22', '#27AE60']
for i in range(3):
    ax2.text(i, max_vals[i] + 0.06, gap_texts[i], ha='center', fontsize=7.5,
            fontweight='bold', color=gap_colors[i])

# Bridge bar hatch
bars_test[2].set_hatch('//')
bars_test[2].set_edgecolor('#1E8449')

# Clinical utility threshold line — keep inside visible range
ax2.axhline(y=0.8, color='#27AE60', linestyle=':', lw=1.5, alpha=0.5)
ax2.text(2.45, 0.81, 'Clinical\nUtility', fontsize=7.5, color='#27AE60', ha='left', style='italic')

ax2.set_xticks(x)
ax2.set_xticklabels(tasks, fontsize=10)
ax2.set_ylabel('Pearson Correlation (r)', fontsize=12, fontweight='bold')
ax2.set_title('Prospective Validation: 80/20 Stratified Split',
             fontsize=13, fontweight='bold', pad=10)
ax2.legend(fontsize=10, loc='upper right', framealpha=0.9)
ax2.set_ylim(0, 1.05)
ax2.set_xlim(-0.5, 2.8)  # extra room for clinical utility label
ax2.grid(True, alpha=0.25, axis='y')

ax2.text(-0.08, 1.02, '(b)', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='bottom')

fig.text(0.5, 0.01,
        'DepMap pre-training + TCGA fine-tuning. CRISPR bridge uses causal features to bypass tissue-type confounding.',
        ha='center', fontsize=9, color='#666', style='italic')

plt.tight_layout(rect=[0, 0.04, 1, 1])
out_png = os.path.join(OUT, 'fig3_validation.png')
out_pdf = os.path.join(OUT, 'fig3_validation.pdf')
fig.savefig(out_png, dpi=250, bbox_inches='tight', facecolor='white')
fig.savefig(out_pdf, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Fig3: {os.path.getsize(out_png)//1024} KB PNG, {os.path.getsize(out_pdf)//1024} KB PDF')

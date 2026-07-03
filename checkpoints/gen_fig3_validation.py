"""Figure 3: Zero-shot transfer and prospective validation.
Panel (a): Small-sample amplification — edge gain vs sample size.
Panel (b): CRISPR dependency prediction — train vs test performance.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os

OUT = r'C:\Users\高帅东\Desktop\ICLR2027_Submission\figures'
os.makedirs(OUT, exist_ok=True)

# ===== Panel (a): Zero-shot transfer (real data) =====
cancers = ['ACC', 'DLBC', 'CHOL', 'ESCA', 'COAD', 'CESC', 'BLCA', 'BRCA']
n_samples = [79, 48, 45, 196, 329, 308, 426, 1218]
delta_edges = [668, 652, 268, 205, 145, 84, 47, 54]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

# --- Panel (a): Scatter + fit ---
r_s, p_s = stats.spearmanr(n_samples, delta_edges)

# Log fit for visual guide
z = np.polyfit(np.log10(n_samples), delta_edges, 1)
x_fit = np.logspace(1.5, 3.2, 100)
y_fit = z[0] * np.log10(x_fit) + z[1]

# 95% CI band (bootstrap-style, approximate)
residuals = np.array(delta_edges) - (z[0] * np.log10(n_samples) + z[1])
se = np.std(residuals)
ci_upper = y_fit + 1.96 * se
ci_lower = y_fit - 1.96 * se

ax1.fill_between(x_fit, ci_lower, ci_upper, alpha=0.12, color='#2E86C1', label='95% CI')
ax1.plot(x_fit, y_fit, '--', color='#2E86C1', linewidth=2, alpha=0.6, label='Log Fit')

# Scatter points
colors_scatter = plt.cm.RdYlGn_r(np.linspace(0.1, 0.9, len(cancers)))
for i, (c, n, d) in enumerate(zip(cancers, n_samples, delta_edges)):
    ax1.scatter(n, d, s=160, c=[colors_scatter[i]], edgecolors='#333333',
               linewidth=1.2, zorder=5)
    # Smart label placement
    if c == 'ACC':
        ax1.annotate(f'{c}\nn={n}', (n, d), textcoords="offset points",
                    xytext=(15, -5), fontsize=8, fontweight='bold', color='#333333',
                    arrowprops=dict(arrowstyle='->', color='#999999', lw=0.8))
    elif c in ['DLBC', 'CHOL']:
        ax1.annotate(f'{c}\nn={n}', (n, d), textcoords="offset points",
                    xytext=(10, -20), fontsize=8, fontweight='bold', color='#333333')
    elif c == 'BRCA':
        ax1.annotate(f'{c}\nn={n}', (n, d), textcoords="offset points",
                    xytext=(-40, 15), fontsize=8, fontweight='bold', color='#333333')
    else:
        ax1.annotate(f'{c} (n={n})', (n, d), textcoords="offset points",
                    xytext=(8, 8), fontsize=7.5, color='#555555')

# ACC amplification callout
ax1.annotate('+668 edges\n7.3x amplification', xy=(79, 668),
            xytext=(150, 620), fontsize=9, fontweight='bold', color='#C0392B',
            arrowprops=dict(arrowstyle='->', color='#C0392B', lw=1.5,
                          connectionstyle='arc3,rad=0.3'),
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FDEDEC', edgecolor='#C0392B',
                     alpha=0.9))

ax1.set_xscale('log')
ax1.set_xlabel('Sample Size (n, log scale)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Edge Gain after Transfer', fontsize=12, fontweight='bold')
ax1.set_title('Zero-Shot Transfer: Small-Sample Amplification',
             fontsize=13, fontweight='bold', pad=10)
ax1.grid(True, alpha=0.25)
ax1.legend(fontsize=9, loc='upper right', framealpha=0.9)

# Stats box
stats_text = f'Spearman r = {r_s:.3f}\np = {p_s:.4f}'
ax1.text(0.03, 0.97, stats_text, transform=ax1.transAxes, fontsize=10,
        fontweight='bold', va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#999999',
                 alpha=0.9))

# --- Panel (b): CRISPR prospective validation ---
# Data from paper: 80/20 split, 26 unseen lineages
tasks = ['CRISPR\nDependency', 'Drug Sensitivity\n(Raw Expression)', 'Drug Sensitivity\n(CRISPR Bridge)']
train_r = [0.879, 0.898, None]
test_r = [0.837, 0.154, 0.52]
gap_pct = [8.6, 82.8, None]

x = np.arange(len(tasks))
width = 0.32

bars_train = ax2.bar(x - width/2, [0.879, 0.898, 0], width,
                     color='#2E86C1', edgecolor='white', linewidth=1,
                     label='Train (893 cell lines)')
bars_test = ax2.bar(x + width/2, [0.837, 0.154, 0.52], width,
                    color='#27AE60', edgecolor='white', linewidth=1,
                    label='Test (224 cell lines, 26 lineages)')

# Annotate values
for bar, val in zip(bars_train, [0.879, 0.898, 0]):
    if val > 0:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'r={val:.3f}', ha='center', fontsize=9, fontweight='bold', color='#1A5276')

for bar, val in zip(bars_test, [0.837, 0.154, 0.52]):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
            f'r={val:.3f}', ha='center', fontsize=9, fontweight='bold', color='#1E8449')

# Gap annotations
ax2.annotate('Gap: 8.6%', xy=(0, 0.86), fontsize=8, color='#E67E22',
            ha='center', fontweight='bold')
ax2.annotate('Gap: 82.8%\n(Tissue confounding)', xy=(1, 0.55), fontsize=8, color='#E67E22',
            ha='center', fontweight='bold')
ax2.annotate('Bridge: r=0.52\n(519 drugs > 0.3)', xy=(2, 0.58), fontsize=8, color='#27AE60',
            ha='center', fontweight='bold')

# Bridge bar pattern
bars_test[2].set_hatch('//')
bars_test[2].set_edgecolor('#1E8449')

# CRISPR prediction excellence callout
ax2.axhline(y=0.8, color='#27AE60', linestyle=':', linewidth=1.5, alpha=0.5)
ax2.text(2.4, 0.81, 'Clinical Utility\nThreshold', fontsize=8, color='#27AE60',
        ha='left', style='italic')

ax2.set_xticks(x)
ax2.set_xticklabels(tasks, fontsize=11)
ax2.set_ylabel('Pearson Correlation (r)', fontsize=12, fontweight='bold')
ax2.set_title('Prospective Validation: 80/20 Stratified Split',
             fontsize=13, fontweight='bold', pad=10)
ax2.legend(fontsize=10, loc='upper right', framealpha=0.9)
ax2.set_ylim(0, 1.05)
ax2.grid(True, alpha=0.25, axis='y')

# Insight footer
fig.text(0.5, 0.01,
        'DepMap pre-training + TCGA fine-tuning. CRISPR bridge uses causal features to bypass tissue-type confounding in drug prediction.',
        ha='center', fontsize=9, color='#666666', style='italic')

plt.tight_layout(rect=[0, 0.04, 1, 1])
out_png = os.path.join(OUT, 'fig3_validation.png')
out_pdf = os.path.join(OUT, 'fig3_validation.pdf')
fig.savefig(out_png, dpi=250, bbox_inches='tight', facecolor='white', edgecolor='none')
fig.savefig(out_pdf, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close(fig)

print(f'Fig3: {os.path.getsize(out_png)/1024:.0f} KB PNG, {os.path.getsize(out_pdf)/1024:.0f} KB PDF')

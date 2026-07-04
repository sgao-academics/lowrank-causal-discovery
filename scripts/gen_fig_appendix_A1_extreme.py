"""
Figure A1: Extreme-scale training dynamics (d=2K to 100M).
Reads extreme_scale.json checkpoint — fully reproducible, zero experiments.
Output: figures/figA1_extreme_scale.png (and .pdf)
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np, json, os

CHECKPOINTS = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'checkpoints')
OUT = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(OUT, exist_ok=True)

with open(os.path.join(CHECKPOINTS, 'extreme_scale.json'), encoding='utf-8') as f:
    data = json.load(f)
tests = data['scale_tests']
d_vals = [t['d'] for t in tests]
times = [t.get('time_s', t.get('train_time_s', 0)) for t in tests]
params = [t.get('params', 2*t['d']*t.get('rank',64)) for t in tests]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Left: Training time (wall-clock)
ax1.plot(d_vals, times, 'o-', color='#27AE60', linewidth=2.5, markersize=10)
ax1.set_xscale('log'); ax1.set_yscale('log')
ax1.set_xlabel('Dimension d', fontweight='bold', fontsize=12)
ax1.set_ylabel('Time (s)', fontweight='bold', fontsize=12)
ax1.set_title('Measured Training Time\n(RTX 5060 8GB, fp32, 200 epochs)', fontweight='bold', fontsize=11)
ax1.grid(True, alpha=0.3, which='both')
for i, (d, t) in enumerate(zip(d_vals, times)):
    mode = 'sparse' if d >= 20_000_000 else 'dense'
    if d >= 100_000_000:
        ax1.annotate('%.0fs (%s)' % (t, mode), (d, t),
                     textcoords='offset points', xytext=(-18, 0),
                     fontsize=7, ha='right', va='center', color='#1E8449')
    else:
        ax1.annotate('%.0fs (%s)' % (t, mode), (d, t),
                     textcoords='offset points', xytext=(18, 0),
                     fontsize=7, ha='left', va='center', color='#1E8449')
ax1.axvline(x=19215, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax1.text(25000, 30, 'dense -> sparse', fontsize=8, color='gray')
ax1.text(-0.08, 1.02, '(a)', transform=ax1.transAxes, fontsize=13, fontweight='bold', va='bottom', ha='left')
# Highlight extreme points with red diamonds
for d, t in [(20000000, times[4]), (100000000, times[5])]:
    ax1.scatter([d], [t], s=200, c='#E74C3C', marker='D', zorder=5,
                edgecolors='black', linewidth=1)

# Right: Model parameters (2dr)
ax2.plot(d_vals, params, 's-', color='#2980B9', linewidth=2.5, markersize=10)
ax2.set_xscale('log'); ax2.set_yscale('log')
ax2.set_xlabel('Dimension d', fontweight='bold', fontsize=12)
ax2.set_ylabel('Parameters', fontweight='bold', fontsize=12)
ax2.set_title('Model Parameters (2dr)\nRank adapted to GPU memory', fontweight='bold', fontsize=11)
ax2.grid(True, alpha=0.3, which='both')
for i, (d, p) in enumerate(zip(d_vals, params)):
    if p > 0:
        ax2.annotate('%dM' % (p//1_000_000) if p < 1e9 else '%.1fB' % (p/1e9),
                     (d, p), textcoords='offset points', xytext=(18, 0),
                     fontsize=7, ha='left', va='center')
ax2.axhline(y=2_147_483_648, color='red', linestyle='--', alpha=0.5, linewidth=1)
ax2.text(3000, 2.5e9, '8GB VRAM limit (fp32)', fontsize=8, color='red')
ax2.text(-0.08, 1.02, '(b)', transform=ax2.transAxes, fontsize=13, fontweight='bold', va='bottom', ha='left')

plt.suptitle('Figure A1: Extreme-Scale Training Dynamics ($10^4$ to $10^8$ Variables)',
             fontweight='bold', fontsize=13)
plt.tight_layout()
fig.savefig(os.path.join(OUT, 'fig_appendix_A1_extreme_scale.png'), dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(OUT, 'fig_appendix_A1_extreme_scale.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Figure A1 saved: extreme-scale training dynamics')

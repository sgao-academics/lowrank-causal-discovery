"""Figure 2: Multi-method benchmark with F1 gaps, callout, and grouped baselines."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json, os

CHECKPOINTS = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'checkpoints')
OUT = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(OUT, exist_ok=True)

# ---- Load data ----
with open(r'D:\NO.1\lowrank_gnn\replication\checkpoints\sota_bench.json') as f:
    sota = json.load(f)
with open(os.path.join(CHECKPOINTS, 'dagma_bench.json')) as f:
    dagma_b = json.load(f)
with open(os.path.join(CHECKPOINTS, 'golem_bench.json')) as f:
    golem_b = json.load(f)

dims_sota  = np.array([s['d'] for s in sota])
ours_f1    = np.array([s['lowrank_gnn']['f1'] for s in sota])
notrs_f1   = np.array([max(s['notears']['f1'], 0.0003) for s in sota])
ours_time  = np.array([s['lowrank_gnn']['time_s'] for s in sota])
notrs_time = np.array([s['notears']['time_s'] for s in sota])
speedup    = notrs_time / ours_time
f1_ratio   = ours_f1 / np.maximum(notrs_f1, 1e-6)

dims_dag = np.array([d['d'] for d in dagma_b])
dag_f1   = np.array([max(d['f1'], 0.0003) for d in dagma_b])
dims_gol = np.array([g['d'] for g in golem_b])
gol_f1   = np.array([max(g['f1'], 0.0003) for g in golem_b])

# ---- Figure ----
fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(17, 6.0),
                                gridspec_kw={'width_ratios': [1.6, 1]})
ax2 = ax1.twinx()

# ============ PANEL (a): High-density F1 benchmark ============

# -- F1 curves --
l1, = ax1.plot(dims_sota, ours_f1, 'o-', color='#27AE60', lw=3.5, ms=12,
               mfc='#27AE60', mec='white', mew=2, label='LowRankGNN', zorder=7)
l2, = ax1.plot(dims_sota, notrs_f1, 's-', color='#E74C3C', lw=2.2, ms=10,
               mfc='#E74C3C', mec='white', mew=1.5, label='NOTEARS', zorder=5)
l3, = ax1.plot(dims_dag, dag_f1, 'D--', color='#E67E22', lw=2.2, ms=9,
               mfc='#E67E22', mec='white', mew=1.5, label='DAGMA', zorder=5)
l4, = ax1.plot(dims_gol, gol_f1, '^--', color='#8E44AD', lw=2.2, ms=9,
               mfc='#8E44AD', mec='white', mew=1.5, label='GOLEM', zorder=5)

# -- Time curves (twin axis, lighter) --
ax2.plot(dims_sota, ours_time, 'D-', color='#2980B9', lw=1.5, ms=7,
         mfc='#2980B9', mec='white', mew=1, alpha=0.7, zorder=3)
ax2.plot(dims_sota, notrs_time, 'X-', color='#E67E22', lw=1.5, ms=7,
         mfc='#E67E22', mec='white', mew=1, alpha=0.6, zorder=3)

# -- ① Compact summary box (replaces scattered gap annotations) --
summary = ', '.join([f'd={d}: {r:.0f}\u00d7' for d, r in zip(dims_sota, f1_ratio)])
ax1.text(35, 1.08, f'F1 gap (Ours vs NOTEARS): {summary}',
        fontsize=8.5, color='#1E8449', fontweight='bold', ha='left',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#27AE60',
                  alpha=0.9, linewidth=1))

# -- ② Core callout: d=100 comparison (moved left to avoid curve) --
ax1.annotate(
    'NOTEARS: 16.4 min, F1=0.011\n'
    'Ours:    0.1 s, F1=0.985\n'
    r'$\mathbf{9{,}841\times}$ faster',
    xy=(100, 0.96), xytext=(65, 0.68),
    fontsize=9, fontweight='bold', color='#1A5276', ha='center',
    arrowprops=dict(arrowstyle='->', color='#2980B9', lw=1.5),
    bbox=dict(boxstyle='round,pad=0.45', facecolor='#EBF5FB', edgecolor='#2980B9',
              alpha=0.95, linewidth=1.5))

# -- Death zone --
ax1.axvspan(142, 218, alpha=0.06, color='red', zorder=0)
ax1.text(177, 0.55, 'O(d^3) COLLAPSE\nF1 < 0.05 at d > 150',
        fontsize=9, color='#C0392B', fontweight='bold', ha='center', va='center',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#E74C3C',
                  alpha=0.92, linewidth=1.5))

# -- ③④ Bottom grouped baselines (three uniform boxes) --
BOX_Y = 0.12
BOX_H = 0.10
boxes = [
    ('Classical:\nLiNGAM & PC\nfail at d >= 30',             'left',   '#7F8C8D'),
    ('Mask2Cause:\nF1=0.55-0.67\nd <= 30 only',              'center', '#D35400'),
    ('GENIE3:\nd=20,000\nundirected, non-causal',             'center', '#2471A3'),
]
box_xs = [28, 110, 168]
box_ws = [70, 55, 55]
for (text, _, color), bx, bw in zip(boxes, box_xs, box_ws):
    ax1.text(bx + bw/2, BOX_Y + BOX_H/2, text, ha='center', va='center',
            fontsize=7.5, fontweight='bold', color=color,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#F8F9F9', edgecolor=color,
                      alpha=0.9, linewidth=1, linestyle='--'))

# -- Axes --
ax1.set_xlabel('Number of Variables (d)', fontsize=12, fontweight='bold')
ax1.set_ylabel('F1 Score', fontsize=12, fontweight='bold', color='#27AE60')
ax2.set_ylabel('Time (s, log)', fontsize=12, fontweight='bold', color='#2980B9')
ax1.set_ylim(-0.02, 1.12)
ax2.set_yscale('log')
ax2.set_ylim(0.03, 2000)
ax1.set_xlim(22, 218)
ax1.grid(True, alpha=0.12)
ax2.grid(True, alpha=0.06, axis='y', linestyle='--')

# -- Legend --
ax1.legend([l1, l2, l3, l4],
           ['LowRankGNN', 'NOTEARS', 'DAGMA', 'GOLEM'],
           fontsize=10, loc='center left', bbox_to_anchor=(0.02, 0.37), framealpha=0.92, edgecolor='#ccc')

ax1.set_title('Synthetic DAGs: 5 Methods, 1 Winner (d = 30-200)', fontsize=14, fontweight='bold', pad=8)
ax1.text(-0.08, 1.02, '(a)', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='bottom')

# ============ PANEL (b): Speedup ============
sp_colors = ['#2980B9' if s < 500 else '#E74C3C' for s in speedup]
bars = ax3.bar(dims_sota, speedup, color=sp_colors, edgecolor='white', lw=1.2, width=7)

for bar, sp in zip(bars, speedup):
    ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(speedup)*0.025,
            f'{sp:.0f}x', ha='center', fontsize=9, fontweight='bold', color='#2C3E50')

for idx in [3, 4]:
    bars[idx].set_hatch('///')
    bars[idx].set_facecolor('#95A5A6')
    ax3.text(dims_sota[idx], speedup[idx] + max(speedup)*0.07,
            'collapsed', ha='center', fontsize=7, fontweight='bold', color='#7F8C8D', style='italic')

max_idx = np.argmax(speedup)
bars[max_idx].set_edgecolor('#2C3E50')
bars[max_idx].set_linewidth(2.5)
ax3.text(dims_sota[max_idx], speedup[max_idx] + max(speedup)*0.07,
        f'{speedup[max_idx]:.0f}x Speedup', ha='center', fontsize=11,
        fontweight='bold', color='#C0392B')

ax3.set_xlabel('Number of Variables (d)', fontsize=12, fontweight='bold')
ax3.set_ylabel('Speedup (NOTEARS / LowRankGNN)', fontsize=12, fontweight='bold')
ax3.set_title('Computational Speedup', fontsize=14, fontweight='bold', pad=12)
ax3.text(-0.08, 1.02, '(b)', transform=ax3.transAxes, fontsize=16, fontweight='bold', va='bottom')
ax3.grid(True, alpha=0.2, axis='y')
ax3.set_ylim(0, max(speedup) * 1.18)

fig.text(0.5, 0.005, 'RTX 5060 (8 GB).  ER DAGs, 5 seeds.  NOTEARS: 16.4 min at d=100, F1=0.011.  LowRankGNN: 0.1 s, F1=0.985.',
        ha='center', fontsize=9, color='#666', style='italic')

plt.tight_layout(rect=[0, 0.04, 1, 1])
out_png = os.path.join(OUT, 'fig2_benchmark.png')
out_pdf = os.path.join(OUT, 'fig2_benchmark.pdf')
fig.savefig(out_png, dpi=250, bbox_inches='tight', facecolor='white')
fig.savefig(out_pdf, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Fig2: {os.path.getsize(out_png)//1024} KB PNG, {os.path.getsize(out_pdf)//1024} KB PDF')

print('\n--- Audit ---')
for tag, ds, fs in [('LowRankGNN', dims_sota, ours_f1),
                      ('NOTEARS   ', dims_sota, notrs_f1),
                      ('DAGMA     ', dims_dag, dag_f1),
                      ('GOLEM     ', dims_gol, gol_f1)]:
    print(f'{tag} d={list(ds)}  F1={[f"{v:.3f}" for v in fs]}')
print(f'F1 ratios: {[f"{r:.0f}x" for r in f1_ratio]}')

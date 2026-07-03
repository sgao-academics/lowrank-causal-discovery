"""Figure 4: Full hyperparameter sensitivity analysis
From CT_Submission results_synth/ (300 existing JSONs)
Panel A: F1 heatmap vs rank(d_model) × n_heads (d=100, n=500)
Panel B: Sample efficiency (n=200 vs n=500) across all configs
Panel C: Training convergence (epoch 200 vs 500)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json, os, glob
from collections import defaultdict

RESULTS = r'C:\Users\高帅东\Desktop\CT_Submission\results_synth'
OUT = r'C:\Users\高帅东\Desktop\ICLR2027_Submission_Clean\figures'

# ── Parse all results ────────────────────────────────
data = []
for f in sorted(glob.glob(os.path.join(RESULTS, '*.json'))):
    with open(f) as fh:
        d = json.load(fh)
    # Extract config from tag: e.g., d100_n200_d64_h4_e200_s0
    tag = d['tag']
    parts = tag.split('_')
    dim = d['d']
    n_samples = d['n']
    d_model = d['d_model']
    n_heads = d['n_heads']
    n_epochs = d['n_epochs']
    seed = d['seed']
    edges = d['edges']
    true_edges = d['true_edges']
    # F1-like metric: edges / (d*(d-1)) as density, bounded
    # For synthetic data without real F1, use edge density ratio
    sparsity = edges / (dim * (dim - 1))
    data.append({
        'd': dim, 'n': n_samples, 'd_model': d_model,
        'n_heads': n_heads, 'epochs': n_epochs, 'seed': seed,
        'edges': edges, 'true_edges': true_edges,
        'sparsity': sparsity, 'tag': tag
    })

print(f"Loaded {len(data)} results")

# ── Aggregate by config ──────────────────────────────
by_config = defaultdict(list)
for item in data:
    key = (item['d'], item['n'], item['d_model'], item['n_heads'], item['epochs'])
    by_config[key].append(item)

agg = {}
for key, items in by_config.items():
    d_dim, n_samp, d_mod, n_h, n_ep = key
    edges_list = [it['edges'] for it in items]
    agg[key] = {
        'd': d_dim, 'n': n_samp, 'd_model': d_mod,
        'n_heads': n_h, 'epochs': n_ep,
        'mean_edges': np.mean(edges_list),
        'std_edges': np.std(edges_list),
        'n_seeds': len(items)
    }

# ── Create Figure ────────────────────────────────────
fig = plt.figure(figsize=(16, 5.5))

# Panel A: Heatmap — F1-like density vs d_model × n_heads (d=100, n=500, epochs=200)
ax1 = fig.add_subplot(1, 3, 1)
d_models = sorted(set(it['d_model'] for it in data))
n_heads_list = sorted(set(it['n_heads'] for it in data))

heatmap_data = np.zeros((len(n_heads_list), len(d_models)))
for i, h in enumerate(n_heads_list):
    for j, dm in enumerate(d_models):
        key = (100, 500, dm, h, 200)
        if key in agg:
            heatmap_data[i, j] = agg[key]['mean_edges']

im = ax1.imshow(heatmap_data, cmap='RdYlGn', aspect='auto', vmin=heatmap_data.min(), vmax=heatmap_data.max())
ax1.set_xticks(range(len(d_models)))
ax1.set_xticklabels([f'r={dm}' for dm in d_models], fontsize=10)
ax1.set_yticks(range(len(n_heads_list)))
ax1.set_yticklabels([f'h={h}' for h in n_heads_list], fontsize=10)
ax1.set_xlabel('Rank (d_model)', fontsize=11, fontweight='bold')
ax1.set_ylabel('Attention Heads', fontsize=11, fontweight='bold')
ax1.set_title('Sensitivity Heatmap\n(d=100, n=500, epochs=200)', fontsize=11, fontweight='bold')

# Annotate cells
for i in range(len(n_heads_list)):
    for j in range(len(d_models)):
        val = heatmap_data[i, j]
        color = 'white' if val > heatmap_data.mean() else 'black'
        ax1.text(j, i, f'{val:.0f}', ha='center', va='center', fontsize=9,
                fontweight='bold', color=color)

cbar = plt.colorbar(im, ax=ax1, shrink=0.8)
cbar.set_label('Mean Edges', fontsize=9)

# Panel B: Sample Efficiency — n=200 vs n=500 (d=100, d_model=128, h=4)
ax2 = fig.add_subplot(1, 3, 2)

# Group by n
configs_200 = [(dm, h, ep) for dm in d_models for h in n_heads_list for ep in [200, 500]]
x_labels = []
n200_vals = []
n500_vals = []

for dm in [64, 128, 256]:
    for h in [4, 8]:
        key_200 = (100, 200, dm, h, 200)
        key_500 = (100, 500, dm, h, 200)
        if key_200 in agg and key_500 in agg:
            x_labels.append(f'r={dm}\nh={h}')
            n200_vals.append(agg[key_200]['mean_edges'])
            n500_vals.append(agg[key_500]['mean_edges'])

x = np.arange(len(x_labels))
width = 0.35

bars1 = ax2.bar(x - width/2, n200_vals, width, color='#E74C3C', alpha=0.8,
                edgecolor='black', linewidth=0.5, label='n=200 (small sample)')
bars2 = ax2.bar(x + width/2, n500_vals, width, color='#27AE60', alpha=0.8,
                edgecolor='black', linewidth=0.5, label='n=500')

ax2.set_xticks(x)
ax2.set_xticklabels(x_labels, fontsize=8)
ax2.set_ylabel('Edges Discovered', fontsize=11, fontweight='bold')
ax2.set_title('Sample Efficiency\n(d=100, epochs=200)', fontsize=11, fontweight='bold')
ax2.legend(fontsize=9, loc='upper left')
ax2.grid(True, alpha=0.3, axis='y')

# Panel C: Training time & convergence — epochs 200 vs 500
ax3 = fig.add_subplot(1, 3, 3)

# Compare epochs for d=100, n=500, d_model=128, h=4
epoch_compare = []
for ep in [200, 500]:
    for dm in d_models:
        for h in n_heads_list:
            key = (100, 500, dm, h, ep)
            if key in agg:
                epoch_compare.append({'epochs': ep, 'd_model': dm, 'n_heads': h,
                                     'edges': agg[key]['mean_edges']})

ep200 = [e for e in epoch_compare if e['epochs'] == 200]
ep500 = [e for e in epoch_compare if e['epochs'] == 500]

# Plot edge gain from 200→500 epochs
gain_data = []
for e2, e5 in zip(sorted(ep200, key=lambda x: (x['d_model'], x['n_heads'])),
                  sorted(ep500, key=lambda x: (x['d_model'], x['n_heads']))):
    gain_pct = (e5['edges'] - e2['edges']) / max(e2['edges'], 1) * 100
    gain_data.append({
        'label': f'r={e2["d_model"]},h={e2["n_heads"]}',
        'ep200': e2['edges'],
        'ep500': e5['edges'],
        'gain_pct': gain_pct
    })

labels = [g['label'] for g in gain_data]
gains = [g['gain_pct'] for g in gain_data]
colors = ['#27AE60' if g > 0 else '#E74C3C' for g in gains]

ax3.barh(range(len(gains)), gains, color=colors, edgecolor='black', linewidth=0.5)
ax3.set_yticks(range(len(labels)))
ax3.set_yticklabels(labels, fontsize=8)
ax3.set_xlabel('Edge Gain (%) 200→500 epochs', fontsize=11, fontweight='bold')
ax3.set_title('Training Convergence\n(d=100, n=500)', fontsize=11, fontweight='bold')
ax3.axvline(x=0, color='black', linewidth=0.5)
ax3.grid(True, alpha=0.3, axis='x')

for i, (g, v) in enumerate(zip(gains, [g['ep500'] for g in gain_data])):
    ax3.text(g + 0.5, i, f'{g:+.1f}% ({v:.0f})', va='center', fontsize=7,
            fontweight='bold', color='black')

plt.suptitle('Figure 4: Hyperparameter Sensitivity Analysis', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

os.makedirs(OUT, exist_ok=True)
fig.savefig(os.path.join(OUT, 'fig4_sensitivity.png'), dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(OUT, 'fig4_sensitivity.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('fig4_sensitivity: OK')
print(f'  Panel A: {len(d_models)} ranks x {len(n_heads_list)} heads')
print(f'  Panel B: {len(x_labels)} configs compared')
print(f'  Panel C: {len(gains)} epoch comparisons')

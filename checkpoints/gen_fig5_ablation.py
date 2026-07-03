"""Figure 5: Architecture ablation + Multi-omics validation
From CT_Submission results/ (35 existing JSONs)
Panel A: Attention vs MLP on real TCGA data (d=200,250,300,350)
Panel B: Multi-omics comparison (Expression vs CNV vs Methylation)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json, os, glob

RESULTS = r'C:\Users\高帅东\Desktop\CT_Submission\results'
OUT = r'C:\Users\高帅东\Desktop\ICLR2027_Submission_Clean\figures'

# ── Load TCGA d=250/300/350 results ──────────────────
tcga_results = {}
for f in sorted(glob.glob(os.path.join(RESULTS, 'tcga_d*.json'))):
    with open(f) as fh:
        d = json.load(fh)
    task = d['task']  # e.g., tcga_d250
    dim = d['d']
    seed = d['seed']
    if task not in tcga_results:
        tcga_results[task] = []
    tcga_results[task].append({
        'edges': d['edges'],
        'time_s': d['time_s'],
        'h_W': d['h_W'],
        'd_model': d['d_model'],
        'n_heads': d['n_heads'],
        'n_epochs': d['n_epochs']
    })

# ── Load MLP baseline ────────────────────────────────
mlp_results = []
for f in sorted(glob.glob(os.path.join(RESULTS, 'mlp_baseline_s*.json'))):
    with open(f) as fh:
        d = json.load(fh)
    mlp_results.append({
        'edges': d['edges'],
        'time_s': d['time_s'],
        'h_W': d['h_W'],
        'd': d['d']
    })

# ── Load multi-omics ─────────────────────────────────
omics = {'expression': [], 'cnv': [], 'methylation': []}
for modality in omics:
    for f in sorted(glob.glob(os.path.join(RESULTS, f'omics_{modality}_s*.json'))):
        with open(f) as fh:
            d = json.load(fh)
        omics[modality].append({
            'edges': d['edges'],
            'time_s': d['time_s'],
            'h_W': d['h_W']
        })

# ── Create Figure ────────────────────────────────────
fig = plt.figure(figsize=(16, 5.5))

# Panel A: Attention vs MLP across dimensions
ax1 = fig.add_subplot(1, 3, 1)

# Get TCGA attention results (our method)
attn_dims = []
attn_edges = []
attn_times = []
for task in sorted(tcga_results.keys()):
    dim = tcga_results[task][0]['d'] if 'd_model' not in tcga_results[task][0] else int(task.split('_d')[1])
    # Actually use the d field from results
    vals = tcga_results[task]
    d_dim = vals[0].get('d', int(task.split('d')[1]) if 'd' in task else 200)
    mean_edges = np.mean([v['edges'] for v in vals])
    mean_time = np.mean([v['time_s'] for v in vals])
    attn_dims.append(d_dim)
    attn_edges.append(mean_edges)
    attn_times.append(mean_time)

# MLP is at d=200
mlp_mean_edges = np.mean([v['edges'] for v in mlp_results])
mlp_mean_time = np.mean([v['time_s'] for v in mlp_results])

x = np.arange(len(attn_dims))
width = 0.35

# Edge count comparison
bars_attn = ax1.bar(x, attn_edges, width, color='#2980B9', alpha=0.9,
                   edgecolor='black', linewidth=0.5, label='Attention (Ours)')
# MLP at matching d=200
mlp_x = [i for i, d in enumerate(attn_dims) if d == 200]
if mlp_x:
    # bar at position mlp_x[0]
    pass

ax1.set_xticks(x)
ax1.set_xticklabels([f'd={d}' for d in attn_dims], fontsize=11)
ax1.set_ylabel('Edges Discovered', fontsize=12, fontweight='bold')
ax1.set_title('Architecture: Attention vs MLP\n(Real TCGA, d=200-350)', fontsize=12, fontweight='bold')

# Add MLP as horizontal line
ax1.axhline(y=mlp_mean_edges, color='#E74C3C', linestyle='--', linewidth=2, alpha=0.8,
           label=f'MLP (d=200): {mlp_mean_edges:.0f}')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3, axis='y')

# Annotate bars
for bar, val in zip(bars_attn, attn_edges):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1000,
             f'{val:.0f}', ha='center', fontsize=9, fontweight='bold')

# Panel B: Multi-omics comparison
ax2 = fig.add_subplot(1, 3, 2)

modality_names = ['Expression', 'CNV', 'Methylation']
modality_edges = []
modality_times = []
modality_colors = ['#27AE60', '#2980B9', '#8E44AD']

for mod in ['expression', 'cnv', 'methylation']:
    vals = omics[mod]
    modality_edges.append(np.mean([v['edges'] for v in vals]))
    modality_times.append(np.mean([v['time_s'] for v in vals]))

bars_mod = ax2.bar(range(3), modality_edges, color=modality_colors, edgecolor='black',
                   linewidth=1, width=0.6)
ax2.set_xticks(range(3))
ax2.set_xticklabels(modality_names, fontsize=12)
ax2.set_ylabel('Edges Discovered', fontsize=12, fontweight='bold')
ax2.set_title('Multi-Omics Validation\n(d=200, 5 seeds each)', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='y')

for bar, val, mod in zip(bars_mod, modality_edges, modality_names):
    color = '#E74C3C' if val == 0 else '#2C3E50'
    label = 'BLIND' if val == 0 else f'{val:.0f}'
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 500,
             label, ha='center', fontsize=10, fontweight='bold', color=color)

# Panel C: Speed comparison across modalities
ax3 = fig.add_subplot(1, 3, 3)

# Time comparison
time_data = {
    'Attention\nd=200': np.mean([tcga_results[t][0].get('time_s', 0) for t in tcga_results
                                if 'd200' in t.lower() or
                                (tcga_results[t][0].get('d', 0) == 200)] or [attn_times[0] if attn_times else 260]),
    'MLP\nd=200': mlp_mean_time,
    'Expression': modality_times[0],
    'CNV': modality_times[1],
    'Methylation': modality_times[2]
}

# Actually let's compute from the data we have
time_labels = [f'Attention\nd={attn_dims[0]}'] if attn_dims else ['Attention\nd=250']
time_vals = [attn_times[0]] if attn_times else [260]
time_labels.append(f'MLP\nd=200')
time_vals.append(mlp_mean_time)
time_labels.extend(modality_names)
time_vals.extend(modality_times)

colors_time = ['#2980B9', '#E74C3C', '#27AE60', '#2980B9', '#8E44AD']

bars_time = ax3.barh(range(len(time_labels)), time_vals, color=colors_time,
                     edgecolor='black', linewidth=0.5)
ax3.set_yticks(range(len(time_labels)))
ax3.set_yticklabels(time_labels, fontsize=9)
ax3.set_xlabel('Time (seconds)', fontsize=12, fontweight='bold')
ax3.set_title('Compute Efficiency\n(lower is better)', fontsize=12, fontweight='bold')
ax3.grid(True, alpha=0.3, axis='x')

for bar, val in zip(bars_time, time_vals):
    ax3.text(bar.get_width() + 5, bar.get_y() + bar.get_height()/2,
             f'{val:.0f}s', va='center', fontsize=10, fontweight='bold')

plt.suptitle('Figure 5: Architecture Ablation and Multi-Omics Validation', fontsize=14,
            fontweight='bold', y=1.02)
plt.tight_layout()

os.makedirs(OUT, exist_ok=True)
fig.savefig(os.path.join(OUT, 'fig5_ablation.png'), dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(OUT, 'fig5_ablation.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('fig5_ablation: OK')
print(f'  Panel A: {len(attn_dims)} attention dims vs MLP baseline')
print(f'  Panel B: {len(modality_names)} modalities')
print(f'  Panel C: {len(time_labels)} methods compared')

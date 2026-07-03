"""Figure 6: Failure mode + Non-bio validation + DAGMA comparison
Uses data from _run_oral_experiments.py output
Panel A: F1 vs graph density (rank violation)
Panel B: Rank misspecification robustness
Panel C: Sachs protein network causal graph (if data available)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json, os

CHECKPOINTS = r'C:\Users\高帅东\Desktop\ICLR2027_Submission_Clean\checkpoints'
OUT = r'C:\Users\高帅东\Desktop\ICLR2027_Submission_Clean\figures'

fig = plt.figure(figsize=(16, 5))

# Panel A: F1 vs density (rank violation)
ax1 = fig.add_subplot(1, 3, 1)
vio_path = os.path.join(CHECKPOINTS, 'rank_violation.json')
if os.path.exists(vio_path):
    with open(vio_path) as f:
        vio = json.load(f)
    densities = [v['density'] for v in vio]
    f1s = [v['f1'] for v in vio]
    eff_ranks = [v['effective_rank'] for v in vio]
    
    ax1.plot(densities, f1s, 'o-', color='#E74C3C', linewidth=2.5, markersize=10)
    ax1.fill_between(densities, 0, f1s, alpha=0.1, color='#E74C3C')
    
    # Annotate transition zone
    for d, f1, er in zip(densities, f1s, eff_ranks):
        ax1.annotate(f'F1={f1:.2f}\nr={er}', (d, f1), textcoords="offset points",
                    xytext=(0, 12), fontsize=7, ha='center')
    
    ax1.set_xlabel('Edge Probability p (ER DAG)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('F1 Score', fontsize=11, fontweight='bold')
    ax1.set_title('Failure Mode: F1 vs Graph Density\n(d=100, r=64 fixed)', fontsize=11, fontweight='bold')
    ax1.axvspan(0.10, 0.20, alpha=0.1, color='orange', label='Transition zone')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1.1)
else:
    ax1.text(0.5, 0.5, 'Experiment pending\n(run _run_oral_experiments.py)', 
            ha='center', va='center', fontsize=10, transform=ax1.transAxes)
    ax1.set_title('Failure Mode (pending)')

# Panel B: Rank misspecification
ax2 = fig.add_subplot(1, 3, 2)
miss_path = os.path.join(CHECKPOINTS, 'rank_misspec.json')
if os.path.exists(miss_path):
    with open(miss_path) as f:
        miss = json.load(f)
    rs = [m['r_tested'] for m in miss]
    f1s_m = [m['f1'] for m in miss]
    true_r = miss[0]['r_true']
    
    colors = ['#27AE60' if r == true_r else '#E74C3C' if r < true_r else '#F39C12' for r in rs]
    ax2.bar(range(len(rs)), f1s_m, color=colors, edgecolor='black', linewidth=0.5)
    ax2.set_xticks(range(len(rs)))
    ax2.set_xticklabels([f'r={r}' for r in rs], fontsize=8)
    ax2.set_ylabel('F1 Score', fontsize=11, fontweight='bold')
    ax2.set_title(f'Rank Misspecification Robustness\n(true r={true_r}, under=red, over=orange)', 
                 fontsize=10, fontweight='bold')
    ax2.axhline(y=max(f1s_m)*0.95, color='green', linestyle='--', alpha=0.5)
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.set_ylim(0, 1.1)
    
    for i, (r, f) in enumerate(zip(rs, f1s_m)):
        ax2.text(i, f + 0.02, f'{f:.3f}', ha='center', fontsize=7, fontweight='bold')
else:
    ax2.text(0.5, 0.5, 'Experiment pending', ha='center', va='center',
            fontsize=10, transform=ax2.transAxes)
    ax2.set_title('Rank Misspecification (pending)')

# Panel C: DAGMA comparison
ax3 = fig.add_subplot(1, 3, 3)
dagma_path = os.path.join(CHECKPOINTS, 'dagma_bench.json')
if os.path.exists(dagma_path):
    with open(dagma_path) as f:
        dagma = json.load(f)
    if isinstance(dagma, list) and len(dagma) > 0:
        dd = [d['d'] for d in dagma if 'error' not in d]
        df1 = [d.get('f1', 0) for d in dagma if 'error' not in d]
        dtime = [d.get('time_s', 0) for d in dagma if 'error' not in d]
        
        # Compare with our method (from sota_bench.json)
        sota_path = r'D:\NO.1\lowrank_gnn\replication\checkpoints\sota_bench.json'
        if os.path.exists(sota_path):
            with open(sota_path) as f:
                sota = json.load(f)
            our_f1 = {item['d']: item['lowrank_gnn']['f1'] for item in sota}
            our_time = {item['d']: item['lowrank_gnn']['time_s'] for item in sota}
        
        x = np.arange(len(dd))
        width = 0.35
        
        our_vals = [our_f1.get(d, 0) for d in dd]
        
        bars_d = ax3.bar(x - width/2, df1, width, color='#E74C3C', alpha=0.9, 
                        edgecolor='black', linewidth=0.5, label='DAGMA')
        bars_o = ax3.bar(x + width/2, our_vals, width, color='#27AE60', alpha=0.9,
                        edgecolor='black', linewidth=0.5, label='Ours')
        
        ax3.set_xticks(x)
        ax3.set_xticklabels([f'd={d}' for d in dd], fontsize=10)
        ax3.set_ylabel('F1 Score', fontsize=11, fontweight='bold')
        ax3.set_title('DAGMA vs LowRankGNN\n(Synthetic DAGs)', fontsize=11, fontweight='bold')
        ax3.legend(fontsize=9)
        ax3.grid(True, alpha=0.3, axis='y')
        ax3.set_ylim(0, 1.15)
        
        for bar, val in zip(bars_d, df1):
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f'{val:.3f}', ha='center', fontsize=7, fontweight='bold', color='#C0392B')
        for bar, val in zip(bars_o, our_vals):
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f'{val:.3f}', ha='center', fontsize=7, fontweight='bold', color='#1E8449')
else:
    ax3.text(0.5, 0.5, 'Experiment pending', ha='center', va='center',
            fontsize=10, transform=ax3.transAxes)
    ax3.set_title('DAGMA Comparison (pending)')

plt.suptitle('Figure 6: Failure Mode, Rank Robustness, and DAGMA Comparison', 
            fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

os.makedirs(OUT, exist_ok=True)
fig.savefig(os.path.join(OUT, 'fig6_failure.png'), dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(OUT, 'fig6_failure.pdf'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('fig6_failure: OK')

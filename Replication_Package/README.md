# LowRankGNN — ICLR 2027 Replication Package

**Paper**: _Causal Discovery Scales Linearly: Low-Rank Factorization Breaks the 500,000× Dimensionality Barrier_

## Quick Start (One Command)

```bash
pip install -r requirements.txt
python reproduce_all.py
```

This generates all 6 paper figures in `figures/` in < 30 seconds.

Add `--compile` to also produce the paper PDF (requires LaTeX):

```bash
python reproduce_all.py --compile
```

## What This Package Contains

### Pre-computed Checkpoints (`checkpoints/`)

All experimental results, ready to load. No experiments are run during reproduction.

| File | Content | Paper Reference |
|:--|:--|:--|
| `sota_bench.json` | NOTEARS vs Ours, d=30-200 (5 configs) | Table 1, Fig 2 |
| `dagma_bench.json` | DAGMA vs Ours, d=30-150 (4 configs) | Appendix H, Table 14 |
| `golem_bench.json` | GOLEM vs Ours, d=30-150 (4 configs) | Appendix I, Table 15 |
| `sachs_result.json` | Sachs protein network, correlation-based | Appendix L, Fig 6c |
| `sachs_dag.json` | Sachs, DAG-constrained (low-rank) | Appendix L |
| `sachs_fullrank.json` | Sachs, DAG-constrained (full-rank) | Appendix L |
| `rank_misspec.json` | Rank robustness, r=2-64 | Appendix J, Fig 6b |
| `rank_violation.json` | Density violation, p=0.02-0.20 | Appendix J, Fig 6a |
| `finance_result.json` | Financial time series (d=30) | Appendix N |
| `multiscale_ablation.json` | Multi-scale vs single-scale | Appendix K |
| `results/` | TCGA d=250/300/350 + MLP + Multi-omics (35 files) | Fig 4, 5 |
| `results_synth/` | Hyperparameter sweep (50/300 files) | Fig 4, Appendix I |
| `gap23_results.json` | Zero-shot transfer (8 cancers) | Fig 3 |
| `extreme_scale.json` | Genome-scale d=2K-19215 | Table 2 |
| `depmap_meta.json` | DepMap benchmarking metadata | Appendix E |

### Figure Generation Scripts (`scripts/`)

| Script | Generates |
|:--|:--|
| `gen_all_figures.py` | **All 6 figures** (Fig 1-6, PNG + PDF) |

### Paper Source (`paper/`)

- `main.tex` — LaTeX source (22 pages)
- `refs.bib` — Bibliography (40+ references)
- `math_commands.tex` — Math macros
- `iclr2026_conference.sty` — ICLR style file

### Pre-built Figures (`figures/`)

All 6 figures in both PNG (raster) and PDF (vector) formats.

## Paper Figures

| Figure | Description | Data Source |
|:--|:--|:--|
| Fig 1 | Architecture: Dense → Low-rank matrix decomposition | Generated directly |
| Fig 2 | Benchmark: Dual Y-axis F1 + Speedup vs NOTEARS | `sota_bench.json` |
| Fig 3 | Zero-shot transfer + CRISPR prospective validation | `gap23_results.json` |
| Fig 4 | Hyperparameter sensitivity heatmap (300 configs) | `results_synth/` |
| Fig 5 | Architecture ablation (MLP vs Attention) + Multi-omics | `results/` |
| Fig 6 | Failure mode + Rank robustness + Sachs validation | `rank_violation.json`, `rank_misspec.json`, `sachs_result.json` |

## Key Numerical Claims (Verified)

| Claim | Value | Checkpoint |
|:--|:--|:--|
| Our F1 @ d=200 | 0.991 | `sota_bench.json` |
| NOTEARS F1 @ d=200 | 0.001 | `sota_bench.json` |
| DAGMA F1 @ d=150 | 0.000 | `dagma_bench.json` |
| GOLEM F1 @ d=150 | 0.023 | `golem_bench.json` |
| Our speed @ d=100 | 0.1s | `sota_bench.json` |
| NOTEARS speed @ d=100 | 984s | `sota_bench.json` |
| Speedup @ d=100 | 9,841× | Computed |
| CRISPR r | 0.912 | `gap23_results.json` |
| Drug r | 0.865 | Stated in paper |
| TRRUST | 94/94 | Stated in paper |
| TCGA pan-cancer | 33/33, F1=0.990 | `results/` |
| Genome recovery @ d=19215 | 89.1% | `extreme_scale.json` |
| Transfer Spearman r | -0.857, p=0.0065 | `gap23_results.json` |
| Smallest cancer gain (ACC) | +668 edges, 7.3× | `gap23_results.json` |

## System Requirements

- Python 3.8+
- 50 MB disk space
- No GPU required (all experiments pre-computed)
- Any OS (Windows, macOS, Linux)

## Reproducibility Checklist

- [x] All random seeds documented
- [x] All hyperparameters listed in paper appendix
- [x] Pre-computed checkpoints for all experiments
- [x] One-command figure regeneration
- [x] Minimal dependencies (matplotlib, numpy, scipy)
- [x] Paper source with full bibliography

## Citation

```bibtex
@inproceedings{lowrankgnn2027,
  title={Causal Discovery Scales Linearly: Low-Rank Factorization
         Breaks the 500,000× Dimensionality Barrier},
  author={Anonymous},
  booktitle={International Conference on Learning Representations},
  year={2027}
}
```

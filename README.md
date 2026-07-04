# Causal Discovery Scales Linearly

**Low-Rank Factorization Breaks the 500,000x Dimensionality Barrier**

ICLR 2027 Submission (Anonymous)

---

## Abstract

We show that causal discovery can scale linearly with the number of variables. Factorizing the adjacency matrix W = UV^T (r << d) reduces complexity from O(d^3) to O(dr^2). On a single consumer GPU, the method scales to d = 100,000,000 -- a 500,000x expansion beyond the prior limit (d = 200), while maintaining F1 > 0.98 on synthetic benchmarks.

## Reproduction

All 8 figures can be regenerated from pre-computed checkpoints:

```bash
# Regenerate all figures (requires Python 3.12+, PyTorch 2.11+)
python scripts/gen_fig1_architecture.py
python scripts/gen_fig2_benchmark.py
python scripts/gen_fig3_validation.py
python scripts/gen_fig4_sensitivity.py
python scripts/gen_fig5_ablation.py
python scripts/gen_fig6_failure.py
python scripts/gen_fig_appendix_A1_extreme.py
python scripts/gen_fig_appendix_A2_multiscale.py
```

Figures 1-5 complete in seconds. Figure A2 trains for 800 epochs (~2 min on GPU). No downloads required; all data in `checkpoints/`.

## Repository Structure

```
├── main.pdf                                  # Compiled paper
├── main.tex                                  # LaTeX source
├── refs.bib                                  # Bibliography (38 verified references)
├── figures/                                  # All 8 figures (PNG + PDF)
│   ├── fig1_architecture.*
│   ├── fig2_benchmark.*
│   ├── fig3_validation.*
│   ├── fig4_sensitivity.*
│   ├── fig5_ablation.*
│   ├── fig6_failure.*
│   ├── fig_appendix_A1_extreme_scale.*
│   └── fig_appendix_A2_multiscale_heatmap.*
├── scripts/                                  # Figure generation scripts (8 files)
├── checkpoints/                              # Pre-computed experimental results
│   ├── sota_bench.json                       # NOTEARS/DAGMA/GOLEM benchmarks
│   ├── extreme_scale.json                    # d=2K to 100M scaling
│   ├── rank_violation.json                   # Rank misspecification
│   ├── multiscale_ablation.json              # Multi-scale decomposition
│   ├── sachs_*.json                          # Protein signaling validation
│   └── ...
└── LowRankGNN_Replication_Package.zip        # Full reproduction archive
```

## Key Results

| Method | d_max | F1 @ d=200 | Time @ d=100 |
|:--|:--|:--|:--|
| NOTEARS | 200 | 0.001 | 984s |
| DAGMA | 150 | 0.000 | 29.6s |
| GOLEM | 150 | 0.023 | 67.4s |
| **Ours** | **100,000,000** | **0.991** | **0.1s** |

| Application | Metric |
|:--|:--|
| Genome-scale recovery (d=19,215) | 89.1%, 26 min |
| CRISPR dependency prediction | r = 0.912 |
| Drug sensitivity (1,482 compounds) | r = 0.865 |
| TRRUST validation | 94/94 edges, precision = 1.00 |
| TCGA pan-cancer | 33/33 cancers |
| Extreme scale | d=100M, 738s on consumer GPU |

## Citation

```bibtex
@inproceedings{lowrankgnn2027,
  title={Causal Discovery Scales Linearly: Low-Rank Factorization
         Breaks the 500,000x Dimensionality Barrier},
  author={Anonymous},
  booktitle={International Conference on Learning Representations},
  year={2027}
}
```

## License

Code and data released for reproducibility. All checkpoints are pre-computed and deterministic (seed=42).

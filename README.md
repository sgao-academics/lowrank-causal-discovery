# A causal network nominates an inflammatory regulatory core in Indian oral cancer

Replication package for the research article

> *An inflammatory regulatory core is activated in Indian oral squamous cell carcinoma:
> a causal-network-guided analysis of two independent patient cohorts.*

The package contains the causal network, the result files behind every number and
figure in the paper, the scripts that generate them, and the external-validation
analyses. Raw third-party datasets are not redistributed here; they are all public
and are listed under **Data sources** below.

---

## What is in here

| Path | Contents |
| :--- | :--- |
| `checkpoints/edge_list.csv` | The released directed network: 28,247 edges over 2,879 genes, with a `source` column giving the provenance of each edge (`discovered`, `TRRUST`, `pathway_*`, `crosstalk`). |
| `checkpoints/*.json` | All numeric results reported in the text and supplement: degree structure, synthetic benchmarks, rank diagnostics, and the MSigDB hallmark enrichment test. |
| `scripts/_build_cs_refs.py` | Builds the numbered reference list in order of first citation. |
| `scripts/_compile_cs.py`, `scripts/_audit_cs.py` | Compile the manuscript and audit it (word counts, float order, page count, undefined citations). |
| `scripts/_bench_partition_notears.py` | Synthetic Erdős–Rényi benchmark: low-rank estimator vs explicit-adjacency methods. |
| `figure_scripts/_gen_figs_oral.py` | Generates Figures 1-3 and Figure S2 from `checkpoints/edge_list.csv` and the cohort matrices. |
| `figures_cs/` | The generated figures (PDF and PNG). |
| `scripts/validation/` | The ten external-validation and robustness analyses (see below). |
| `validation_outputs/` | The raw console output of each validation script, as run. |

### Validation scripts

| Script | Analysis reported in |
| :--- | :--- |
| `step1_two_cohort_core.py` | Benjamini–Hochberg recomputation of the two Indian cohorts, and the both-cohort core |
| `step2_tcga_hnsc_replication.py` | Section S8 — replication in an independent non-Indian cohort (TCGA-HNSC) |
| `step3_purity_deconvolution.py` | Section S9 — the core survives adjustment for immune and stromal content |
| `step4a_singlecell_fetch.py`, `step4b_singlecell_localisation.py` | Section S10 — single-cell localisation of the core (GSE103322) |
| `step5a_drug_response_gdsc2.py` | Section S11 — core score vs drug response (GDSC2) |
| `step5b_dependency_depmap.py` | Section S11 — CRISPR gene dependency (DepMap), reported as a negative result |
| `step6_hallmark_enrichment.py` | Section S3 — MSigDB hallmark enrichment of the estimator-proposed edges, against uniform and degree-matched nulls |
| `step7_degree_vs_differential_expression.py` | Section S6 — what the degree ranking does and does not predict |
| `step8_gene_set_enrichment.py` | Section S4 — gene-set enrichment with odds ratios and confidence intervals |
| `step9_pancancer_specificity.py` | Section S12 — the core across all 21 TCGA cohorts with adjacent normal tissue |
| `step10_multiomics_cnv_methylation.py` | Section S13 — GISTIC copy number and 450k promoter methylation of the core genes (TCGA-HNSC) |

---

## Data sources

Download these yourself; none is redistributed here.

| Resource | Where |
| :--- | :--- |
| GSE85195 (34 gingivobuccal carcinomas, 15 leukoplakias; GPL6480) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85195 |
| GSE23558 (27 oral tumours, 5 normal; GPL6480) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23558 |
| GPL6480 platform annotation | GEO, `GPL6480.annot.gz` |
| GSE103322 (head-and-neck single-cell RNA-seq) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322 |
| DepMap 24Q2 (CRISPR Chronos gene effect, expression, model table) | https://depmap.org/portal/data_page/ |
| GDSC2 drug-response dataset | https://www.cancerrxgene.org/downloads |
| MSigDB hallmark gene sets | https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp |
| TCGA-HNSC expression (pan-cancer normalised release) | https://xenabrowser.net/datapages/ |
| TCGA-HNSC gene-level GISTIC copy number and HumanMethylation450 beta values | UCSC Xena, `TCGA.HNSC.sampleMap` (https://xenabrowser.net/datapages/) |
| Illumina HumanMethylation450 manifest v1.2 (probe-to-gene and region annotation) | https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv |
| TRRUST v2 | https://www.grnpedia.org/trrust/ |

---

## Setting paths

The scripts use two placeholders so that no local path is hard-coded:

* `{DATA_ROOT}` — the directory holding the downloaded public datasets. The scripts
  expect a layout such as
  `{DATA_ROOT}/cancer_application/data/indian_oral/` (the two cohort series matrices
  and `indian_genomewide_de.json`), `{DATA_ROOT}/single_cell/`, and so on.
  `step9_pancancer_specificity.py` additionally expects the 33 per-cohort
  `TCGA_<CODE>_HiSeqV2.tsv` matrices under `{DATA_ROOT}/cancer_application/data/`.
  `step10_multiomics_cnv_methylation.py` downloads its two Xena files and the
  Illumina manifest into `{DATA_ROOT}/cancer_application/data/hnsc_multiomics/`
  on first run.
* Paths beginning with `.` are relative to this repository.

Run every script from the repository root.

---

## Requirements

Python 3.12 with NumPy, SciPy, pandas, matplotlib and PyTorch. See
`requirements.txt`. The network itself was estimated on a single workstation
(AMD Ryzen 9 8945HX, 32 GB RAM, one NVIDIA RTX 5060 with 8 GB); none of the
validation scripts needs a GPU.

---

## Notes on reproducibility

* The network edge list and the result files are included, so every number in the
  paper can be checked without refitting the network.
* Synthetic benchmarks in `scripts/_bench_partition_notears.py` use fixed seeds.
* The exact plot style used for the figures is defined inside
  `figure_scripts/_gen_figs_oral.py`.
* `.gitignore` excludes publisher PDFs of the cited literature; those are not
  redistributed for copyright reasons.

---

## Licence

The code and derived data in this repository are released for academic reuse.
Third-party datasets remain under the terms of their original providers.

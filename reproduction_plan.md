# Reproduction Plan: scRNA-seq Analysis of Advanced NSCLC

<!-- Badges -->
![Snakemake](https://img.shields.io/badge/workflow-Snakemake-brightgreen)
![Python](https://img.shields.io/badge/python-%3E%3D3.9-blue)
![Platform](https://img.shields.io/badge/platform-HPC%20%7C%20SLURM-orange)
![Status](https://img.shields.io/badge/status-in%20progress-yellow)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

---

## Paper Reference

> **Wu et al., 2021** — *Single-cell profiling of tumor heterogeneity and the microenvironment in advanced non-small cell lung cancer*  
> *Nature Communications* 12, 2540 (2021)  
> DOI: [10.1038/s41467-021-22801-0](https://doi.org/10.1038/s41467-021-22801-0)  
> GEO Accession: [GSE148071](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148071)

---

## Quickstart

```bash
# 1. Clone the repository
git clone https://github.com/d1stadeyemi/scrna_GSE148071_project.git
cd scrna_GSE148071_project

# 2. Set up the base conda environment
conda env create -f workflow/envs/base.yaml
conda activate scrna-nsclc

# 3. Configure cluster profile
cp config/cluster_template.yaml config/slurm/config.yaml

# 4. Dry run to validate the workflow
snakemake --profile config/slurm/config.yaml -n

# 5. Execute the full pipeline
snakemake --profile config/slurm/config.yaml --cores all
```

> **Note:** GEO data download is handled automatically by the pipeline (Stage 1). No manual data preparation is required.

---

## Project Overview

This repository reproduces the core single-cell RNA-seq analyses from Wu et al. (2021), a study profiling tumor heterogeneity and the tumor microenvironment (TME) in advanced non-small cell lung cancer (NSCLC).

The pipeline is implemented from scratch in Python using the **Scanpy ecosystem**, translating the original Seurat-based workflow. It is fully automated via **Snakemake** and designed for execution on an **HPC cluster**.

This project serves as both a rigorous scientific reproduction and a portfolio demonstration of production-grade single-cell analysis.

---

## Objectives

### Primary Goals
- Reproduce the following core analyses from the paper:
  - Quality control and preprocessing
  - Dimensionality reduction (PCA, UMAP)
  - Cell clustering (Leiden algorithm)
  - Cell type annotation using marker genes
  - Tumor heterogeneity analysis via CNV inference
  - Immune subpopulation analysis (T cells, myeloid cells)
  - Trajectory / pseudotime analysis
  - Cell-cell interaction analysis *(conditional — see Stage 10)*

### Secondary Goals
- Provide a clean Scanpy-based alternative to the original Seurat workflow
- Build a modular, documented, and reusable pipeline template
- Produce a portfolio-ready GitHub repository

---

## Reproducibility Scope

| Included | Excluded |
|---|---|
| All analyses from GEO-provided processed data | Raw FASTQ processing |
| Full Scanpy reimplementation of the analysis | Alignment and CellRanger steps |
| Automated pipeline from GEO download to figures | Proprietary or restricted data |
| Parameter documentation and justification | Exact numerical reproduction (Seurat vs Scanpy differences expected) |

> Exact numerical reproduction is not expected due to differences between Seurat (original) and Scanpy implementations. Biological concordance with reported findings is the target criterion.

---

## Data Sources

**Primary Dataset**
- GEO Accession: [GSE148071](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148071)
- Files used:
  - `GSE148071_RAW.tar` — processed per-sample expression matrices (TXT format)
  - Series Matrix File — sample-level metadata
- Raw FASTQ files are not publicly available

**Notes on Data Format**
- Per-sample expression matrices require format standardization before merging
- Metadata requires reconstruction from Series Matrix and/or supplementary tables
- Cell barcode and sample identity must be explicitly tracked across the merge

---

## Project Structure

```
project/
├── config/
│   ├── config.yaml              # Pipeline parameters
│   ├── cluster.yaml             # HPC SLURM profile (user-configured)
│   └── cluster_template.yaml   # Template for cluster configuration
├── workflow/
│   ├── Snakefile                # Main workflow entry point
│   ├── rules/                   # Modular rule files per stage
│   │   ├── download.smk
│   │   ├── parse.smk
│   │   ├── preprocess.smk
│   │   ├── dimred.smk
│   │   ├── cluster.smk
│   │   ├── annotate.smk
│   │   ├── cnv.smk
│   │   ├── subpopulation.smk
│   │   ├── trajectory.smk
│   │   └── interactions.smk
│   └── envs/                    # Per-rule conda environments
│       ├── base.yaml
│       ├── infercnvpy.yaml
│       └── scvelo.yaml
├── scripts/                     # Python scripts called by rules
├── data/
│   ├── raw/                     # Downloaded GEO files (unmodified)
│   ├── processed/               # Intermediate AnnData objects (.h5ad)
│   └── metadata/                # Parsed and reconstructed metadata
├── results/
│   ├── figures/                 # All output plots
│   └── tables/                  # Marker gene tables, cell type assignments
├── logs/                        # Per-rule HPC logs
├── resources/                   # Reference gene lists, marker panels
├── notebooks/                   # Exploratory analysis and QC review
└── README.md
```

---

## Computational Environment

| Component | Specification |
|---|---|
| Execution platform | HPC cluster (SLURM scheduler) |
| Workflow manager | Snakemake ≥ 7.0 |
| Language | Python ≥ 3.9 |
| Environment manager | Conda (per-rule environments) |

**Core dependencies**

| Package | Purpose |
|---|---|
| `scanpy` | Core scRNA-seq analysis |
| `anndata` | Data container |
| `infercnvpy` | CNV inference for tumor cell identification |
| `scvi-tools` | Optional: deep generative models / batch correction |
| `scvelo` | RNA velocity and trajectory analysis |
| `squidpy` | Optional: spatial/cell-cell interaction analysis |
| `liana` | Ligand-receptor interaction analysis |
| `matplotlib` / `seaborn` | Visualization |

---

## Workflow Stages

### Stage 1 — Data Acquisition
- Download `GSE148071_RAW.tar` and Series Matrix from GEO
- Verify file integrity (MD5 checksums)
- Extract archive and generate file inventory
- **Expected output:** Raw data in `data/raw/`, checksum log

### Stage 2 — Data Parsing
- Identify and inspect all per-sample expression matrices
- Standardize to genes × cells orientation
- Convert to AnnData (`.h5ad`) format
- Construct unified metadata table with sample annotations
- **Expected output:** `data/processed/merged_raw.h5ad`

### Stage 3 — Preprocessing
- Cell-level QC: filter on n_genes, n_counts, % mitochondrial reads
- Gene-level QC: remove low-expressed genes
- Normalization (scran or scanpy default)
- Log1p transformation
- Highly variable gene (HVG) selection
- **Expected output:** `data/processed/preprocessed.h5ad`, QC violin plots

### Stage 4 — Dimensionality Reduction
- PCA (top 50 PCs)
- Batch-aware neighborhood graph (Harmony or scVI if batch effects detected)
- UMAP embedding
- **Expected output:** `data/processed/dimred.h5ad`, PCA elbow plot, UMAP

### Stage 5 — Clustering
- Leiden clustering across multiple resolutions (0.3 – 1.5)
- Silhouette score and cluster stability assessment
- Select final resolution matching paper's reported cluster count
- **Expected output:** Cluster UMAP, resolution sweep summary

### Stage 6 — Cell Type Annotation
- Marker gene-based annotation using paper's Supplementary Table
- Dotplot and violin plots per cluster
- Map to paper-defined cell type labels
- **Expected output:** Annotated UMAP, marker dotplots, annotation table

### Stage 7 — Tumor Heterogeneity Analysis
- CNV inference using `infercnvpy` (epithelial/immune reference)
- Classify malignant vs non-malignant epithelial cells
- CNV score distribution across clusters
- **Expected output:** CNV heatmap, malignant cell annotation, updated `.h5ad`

### Stage 8 — Subpopulation Analysis
- Extract and recluster T cells, NK cells, myeloid cells separately
- Annotate subpopulations (CD8+ T effector, Treg, TAM, etc.)
- Compute cell type proportions per patient and treatment group
- **Expected output:** Subpopulation UMAPs, proportion barplots

### Stage 9 — Trajectory Analysis
- Pseudotime inference (Scanpy PAGA + diffusion pseudotime)
- Optional: RNA velocity with `scVelo` if spliced/unspliced counts available
- Compare inferred trajectories to paper-reported differentiation paths
- **Expected output:** Trajectory plots, pseudotime-ordered gene expression

### Stage 10 — Cell-Cell Interaction Analysis *(conditional)*
- Ligand-receptor analysis using `LIANA`
- Visualize interaction networks between major cell types
- **Condition:** Implemented only if Stage 6–8 annotations meet quality threshold
- **Expected output:** Interaction heatmaps, chord diagrams

---

## Resource Estimates (HPC)

| Stage | CPUs | Memory | Estimated Wall Time |
|---|---|---|---|
| Data Parsing | 4 | 32 GB | ~30 min |
| Preprocessing | 8 | 64 GB | ~1 hr |
| Dimensionality Reduction | 8 | 64 GB | ~1–2 hr |
| Clustering | 4 | 32 GB | ~30 min |
| CNV Inference (Stage 7) | 16 | 128 GB | ~3–6 hr |
| Trajectory Analysis | 8 | 64 GB | ~1–2 hr |

> Estimates based on ~200,000 cells (approximate dataset size). Adjust `config/cluster.yaml` accordingly.

---

## Target Outputs

| Output | Description |
|---|---|
| `results/figures/umap_annotated.pdf` | Full UMAP with cell type annotations |
| `results/figures/umap_subpopulations/` | Per-lineage subpopulation UMAPs |
| `results/figures/cnv_heatmap.pdf` | CNV inference heatmap |
| `results/figures/cell_proportions.pdf` | Cell type proportions per sample |
| `results/figures/trajectory.pdf` | Pseudotime trajectory plot |
| `results/tables/cluster_markers.csv` | Per-cluster marker genes |
| `results/tables/cell_type_assignments.csv` | Final cell → type mapping |
| `data/processed/final_annotated.h5ad` | Final AnnData object |

---

## Success Criteria

The reproduction is considered successful when:

- [ ] Major cell types from the paper are recovered (T cells, myeloid, epithelial, fibroblasts, endothelial)
- [ ] Tumor vs non-tumor epithelial cells are distinguishable via CNV
- [ ] Known immune subpopulations are detected (CD8+ T, Treg, exhausted T, TAM)
- [ ] Pipeline runs end-to-end from a fresh GEO download with a single Snakemake command
- [ ] Results are reproducible on a clean conda environment
- [ ] All figures are saved to `results/figures/` with corresponding config parameters logged

---

## Risks and Mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| Incomplete / inconsistent GEO metadata | Medium | Manual curation from paper supplementary; document all assumptions |
| Seurat → Scanpy implementation differences | High | Target biological concordance, not numerical identity; document parameter choices |
| Batch effects across patients | High | Test Harmony and scVI integration; compare with and without correction |
| CNV inference sensitivity | Medium | Tune reference cell population; validate against known tumor markers |
| Missing spliced/unspliced counts for scVelo | High | Fall back to PAGA + diffusion pseudotime |
| Parameter sensitivity in clustering | Medium | Sweep resolutions; report stability metrics |

---

## Future Extensions

- Integrate ML-based cell classifiers (scANVI, CellTypist) as an alternative annotation strategy
- Apply to Fan et al., 2025 dataset for cross-study comparison
- Extend with spatial transcriptomics context if publicly available data exists
- Add scVI/deep generative model integration as a dedicated optional stage

---

## Design Principles

- **Reproducibility over convenience** — every parameter is explicit and logged
- **Modular by design** — each stage is an independent Snakemake rule
- **Fail loudly** — QC checkpoints at each stage with logged thresholds
- **Separation of concerns** — data, code, and results are strictly isolated
- **Documented assumptions** — all deviations from the original paper are noted inline

---

## Citation

If you use or adapt this pipeline, please cite the original study:

```
Wu, F. et al. Single-cell profiling of tumor heterogeneity and the microenvironment 
in advanced non-small cell lung cancer. Nat Commun 12, 2540 (2021). 
https://doi.org/10.1038/s41467-021-22801-0
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.
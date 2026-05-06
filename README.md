# scRNA-seq Reproduction: Wu et al. 2021 — Advanced NSCLC

Reproduction of key analyses from Wu et al. (2021), 
*Nature Communications*, using the Scanpy ecosystem on HPC.

GEO: GSE148071 | DOI: 10.1038/s41467-021-22801-0

## Status
🚧 In progress — pipeline under construction

## What this repo will contain
- Full Snakemake pipeline (GEO download → figures)
- Scanpy reimplementation of the original Seurat workflow
- CNV inference, subpopulation analysis, trajectory analysis

## Quickstart
_Coming soon — environment and pipeline not yet finalized_

## Repo Structure
_Coming soon_

## See also
- [reproduction_plan.md](reproduction_plan.md) — detailed analysis plan and scope
## Required Reference Files

The following files must be downloaded manually before running Stage 7 (CNV inference):

```bash
mkdir -p resources/genome
wget -O resources/genome/Homo_sapiens.GRCh38.92.gtf.gz \
    https://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz
gunzip resources/genome/Homo_sapiens.GRCh38.92.gtf.gz
```

This is Ensembl GRCh38 v92 — the same annotation used by Wu et al. 2021.

"""
Stage 3 — Quality Control and Filtering

Input:  data/processed/merged_raw.h5ad
Output: data/processed/qc_filtered.h5ad
        results/figures/qc/

QC strategy (following Wu et al. 2021 and standard Scanpy best practices):
- Calculate per-cell: n_genes, n_counts, pct_counts_mt
- Visualize distributions per sample before filtering
- Filter cells: min/max genes, max % mitochondrial reads
- Filter genes: minimum cells expressing the gene
- Save filtered AnnData and all QC figures
"""

import os
import argparse
import logging
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use("Agg")   # non-interactive backend for HPC
import matplotlib.pyplot as plt
import seaborn as sns

# ── Logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s [%(levelname)s] %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger(__name__)


def calculate_qc_metrics(adata):
    """Calculate standard QC metrics and add to adata.obs and adata.var."""
    log.info("Calculating QC metrics...")

    # Flag mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars        = ["mt"],
        percent_top    = None,
        log1p          = False,
        inplace        = True
    )

    log.info(f"  Median genes per cell:   {adata.obs['n_genes_by_counts'].median():.0f}")
    log.info(f"  Median counts per cell:  {adata.obs['total_counts'].median():.0f}")
    log.info(f"  Median MT%:              {adata.obs['pct_counts_mt'].median():.2f}%")
    return adata


def plot_qc_metrics(adata, figures_dir):
    """Generate QC violin and scatter plots before filtering."""
    log.info("Generating pre-filtering QC plots...")
    os.makedirs(figures_dir, exist_ok=True)

    # ── Violin plots — overall distributions ─────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    metrics = [
        ("n_genes_by_counts",  "Genes per cell",         "steelblue"),
        ("total_counts",       "UMI counts per cell",    "seagreen"),
        ("pct_counts_mt",      "% Mitochondrial reads",  "tomato"),
    ]
    
    for ax, (metric, label, color) in zip(axes, metrics):
        ax.violinplot(adata.obs[metric], showmedians=True)
        ax.set_title(label)
        ax.set_ylabel(label)
        ax.set_xticks([])
    
    plt.suptitle("Pre-filtering QC metrics (all cells)", fontsize=13)
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "qc_violin_overall.pdf"), dpi=150)
    plt.close()

    # ── Violin plots — per patient ────────────────────────────────────────────
    fig, axes = plt.subplots(3, 1, figsize=(20, 12))
    patients  = sorted(adata.obs["patient_id"].unique(),
                       key=lambda x: int(x.replace("P", "")))

    for ax, (metric, label, color) in zip(axes, metrics):
        data_per_patient = [
            adata.obs.loc[adata.obs["patient_id"] == p, metric].values
            for p in patients
        ]
        parts = ax.violinplot(data_per_patient, showmedians=True)
        for pc in parts["bodies"]:
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
        ax.set_xticks(range(1, len(patients) + 1))
        ax.set_xticklabels(patients, rotation=90, fontsize=7)
        ax.set_ylabel(label)
        ax.set_title(f"{label} per patient")

    plt.suptitle("Pre-filtering QC metrics per patient", fontsize=13)
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "qc_violin_per_patient.pdf"), dpi=150)
    plt.close()

    # ── Scatter: counts vs genes, colored by MT% ──────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    sc_kwargs = dict(
        x      = adata.obs["total_counts"],
        c      = adata.obs["pct_counts_mt"],
        cmap   = "RdYlGn_r",
        alpha  = 0.3,
        s      = 1,
    )

    sc1 = axes[0].scatter(y=adata.obs["n_genes_by_counts"], **sc_kwargs)
    axes[0].set_xlabel("Total counts (UMI)")
    axes[0].set_ylabel("Genes detected")
    axes[0].set_title("Counts vs Genes (colored by MT%)")
    plt.colorbar(sc1, ax=axes[0], label="MT%")

    axes[1].scatter(
        x     = adata.obs["total_counts"],
        y     = adata.obs["pct_counts_mt"],
        alpha = 0.3,
        s     = 1,
        color = "tomato"
    )
    axes[1].set_xlabel("Total counts (UMI)")
    axes[1].set_ylabel("% Mitochondrial reads")
    axes[1].set_title("Counts vs MT%")

    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "qc_scatter_counts_genes_mt.pdf"), dpi=150)
    plt.close()

    # ── Cell count per patient (bar) ──────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(14, 4))
    counts = adata.obs["patient_id"].value_counts()
    counts = counts.reindex(patients)
    ax.bar(patients, counts.values, color="steelblue", edgecolor="white")
    ax.set_xlabel("Patient")
    ax.set_ylabel("Cell count")
    ax.set_title("Cells per patient (pre-filtering)")
    ax.axhline(500, color="red", linestyle="--", linewidth=1, label="500 cell threshold")
    plt.xticks(rotation=90, fontsize=7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "qc_cells_per_patient_prefilt.pdf"), dpi=150)
    plt.close()

    log.info(f"  Saved QC plots to {figures_dir}")


def filter_cells_and_genes(adata, min_genes, max_genes, max_pct_mt, min_cells):
    """Apply QC filters to cells and genes."""
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    log.info("Filtering cells...")
    log.info(f"  Min genes per cell:  {min_genes}")
    log.info(f"  Max genes per cell:  {max_genes}")
    log.info(f"  Max MT%:             {max_pct_mt}%")

    sc.pp.filter_cells(adata, min_genes=min_genes)
    log.info(f"  After min_genes filter: {adata.n_obs:,} cells remaining")

    adata = adata[adata.obs["n_genes_by_counts"] < max_genes].copy()
    log.info(f"  After max_genes filter: {adata.n_obs:,} cells remaining")

    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mt].copy()
    log.info(f"  After MT% filter:       {adata.n_obs:,} cells remaining")

    log.info("Filtering genes...")
    sc.pp.filter_genes(adata, min_cells=min_cells)
    log.info(f"  After min_cells filter: {adata.n_vars:,} genes remaining")

    log.info(f"Cells removed: {n_cells_before - adata.n_obs:,} "
             f"({(n_cells_before - adata.n_obs) / n_cells_before * 100:.1f}%)")
    log.info(f"Genes removed: {n_genes_before - adata.n_vars:,} "
             f"({(n_genes_before - adata.n_vars) / n_genes_before * 100:.1f}%)")

    return adata


def plot_qc_post_filter(adata, figures_dir):
    """Plot cell counts per patient after filtering."""
    log.info("Generating post-filtering QC plots...")

    patients = sorted(adata.obs["patient_id"].unique(),
                      key=lambda x: int(x.replace("P", "")))

    fig, ax = plt.subplots(figsize=(14, 4))
    counts = adata.obs["patient_id"].value_counts().reindex(patients)
    ax.bar(patients, counts.values, color="seagreen", edgecolor="white")
    ax.set_xlabel("Patient")
    ax.set_ylabel("Cell count")
    ax.set_title("Cells per patient (post-filtering)")
    plt.xticks(rotation=90, fontsize=7)
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "qc_cells_per_patient_postfilt.pdf"), dpi=150)
    plt.close()


def main(input_path, output_path, figures_dir,
         min_genes, max_genes, max_pct_mt, min_cells):

    log.info(f"Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    log.info(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # QC metrics
    adata = calculate_qc_metrics(adata)

    # Pre-filter plots
    plot_qc_metrics(adata, figures_dir)

    # Filter
    adata = filter_cells_and_genes(
        adata,
        min_genes  = min_genes,
        max_genes  = max_genes,
        max_pct_mt = max_pct_mt,
        min_cells  = min_cells,
    )

    # Post-filter plot
    plot_qc_post_filter(adata, figures_dir)

    # Log final state
    log.info(f"Final AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    log.info(f"Patients retained: {adata.obs['patient_id'].nunique()}")

    # Save
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    log.info(f"Saved filtered AnnData to {output_path}")
    log.info("Stage 3 complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="QC and filter scRNA-seq AnnData")
    parser.add_argument("--input",       default="data/processed/merged_raw.h5ad")
    parser.add_argument("--output",      default="data/processed/qc_filtered.h5ad")
    parser.add_argument("--figures_dir", default="results/figures/qc")
    parser.add_argument("--min_genes",   type=int,   default=200)
    parser.add_argument("--max_genes",   type=int,   default=6000)
    parser.add_argument("--max_pct_mt",  type=float, default=25.0)
    parser.add_argument("--min_cells",   type=int,   default=3)
    args = parser.parse_args()

    main(
        input_path  = args.input,
        output_path = args.output,
        figures_dir = args.figures_dir,
        min_genes   = args.min_genes,
        max_genes   = args.max_genes,
        max_pct_mt  = args.max_pct_mt,
        min_cells   = args.min_cells,
    )
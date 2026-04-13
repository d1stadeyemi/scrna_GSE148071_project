"""
Stage 4a — Normalization and HVG Selection

Input:  data/processed/qc_filtered.h5ad
Output: data/processed/normalized.h5ad

Steps:
- Normalize total counts per cell (target sum 10,000)
- Log1p transform
- Select highly variable genes
- PCA on HVG space
- Store raw counts in adata.layers['counts'] for scVI
"""

import os
import argparse
import logging
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s [%(levelname)s] %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger(__name__)


def main(input_path, output_path, figures_dir, n_hvgs, n_pcs):

    log.info(f"Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    log.info(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ── Store raw counts for scVI (requires integer counts) ──────────────────
    # scVI needs raw counts — store before normalization
    adata.layers["counts"] = adata.X.copy()
    log.info("Stored raw counts in adata.layers['counts']")

    # ── Normalize ────────────────────────────────────────────────────────────
    log.info("Normalizing total counts per cell (target sum = 10,000)...")
    sc.pp.normalize_total(adata, target_sum=1e4)

    log.info("Log1p transforming...")
    sc.pp.log1p(adata)

    # Store normalized+log1p in a layer for reference
    adata.layers["log1p_norm"] = adata.X.copy()

    # ── Highly variable genes ─────────────────────────────────────────────────
    log.info(f"Selecting top {n_hvgs} highly variable genes...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes = n_hvgs,
        batch_key   = "patient_id",   # account for patient batch during HVG selection
        flavor      = "seurat_v3",
        layer       = "counts",       # use raw counts for seurat_v3 method
    )

    n_hvg_found = adata.var["highly_variable"].sum()
    log.info(f"  HVGs selected: {n_hvg_found:,}")

    # ── Plot HVG selection summary ────────────────────────────────────────────
    os.makedirs(figures_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    hvg_mask = adata.var["highly_variable"]
    colors = hvg_mask.map({True: "tomato", False: "lightgrey"})
    ax.scatter(
        range(adata.n_vars),
        adata.var["mean_counts"].values,
        c     = colors.values,
        s     = 1,
        alpha = 0.5,
    )
    ax.set_xlabel("Gene rank (by mean expression)")
    ax.set_ylabel("Mean counts")
    ax.set_title(f"HVG selection: {hvg_mask.sum():,} HVGs (red) / {adata.n_vars:,} total genes")
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0],[0], marker='o', color='w', markerfacecolor='tomato',
            markersize=6, label=f'HVG (n={hvg_mask.sum():,})'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='lightgrey',
            markersize=6, label='Other genes'),
    ])
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "hvg_selection.pdf"), dpi=150)
    plt.close()

    # ── PCA ──────────────────────────────────────────────────────────────────
    log.info(f"Running PCA (n_comps={n_pcs}) on HVG subset...")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True)

    # Elbow plot
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(
        range(1, n_pcs + 1),
        adata.uns["pca"]["variance_ratio"],
        marker = "o",
        ms     = 3,
    )
    ax.set_xlabel("PC")
    ax.set_ylabel("Variance ratio")
    ax.set_title("PCA elbow plot")
    ax.axvline(30, color="red", linestyle="--", linewidth=1, label="n_pcs=30")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "pca_elbow.pdf"), dpi=150)
    plt.close()

    log.info(f"Final AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    log.info(f"  HVG flag stored in adata.var['highly_variable']")
    log.info(f"  PCA stored in adata.obsm['X_pca']")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    log.info(f"Saved to {output_path}")
    log.info("Stage 4a complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",       default="data/processed/qc_filtered.h5ad")
    parser.add_argument("--output",      default="data/processed/normalized.h5ad")
    parser.add_argument("--figures_dir", default="results/figures/hvg_pca")
    parser.add_argument("--n_hvgs",      type=int, default=3000)
    parser.add_argument("--n_pcs",       type=int, default=50)
    args = parser.parse_args()
    main(args.input, args.output, args.figures_dir, args.n_hvgs, args.n_pcs)
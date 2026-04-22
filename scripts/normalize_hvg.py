"""
Stage 4a — Normalization, HVG Selection, and PCA

Shared by both tracks. Track-specific parameters passed as arguments:
  paper_track:       600 HVGs, seurat flavor, 20 PCs
  alternative_track: 3000 HVGs, seurat_v3 flavor, 50 PCs

Biological steps:
1. Store raw counts before normalization (scVI requires raw counts)
2. Normalize to 10,000 counts per cell (CPM-like, standard)
3. Log1p transform (stabilizes variance, approximates normal distribution)
4. Select HVGs (genes informative across cell types and patients)
5. PCA on HVG subset (linear dimensionality reduction)

HVG flavor choice:
  seurat:    Fano-factor dispersion (Seurat 2.3 method, paper-faithful)
  seurat_v3: Variance-stabilizing, batch-aware (modern, alternative track)
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
from matplotlib.lines import Line2D

logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s [%(levelname)s] %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger(__name__)


def main(input_path, output_path, figures_dir, n_hvgs, n_pcs,
         hvg_flavor, track):

    log.info(f"[{track.upper()} TRACK] Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    log.info(f"  Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ── Store raw counts ──────────────────────────────────────────────────────
    # scVI requires integer counts — store before any normalization.
    # Also useful as reference for differential expression.
    adata.layers["counts"] = adata.X.copy()
    log.info("  Raw counts stored in adata.layers['counts']")

    # ── Normalize ────────────────────────────────────────────────────────────
    # Target sum = 10,000: standard CPM-like normalization.
    # Makes cells comparable regardless of sequencing depth.
    log.info("  Normalizing to 10,000 counts per cell...")
    sc.pp.normalize_total(adata, target_sum=1e4)

    log.info("  Log1p transforming...")
    sc.pp.log1p(adata)
    adata.layers["log1p_norm"] = adata.X.copy()

    # ── Highly Variable Genes ─────────────────────────────────────────────────
    log.info(f"  Selecting top {n_hvgs} HVGs using '{hvg_flavor}' method...")

    if hvg_flavor == "seurat_v3":
        # seurat_v3: variance-stabilizing, batch-aware
        # Requires raw counts (not log-normalized) to model mean-variance
        # batch_key ensures HVGs are informative across all patients
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes = n_hvgs,
            flavor      = hvg_flavor,
            batch_key   = "patient_id",
            layer       = "counts",
        )
    else:
        # seurat/cell_ranger: Fano-factor based (paper-faithful)
        # Does not require batch_key — simpler, older method
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes = n_hvgs,
            flavor      = hvg_flavor,
        )

    n_hvg_found = adata.var["highly_variable"].sum()
    log.info(f"  HVGs selected: {n_hvg_found:,}")

    # ── HVG visualization ─────────────────────────────────────────────────────
    os.makedirs(figures_dir, exist_ok=True)
    _plot_hvg(adata, n_hvg_found, hvg_flavor, track, figures_dir)

    # ── PCA ──────────────────────────────────────────────────────────────────
    # Computed on HVG subset only.
    # Scale first (center, cap at 10) to prevent high-expression genes dominating.
    log.info(f"  Running PCA ({n_pcs} components) on HVG subset...")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True)

    _plot_pca_elbow(adata, n_pcs, track, figures_dir)

    log.info(f"  Final: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    log.info(f"  X_pca shape: {adata.obsm['X_pca'].shape}")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    log.info(f"  Saved to {output_path}")
    log.info(f"Stage 4a [{track}] complete.")


def _plot_hvg(adata, n_hvg_found, hvg_flavor, track, figures_dir):
    """Scatter of mean expression vs HVG status."""
    fig, ax = plt.subplots(figsize=(8, 5))
    hvg_mask = adata.var["highly_variable"]
    colors = hvg_mask.map({True: "tomato", False: "lightgrey"})

    ax.scatter(
        range(adata.n_vars),
        adata.var["mean_counts"].values if "mean_counts" in adata.var.columns
        else np.zeros(adata.n_vars),
        c=colors.values, s=1, alpha=0.5, rasterized=True
    )
    ax.set_xlabel("Gene rank", fontsize=10)
    ax.set_ylabel("Mean counts", fontsize=10)
    ax.set_title(
        f"[{track}] HVG selection — {n_hvg_found:,} HVGs / {adata.n_vars:,} genes\n"
        f"Method: {hvg_flavor}",
        fontsize=11, fontweight="bold"
    )
    ax.legend(handles=[
        Line2D([0],[0], marker="o", color="w", markerfacecolor="tomato",
               markersize=6, label=f"HVG (n={n_hvg_found:,})"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor="lightgrey",
               markersize=6, label="Other"),
    ])
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "hvg_selection.pdf"),
                dpi=200, bbox_inches="tight")
    plt.close()


def _plot_pca_elbow(adata, n_pcs, track, figures_dir):
    """PCA elbow plot to assess dimensionality."""
    fig, ax = plt.subplots(figsize=(8, 4))
    var_ratio = adata.uns["pca"]["variance_ratio"]
    ax.plot(range(1, n_pcs + 1), var_ratio, "o-", ms=3, color="steelblue")
    ax.set_xlabel("Principal Component", fontsize=10)
    ax.set_ylabel("Variance Explained", fontsize=10)
    ax.set_title(f"[{track}] PCA elbow plot", fontsize=11, fontweight="bold")
    # Mark the selected number of PCs
    ax.axvline(n_pcs, color="red", linestyle="--", alpha=0.7,
               label=f"Selected: {n_pcs} PCs")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "pca_elbow.pdf"),
                dpi=200, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Normalization, HVG selection, and PCA"
    )
    parser.add_argument("--input",       required=True)
    parser.add_argument("--output",      required=True)
    parser.add_argument("--figures_dir", required=True)
    parser.add_argument("--n_hvgs",      type=int, default=600)
    parser.add_argument("--n_pcs",       type=int, default=20)
    parser.add_argument("--hvg_flavor",  default="seurat",
                        choices=["seurat", "seurat_v3", "cell_ranger"])
    parser.add_argument("--track",       default="paper",
                        choices=["paper", "alternative"])
    args = parser.parse_args()
    main(args.input, args.output, args.figures_dir,
         args.n_hvgs, args.n_pcs, args.hvg_flavor, args.track)

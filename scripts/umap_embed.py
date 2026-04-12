"""
Stage 4c — Neighborhood Graph and UMAP Embedding

Input:  data/processed/scvi_integrated.h5ad
Output: data/processed/embedded.h5ad
        results/figures/umap/

Steps:
- Build neighborhood graph on scVI latent space
- Compute UMAP embedding
- Plot UMAP colored by patient, sample, QC metrics
"""

import os
import argparse
import logging
import warnings
warnings.filterwarnings("ignore")

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


def main(input_path, output_path, figures_dir, n_neighbors, n_pcs):

    log.info(f"Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    log.info(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ── Neighborhood graph on scVI latent ────────────────────────────────────
    log.info(f"Building neighborhood graph (n_neighbors={n_neighbors})...")
    sc.pp.neighbors(
        adata,
        use_rep    = "X_scvi",   # use scVI latent, not PCA
        n_neighbors = n_neighbors,
        n_pcs       = n_pcs,
    )

    # ── UMAP ──────────────────────────────────────────────────────────────────
    log.info("Computing UMAP embedding...")
    sc.tl.umap(adata)
    log.info("UMAP complete.")

    # ── Plots ─────────────────────────────────────────────────────────────────
    os.makedirs(figures_dir, exist_ok=True)
    sc.settings.figdir = figures_dir

    # UMAP colored by patient
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color  = "patient_id",
        title  = "UMAP — colored by patient",
        legend_loc = "on data",
        legend_fontsize = 6,
        ax     = ax,
        show   = False,
    )
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "umap_by_patient.pdf"), dpi=150)
    plt.close()

    # UMAP colored by QC metrics
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, color in zip(axes, ["n_genes_by_counts", "total_counts", "pct_counts_mt"]):
        sc.pl.umap(adata, color=color, ax=ax, show=False, vmax="p99")
    plt.suptitle("UMAP — QC metrics", fontsize=13)
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "umap_qc_metrics.pdf"), dpi=150)
    plt.close()

    log.info(f"Figures saved to {figures_dir}")

    # ── Save ──────────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    log.info(f"Saved to {output_path}")
    log.info("Stage 4c complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",       default="data/processed/scvi_integrated.h5ad")
    parser.add_argument("--output",      default="data/processed/embedded.h5ad")
    parser.add_argument("--figures_dir", default="results/figures/umap")
    parser.add_argument("--n_neighbors", type=int, default=15)
    parser.add_argument("--n_pcs",       type=int, default=30)
    args = parser.parse_args()
    main(args.input, args.output, args.figures_dir, args.n_neighbors, args.n_pcs)
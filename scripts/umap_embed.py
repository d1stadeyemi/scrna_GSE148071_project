"""
Stage 4c — Neighborhood Graph and UMAP Embedding

Shared by both tracks. Key difference is the representation used:
  paper_track:       use_rep="X_pca"   (no batch correction)
  alternative_track: use_rep="X_scvi"  (scVI latent space)

This single parameter change determines whether patient-level variation
is preserved (X_pca) or harmonized (X_scvi).

Biological interpretation of outputs:
- UMAP by patient: should show patient-specific cancer islands (paper track)
  or harmonized cloud (alternative track)
- UMAP by QC metrics: smooth gradients confirm no technical artefacts
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


def main(input_path, output_path, figures_dir, n_neighbors, n_pcs,
         use_rep, track):

    log.info(f"[{track.upper()} TRACK] Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    log.info(f"  Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    log.info(f"  Using representation: {use_rep}")

    # ── Neighborhood graph ────────────────────────────────────────────────────
    # Build kNN graph in the reduced space.
    # n_pcs argument is only used if use_rep="X_pca"
    log.info(f"  Building neighborhood graph (n_neighbors={n_neighbors})...")

    if use_rep == "X_pca":
        sc.pp.neighbors(
            adata,
            use_rep     = "X_pca",
            n_pcs       = n_pcs,
            n_neighbors = n_neighbors,
        )
    elif use_rep == "X_scvi":
        sc.pp.neighbors(
            adata,
            use_rep     = "X_scvi",
            n_neighbors = n_neighbors,
            # n_pcs not used for scVI latent space
        )
    else:
        raise ValueError(f"Unknown use_rep: {use_rep}")

    # ── UMAP ──────────────────────────────────────────────────────────────────
    log.info("  Computing UMAP...")
    sc.tl.umap(adata, min_dist=0.3)
    log.info(f"  UMAP shape: {adata.obsm['X_umap'].shape}")

    # ── Figures ───────────────────────────────────────────────────────────────
    os.makedirs(figures_dir, exist_ok=True)

    _plot_umap_by_patient(adata, track, figures_dir)
    _plot_umap_qc_metrics(adata, track, figures_dir)

    # ── Save ──────────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    log.info(f"  Saved to {output_path}")
    log.info(f"Stage 4c [{track}] complete.")


def _plot_umap_by_patient(adata, track, figures_dir):
    """UMAP colored by patient — shows patient mixing / separation."""
    patients = sorted(
        adata.obs["patient_id"].unique(),
        key=lambda x: int(x.replace("P", ""))
    )
    # Generate 42 visually distinct colors
    cmap1 = plt.cm.get_cmap("tab20", 20)
    cmap2 = plt.cm.get_cmap("tab20b", 20)
    cmap3 = plt.cm.get_cmap("tab20c", 20)
    pool = ([cmap1(i) for i in range(20)] +
            [cmap2(i) for i in range(20)] +
            [cmap3(i) for i in range(4)])
    pcolors = {p: pool[i] for i, p in enumerate(patients)}

    fig, ax = plt.subplots(figsize=(10, 8))
    for p in patients:
        mask = adata.obs["patient_id"] == p
        ax.scatter(
            adata.obsm["X_umap"][mask, 0],
            adata.obsm["X_umap"][mask, 1],
            c=[pcolors[p]], s=0.3, alpha=0.6, rasterized=True
        )
    # Legend in 3 columns
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=pcolors[p], label=p) for p in patients]
    ax.legend(handles=handles, loc="upper left",
              bbox_to_anchor=(1.01, 1), fontsize=7,
              ncol=2, frameon=True, title="Patient", title_fontsize=9)
    ax.set_xlabel("UMAP1", fontsize=11)
    ax.set_ylabel("UMAP2", fontsize=11)
    ax.set_title(f"[{track}] UMAP — colored by patient\n"
                 f"({'No batch correction' if 'pca' in str(adata.uns.get('neighbors', {}).get('params', {}).get('use_rep', '')) else 'scVI batch corrected'})",
                 fontsize=12, fontweight="bold")
    ax.set_aspect("equal")
    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir, f"umap_by_patient.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()


def _plot_umap_qc_metrics(adata, track, figures_dir):
    """UMAP colored by QC metrics — checks for technical artefacts."""
    qc_cols = [c for c in
               ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
               if c in adata.obs.columns]

    if not qc_cols:
        log.warning("  No QC columns found for QC metric UMAP")
        return

    fig, axes = plt.subplots(1, len(qc_cols),
                             figsize=(6 * len(qc_cols), 5))
    if len(qc_cols) == 1:
        axes = [axes]

    for ax, col in zip(axes, qc_cols):
        sc.pl.umap(adata, color=col, ax=ax, show=False, vmax="p99")
        ax.set_title(col, fontsize=9)

    plt.suptitle(f"[{track}] UMAP — QC metrics", fontsize=12,
                 fontweight="bold")
    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir, f"umap_qc_metrics.{ext}"),
                    dpi=200, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build neighborhood graph and compute UMAP"
    )
    parser.add_argument("--input",       required=True)
    parser.add_argument("--output",      required=True)
    parser.add_argument("--figures_dir", required=True)
    parser.add_argument("--n_neighbors", type=int, default=15)
    parser.add_argument("--n_pcs",       type=int, default=20)
    parser.add_argument("--use_rep",     default="X_pca",
                        choices=["X_pca", "X_scvi"],
                        help="Representation for neighborhood graph. "
                             "X_pca = no batch correction (paper track). "
                             "X_scvi = scVI latent (alternative track).")
    parser.add_argument("--track",       default="paper",
                        choices=["paper", "alternative"])
    args = parser.parse_args()
    main(args.input, args.output, args.figures_dir,
         args.n_neighbors, args.n_pcs, args.use_rep, args.track)

"""
Stage 3 — Quality Control Filtering

Applies QC thresholds from Wu et al. 2021 Methods:
  "We removed cells that had either lower than 200 or higher than 5000
   expressed genes. Furthermore, we discarded cells with more than 30,000
   UMIs and mitochondria content higher than 30%."

Generates QC visualization plots before and after filtering:
  - Violin plots per patient (n_genes, total_counts, pct_counts_mt)
  - Scatter plots (n_genes vs total_counts, colored by pct_counts_mt)
  - Pre/post filtering cell count comparison

Biological interpretation of thresholds:
  min_genes=200: removes empty droplets and completely dead cells
  max_genes=5000: removes likely doublets (two cells captured together)
  max_pct_mt=30%: removes apoptotic/dying cells (high mitochondrial activity)
  min_cells=3:    removes ultra-rare genes (likely noise or alignment errors)
"""

import os
import argparse
import logging
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s [%(levelname)s] %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger(__name__)


def main(input_path, output_path, figures_dir,
         min_genes, max_genes, max_pct_mt, min_cells):

    log.info(f"Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    log.info(f"  Initial: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ── Compute QC metrics ────────────────────────────────────────────────────
    # Identify mitochondrial genes (MT- prefix in human)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    n_mt = adata.var["mt"].sum()
    log.info(f"  Mitochondrial genes: {n_mt}")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars    = ["mt"],
        percent_top = None,
        log1p      = False,
        inplace    = True,
    )

    # ── Pre-filter statistics ─────────────────────────────────────────────────
    log.info("Pre-filter statistics:")
    log.info(f"  Median genes/cell:  {adata.obs['n_genes_by_counts'].median():.0f}")
    log.info(f"  Median UMIs/cell:   {adata.obs['total_counts'].median():.0f}")
    log.info(f"  Median MT%:         {adata.obs['pct_counts_mt'].median():.1f}%")

    # ── Generate pre-filter figures ───────────────────────────────────────────
    os.makedirs(figures_dir, exist_ok=True)
    _plot_qc_overview(adata, figures_dir, suffix="pre_filter",
                      min_genes=min_genes, max_genes=max_genes,
                      max_pct_mt=max_pct_mt)

    # ── Apply filters ─────────────────────────────────────────────────────────
    n_cells_init = adata.n_obs
    n_genes_init = adata.n_vars

    # Cell filters (in order matching paper)
    mask_min  = adata.obs["n_genes_by_counts"] >= min_genes
    mask_max  = adata.obs["n_genes_by_counts"] <= max_genes
    mask_mt   = adata.obs["pct_counts_mt"]     <= max_pct_mt
    mask_all  = mask_min & mask_max & mask_mt

    log.info("Applying cell filters:")
    log.info(f"  min_genes >= {min_genes}:   {(~mask_min).sum():,} cells removed")
    log.info(f"  max_genes <= {max_genes}:   {(~mask_max).sum():,} cells removed")
    log.info(f"  max_pct_mt <= {max_pct_mt}%: {(~mask_mt).sum():,} cells removed")
    log.info(f"  Combined (unique removals): "
             f"{(~mask_all).sum():,} cells removed")

    adata = adata[mask_all].copy()

    # Gene filter
    sc.pp.filter_genes(adata, min_cells=min_cells)
    log.info(f"  min_cells >= {min_cells}: "
             f"{n_genes_init - adata.n_vars:,} genes removed")

    # ── Post-filter statistics ────────────────────────────────────────────────
    log.info("Post-filter statistics:")
    log.info(f"  Cells retained: {adata.n_obs:,} / {n_cells_init:,} "
             f"({adata.n_obs/n_cells_init*100:.1f}%)")
    log.info(f"  Genes retained: {adata.n_vars:,} / {n_genes_init:,}")
    log.info(f"  Patients retained: {adata.obs['patient_id'].nunique()}")
    log.info(f"  Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
    log.info(f"  Avg genes/cell: {adata.obs['n_genes_by_counts'].mean():.0f}")
    log.info(f"  Avg UMIs/cell:  {adata.obs['total_counts'].mean():.0f}")

    # Per-patient cell counts
    log.info("Cells per patient after QC:")
    for pid in sorted(adata.obs["patient_id"].unique(),
                      key=lambda x: int(x.replace("P", ""))):
        n = (adata.obs["patient_id"] == pid).sum()
        log.info(f"  {pid}: {n:,}")

    # ── Post-filter figures ───────────────────────────────────────────────────
    _plot_qc_overview(adata, figures_dir, suffix="post_filter",
                      min_genes=min_genes, max_genes=max_genes,
                      max_pct_mt=max_pct_mt)
    _plot_cell_counts_per_patient(adata, n_cells_init, figures_dir)

    # ── Save ──────────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    log.info(f"Saved to {output_path}")
    log.info("Stage 3 complete.")


def _plot_qc_overview(adata, figures_dir, suffix,
                      min_genes, max_genes, max_pct_mt):
    """Violin + scatter QC plots."""
    patients = sorted(adata.obs["patient_id"].unique(),
                      key=lambda x: int(x.replace("P", "")))
    n_patients = len(patients)

    fig = plt.figure(figsize=(20, 12))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.3)

    # Row 1: Violin plots
    metrics = [
        ("n_genes_by_counts", "Genes per cell", min_genes, max_genes),
        ("total_counts",       "UMIs per cell",  None, None),
        ("pct_counts_mt",      "% MT reads",     None, max_pct_mt),
    ]
    for col, (metric, ylabel, low_thresh, high_thresh) in enumerate(metrics):
        ax = fig.add_subplot(gs[0, col])

        # Subsample to max 30 patients to keep violin readable
        per_patient = [
            adata.obs.loc[adata.obs["patient_id"] == p, metric].values
            for p in patients[:30]
        ]
        ax.violinplot(per_patient, positions=range(len(per_patient)),
                      showmedians=True, showextrema=False)
        ax.set_xlabel("Patient", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(ylabel, fontsize=10, fontweight="bold")
        ax.set_xticks([])

        if low_thresh:
            ax.axhline(low_thresh, color="red", linestyle="--",
                       alpha=0.7, linewidth=1, label=f"min={low_thresh}")
        if high_thresh:
            ax.axhline(high_thresh, color="red", linestyle="--",
                       alpha=0.7, linewidth=1, label=f"max={high_thresh}")
        if low_thresh or high_thresh:
            ax.legend(fontsize=7)

    # Row 2: Scatter plots
    ax4 = fig.add_subplot(gs[1, 0])
    sc_obj = ax4.scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        c=adata.obs["pct_counts_mt"],
        s=0.3, alpha=0.3, cmap="viridis_r", rasterized=True
    )
    plt.colorbar(sc_obj, ax=ax4, label="% MT")
    ax4.set_xlabel("Total UMIs", fontsize=9)
    ax4.set_ylabel("Genes detected", fontsize=9)
    ax4.set_title("UMIs vs Genes (colored by MT%)", fontsize=10,
                  fontweight="bold")
    if high_thresh := max_genes:
        ax4.axhline(high_thresh, color="red", linestyle="--", alpha=0.7)
    if low_thresh := min_genes:
        ax4.axhline(low_thresh, color="orange", linestyle="--", alpha=0.7)

    ax5 = fig.add_subplot(gs[1, 1])
    ax5.hist(adata.obs["n_genes_by_counts"], bins=100,
             color="steelblue", alpha=0.7)
    ax5.axvline(min_genes, color="red", linestyle="--",
                label=f"min={min_genes}")
    ax5.axvline(max_genes, color="red", linestyle="--",
                label=f"max={max_genes}")
    ax5.set_xlabel("Genes per cell", fontsize=9)
    ax5.set_ylabel("Number of cells", fontsize=9)
    ax5.set_title("Gene count distribution", fontsize=10, fontweight="bold")
    ax5.legend(fontsize=8)

    ax6 = fig.add_subplot(gs[1, 2])
    ax6.hist(adata.obs["pct_counts_mt"], bins=100,
             color="tomato", alpha=0.7)
    ax6.axvline(max_pct_mt, color="darkred", linestyle="--",
                label=f"max={max_pct_mt}%")
    ax6.set_xlabel("% Mitochondrial reads", fontsize=9)
    ax6.set_ylabel("Number of cells", fontsize=9)
    ax6.set_title("MT% distribution", fontsize=10, fontweight="bold")
    ax6.legend(fontsize=8)

    fig.suptitle(
        f"QC Overview — {suffix.replace('_', ' ').title()}\n"
        f"{adata.n_obs:,} cells × {adata.n_vars:,} genes",
        fontsize=13, fontweight="bold"
    )
    plt.savefig(os.path.join(figures_dir, f"qc_{suffix}.pdf"),
                dpi=200, bbox_inches="tight")
    plt.savefig(os.path.join(figures_dir, f"qc_{suffix}.png"),
                dpi=200, bbox_inches="tight")
    plt.close()
    log.info(f"  Saved qc_{suffix}.pdf/png")


def _plot_cell_counts_per_patient(adata, n_cells_before, figures_dir):
    """Bar chart of cell counts per patient after QC."""
    patients = sorted(adata.obs["patient_id"].unique(),
                      key=lambda x: int(x.replace("P", "")))
    counts = [
        (adata.obs["patient_id"] == p).sum()
        for p in patients
    ]

    fig, ax = plt.subplots(figsize=(16, 4))
    ax.bar(patients, counts, color="steelblue", alpha=0.8)
    ax.set_xlabel("Patient", fontsize=10)
    ax.set_ylabel("Cells retained", fontsize=10)
    ax.set_title(
        f"Cells per patient after QC — "
        f"{sum(counts):,} total ({sum(counts)/n_cells_before*100:.1f}% retained)",
        fontsize=11, fontweight="bold"
    )
    plt.xticks(rotation=90, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "cell_counts_per_patient.pdf"),
                dpi=200, bbox_inches="tight")
    plt.savefig(os.path.join(figures_dir, "cell_counts_per_patient.png"),
                dpi=200, bbox_inches="tight")
    plt.close()
    log.info("  Saved cell_counts_per_patient.pdf/png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="QC filtering following Wu et al. 2021 thresholds"
    )
    parser.add_argument("--input",       required=True)
    parser.add_argument("--output",      required=True)
    parser.add_argument("--figures_dir", required=True)
    parser.add_argument("--min_genes",   type=int,   default=200)
    parser.add_argument("--max_genes",   type=int,   default=5000)
    parser.add_argument("--max_pct_mt",  type=float, default=30.0)
    parser.add_argument("--min_cells",   type=int,   default=3)
    args = parser.parse_args()
    main(args.input, args.output, args.figures_dir,
         args.min_genes, args.max_genes, args.max_pct_mt, args.min_cells)

"""
Track Comparison — Paper Track vs Alternative Track

Generates a side-by-side report quantifying the biological impact of
scVI batch correction on:
1. Cell type proportions
2. Cancer cell patient specificity
3. Cluster count and resolution
4. UMAP topology (patient mixing)

This comparison justifies the choice of paper track as primary analysis
for studying intertumoral heterogeneity.
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
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec

logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s [%(levelname)s] %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger(__name__)

CELL_TYPE_COLORS = {
    "Cancer":      "#E6851E",
    "Myeloid":     "#E84F9E",
    "Fibroblast":  "#E8D419",
    "T_cell":      "#4A90C4",
    "B_cell":      "#7B5EA7",
    "Neutrophil":  "#8B1A1A",
    "Alveolar":    "#E8401C",
    "Epithelial":  "#9B59B6",
    "Endothelial": "#A8D44B",
    "Mast_cell":   "#1A6B3C",
    "fDC":         "#0D4A3A",
    "Unknown":     "#CCCCCC",
}


def main(paper_path, alt_path, figures_dir, output_path):

    log.info("Loading paper track...")
    paper = sc.read_h5ad(paper_path)
    log.info(f"  Paper track: {paper.n_obs:,} cells")

    log.info("Loading alternative track...")
    alt = sc.read_h5ad(alt_path)
    log.info(f"  Alternative track: {alt.n_obs:,} cells")

    os.makedirs(figures_dir, exist_ok=True)

    # ── Figure 1: Side-by-side UMAP by cell type ──────────────────────────────
    _plot_umap_comparison(paper, alt, figures_dir)

    # ── Figure 2: Cell type proportion comparison ─────────────────────────────
    _plot_proportion_comparison(paper, alt, figures_dir)

    # ── Figure 3: Cancer patient specificity ─────────────────────────────────
    _plot_cancer_patient_specificity(paper, alt, figures_dir)

    # ── Figure 4: Summary table ───────────────────────────────────────────────
    _plot_summary_table(paper, alt, figures_dir)

    # ── Merge all figures into single PDF ─────────────────────────────────────
    _merge_to_report(figures_dir, output_path)

    log.info("Track comparison complete.")


def _plot_umap_comparison(paper, alt, figures_dir):
    """Side-by-side UMAP colored by cell type — paper vs alternative."""
    cell_types_paper = [ct for ct in CELL_TYPE_COLORS
                        if ct in paper.obs["cell_type"].cat.categories]
    cell_types_alt   = [ct for ct in CELL_TYPE_COLORS
                        if ct in alt.obs["cell_type"].cat.categories]

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    titles = ["Paper Track\n(No batch correction, 600 HVGs, 20 PCs)",
              "Alternative Track\n(scVI batch correction, 3000 HVGs)"]

    for ax, adata, cts, title in zip(
            axes, [paper, alt],
            [cell_types_paper, cell_types_alt], titles):
        for ct in cts:
            mask = adata.obs["cell_type"] == ct
            ax.scatter(
                adata.obsm["X_umap"][mask, 0],
                adata.obsm["X_umap"][mask, 1],
                c=CELL_TYPE_COLORS[ct], s=0.3, alpha=0.7, rasterized=True
            )
        ax.set_xlabel("UMAP1", fontsize=10)
        ax.set_ylabel("UMAP2", fontsize=10)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_aspect("equal")

    handles = [Patch(facecolor=CELL_TYPE_COLORS[ct], label=ct)
               for ct in cell_types_paper]
    axes[1].legend(handles=handles, loc="upper left",
                   bbox_to_anchor=(1.01, 1), fontsize=8,
                   title="Cell type", title_fontsize=9)
    plt.suptitle("Track Comparison — UMAP by cell type", fontsize=13,
                 fontweight="bold")
    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir, f"compare_umap_celltype.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log.info("  Saved compare_umap_celltype")


def _plot_proportion_comparison(paper, alt, figures_dir):
    """Bar chart comparing overall cell type proportions."""
    # Get all cell types present in either track
    all_ct = [ct for ct in CELL_TYPE_COLORS
              if ct in paper.obs["cell_type"].cat.categories
              or ct in alt.obs["cell_type"].cat.categories
              if ct != "Unknown"]

    paper_pct = paper.obs["cell_type"].value_counts(normalize=True) * 100
    alt_pct   = alt.obs["cell_type"].value_counts(normalize=True) * 100

    paper_vals = [paper_pct.get(ct, 0) for ct in all_ct]
    alt_vals   = [alt_pct.get(ct, 0)   for ct in all_ct]

    x = np.arange(len(all_ct))
    width = 0.35

    fig, ax = plt.subplots(figsize=(14, 5))
    bars1 = ax.bar(x - width/2, paper_vals, width,
                   color=[CELL_TYPE_COLORS[ct] for ct in all_ct],
                   alpha=0.9, label="Paper track", edgecolor="k",
                   linewidth=0.5)
    bars2 = ax.bar(x + width/2, alt_vals, width,
                   color=[CELL_TYPE_COLORS[ct] for ct in all_ct],
                   alpha=0.5, label="Alternative track", edgecolor="k",
                   linewidth=0.5, hatch="//")

    ax.set_xticks(x)
    ax.set_xticklabels(all_ct, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Percentage of cells (%)", fontsize=11)
    ax.set_title("Cell type proportions: Paper track vs Alternative track",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir,
                                 f"compare_proportions.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log.info("  Saved compare_proportions")


def _plot_cancer_patient_specificity(paper, alt, figures_dir):
    """
    Quantify how much cancer cells cluster by patient.

    High patient specificity (paper track) = large variance in cancer
    cell UMAP coordinates per patient = preserved heterogeneity.
    Low patient specificity (alt track) = harmonized, heterogeneity lost.
    """
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    titles = ["Paper Track — Cancer cells by patient",
              "Alternative Track — Cancer cells by patient"]

    patients = sorted(
        paper.obs["patient_id"].unique(),
        key=lambda x: int(x.replace("P", ""))
    )
    cmap1 = plt.cm.get_cmap("tab20", 20)
    cmap2 = plt.cm.get_cmap("tab20b", 20)
    cmap3 = plt.cm.get_cmap("tab20c", 20)
    pool = ([cmap1(i) for i in range(20)] +
            [cmap2(i) for i in range(20)] +
            [cmap3(i) for i in range(4)])
    pcolors = {p: pool[i] for i, p in enumerate(patients)}

    for ax, adata, title in zip(axes, [paper, alt], titles):
        cancer_mask = adata.obs["cell_type"] == "Cancer"
        cancer_obs  = adata.obs[cancer_mask]
        cancer_umap = adata.obsm["X_umap"][cancer_mask]

        for p in patients:
            pmask = cancer_obs["patient_id"] == p
            if pmask.sum() == 0:
                continue
            ax.scatter(
                cancer_umap[pmask.values, 0],
                cancer_umap[pmask.values, 1],
                c=[pcolors[p]], s=0.5, alpha=0.7, rasterized=True
            )

        n_cancer = cancer_mask.sum()
        ax.set_title(f"{title}\nn={n_cancer:,} cancer cells",
                     fontsize=11, fontweight="bold")
        ax.set_xlabel("UMAP1", fontsize=10)
        ax.set_ylabel("UMAP2", fontsize=10)
        ax.set_aspect("equal")

    plt.suptitle("Cancer cell patient specificity — key biological validation",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir,
                                 f"compare_cancer_patient.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log.info("  Saved compare_cancer_patient")


def _plot_summary_table(paper, alt, figures_dir):
    """Summary table of key metrics across both tracks."""

    paper_ann = paper.uns.get("annotation", {})
    alt_ann   = alt.uns.get("annotation", {})

    paper_ct = paper.obs["cell_type"].value_counts(normalize=True) * 100
    alt_ct   = alt.obs["cell_type"].value_counts(normalize=True) * 100

    rows = []
    rows.append(["Metric", "Paper Track", "Alternative Track", "Wu et al. 2021"])
    rows.append(["─" * 20] * 4)
    rows.append(["Cells",
                 f"{paper.n_obs:,}",
                 f"{alt.n_obs:,}",
                 "90,406"])
    rows.append(["HVGs",
                 "600",
                 "3000",
                 "600"])
    rows.append(["Batch correction",
                 "None",
                 "scVI",
                 "None"])
    rows.append(["PCA / Latent dims",
                 "20 PCs",
                 "30 scVI dims",
                 "20 PCs"])
    rows.append(["Leiden resolution",
                 str(paper_ann.get("resolution", "?")),
                 str(alt_ann.get("resolution", "?")),
                 "1.0"])
    rows.append(["N clusters",
                 str(paper_ann.get("n_clusters", "?")),
                 str(alt_ann.get("n_clusters", "?")),
                 "37"])
    rows.append(["Silhouette",
                 str(paper_ann.get("silhouette", "?")),
                 str(alt_ann.get("silhouette", "?")),
                 "N/A"])
    rows.append(["─" * 20] * 4)
    rows.append(["Cancer %",
                 f"{paper_ct.get('Cancer', 0):.1f}%",
                 f"{alt_ct.get('Cancer', 0):.1f}%",
                 "~35-40%"])
    rows.append(["T_cell %",
                 f"{paper_ct.get('T_cell', 0):.1f}%",
                 f"{alt_ct.get('T_cell', 0):.1f}%",
                 "~18-20%"])
    rows.append(["Myeloid %",
                 f"{paper_ct.get('Myeloid', 0):.1f}%",
                 f"{alt_ct.get('Myeloid', 0):.1f}%",
                 "~15%"])
    rows.append(["Neutrophil %",
                 f"{paper_ct.get('Neutrophil', 0):.1f}%",
                 f"{alt_ct.get('Neutrophil', 0):.1f}%",
                 "~5-7%"])
    rows.append(["fDC %",
                 f"{paper_ct.get('fDC', 0):.2f}%",
                 f"{alt_ct.get('fDC', 0):.2f}%",
                 "<2%"])

    fig, ax = plt.subplots(figsize=(12, len(rows) * 0.5 + 1))
    ax.axis("off")
    tbl = ax.table(
        cellText  = rows,
        cellLoc   = "left",
        loc       = "center",
        bbox      = [0, 0, 1, 1],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_facecolor("#2C3E50")
            cell.set_text_props(color="white", fontweight="bold")
        elif rows[row][0].startswith("─"):
            cell.set_facecolor("#ECF0F1")
        elif row % 2 == 0:
            cell.set_facecolor("#FDFEFE")
        cell.set_edgecolor("none")

    plt.title("Pipeline Track Comparison Summary",
              fontsize=13, fontweight="bold", pad=15)
    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir, f"compare_summary_table.{ext}"),
                    dpi=200, bbox_inches="tight")
    plt.close()
    log.info("  Saved compare_summary_table")


def _merge_to_report(figures_dir, output_path):
    """Combine all comparison figures into a single PDF report."""
    try:
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.image as mpimg

        figures = [
            "compare_umap_celltype.png",
            "compare_proportions.png",
            "compare_cancer_patient.png",
            "compare_summary_table.png",
        ]

        with PdfPages(output_path) as pdf:
            for fig_name in figures:
                fig_path = os.path.join(figures_dir, fig_name)
                if not os.path.exists(fig_path):
                    continue
                img = mpimg.imread(fig_path)
                fig, ax = plt.subplots(figsize=(16, 9))
                ax.imshow(img)
                ax.axis("off")
                pdf.savefig(fig, bbox_inches="tight")
                plt.close()

        log.info(f"  Report saved to {output_path}")

    except Exception as e:
        log.warning(f"  Could not merge PDF report: {e}")
        log.info("  Individual figures are available in " + figures_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare paper and alternative analysis tracks"
    )
    parser.add_argument("--paper_h5ad",  required=True)
    parser.add_argument("--alt_h5ad",    required=True)
    parser.add_argument("--figures_dir", required=True)
    parser.add_argument("--output",      required=True)
    args = parser.parse_args()
    main(args.paper_h5ad, args.alt_h5ad, args.figures_dir, args.output)

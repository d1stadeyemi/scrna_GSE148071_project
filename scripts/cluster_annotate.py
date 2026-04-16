"""
Stage 5 — Leiden Clustering and Cell Type Annotation

Input:  data/processed/embedded.h5ad
Output: data/processed/annotated.h5ad
        results/figures/annotation/
        results/tables/cluster_markers.csv
        results/tables/cell_type_assignments.csv

Steps:
- Leiden clustering at multiple resolutions
- Marker gene scoring per cluster
- Cell type annotation using paper's canonical markers
- Reproduce Fig. 1b, 1c, 1d, 1e.
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
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import seaborn as sns

logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s [%(levelname)s] %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger(__name__)

# ── Cell type colors matching publication ─────────────────────────────────────
CELL_TYPE_COLORS = {
    "Cancer":      "#E6851E",   # orange
    "Myeloid":     "#E84F9E",   # magenta/pink
    "Fibroblast":  "#E8D419",   # yellow
    "T_cell":      "#4A90C4",   # steel blue
    "B_cell":      "#7B5EA7",   # purple
    "Neutrophil":  "#8B1A1A",   # dark red
    "Alveolar":    "#E8401C",   # red-orange
    "Epithelial":  "#9B59B6",   # medium purple
    "Endothelial": "#A8D44B",   # yellow-green
    "Mast_cell":   "#1A6B3C",   # dark green
    "fDC":         "#0D4A3A",   # very dark green
    "Unknown":     "#CCCCCC",   # grey
}

# ── Canonical marker genes from the paper (Supplementary Table 2) ─────────────
MARKER_GENES = {
    "Endothelial": ["CLDN5", "VWF", "PECAM1", "FLT1", "CDH5"],
    "Epithelial":  ["CAPS", "SNTN", "TPPP3", "FOXJ1"],
    "Alveolar":    ["CLDN18", "AQP4", "SFTPA1", "SFTPC", "ABCA3"],
    "Fibroblast":  ["COL1A1", "COL1A2", "DCN", "LUM", "FAP"],
    "T_cell":      ["CD2", "CD3D", "CD3E", "CD3G", "TRAC"],
    "B_cell":      ["CD79A", "CD79B", "MS4A1", "IGLC3", "MZB1"],
    "Myeloid":     ["CD14", "LYZ", "CD68", "FCGR3A", "MNDA"],
    "Neutrophil":  ["CSF3R", "S100A8", "S100A9", "FCGR3B", "CXCR2"],
    "fDC":         ["FDCSP", "CR2", "CLU"],
    "Mast_cell":   ["GATA2", "TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "Cancer":      ["EPCAM", "KRT18", "KRT19", "MKI67", "TP63"],
}

# Flat list of all markers for heatmap
ALL_MARKERS = [g for genes in MARKER_GENES.values() for g in genes]


def run_leiden(adata, resolution):
    """Run Leiden clustering at given resolution."""
    log.info(f"Running Leiden clustering (resolution={resolution})...")
    sc.tl.leiden(adata, resolution=resolution, key_added=f"leiden_{resolution}")
    n_clusters = adata.obs[f"leiden_{resolution}"].nunique()
    log.info(f"  Found {n_clusters} clusters at resolution {resolution}")
    return n_clusters


def score_and_annotate(adata, resolution):
    """Score each cluster by marker genes and assign cell type."""
    cluster_key = f"leiden_{resolution}"
    log.info(f"Scoring clusters for cell type annotation...")

    # Score each cell type using scanpy gene scoring
    for cell_type, markers in MARKER_GENES.items():
        valid_markers = [g for g in markers if g in adata.var_names]
        if len(valid_markers) == 0:
            log.warning(f"  No markers found for {cell_type}")
            continue
        sc.tl.score_genes(
            adata,
            gene_list  = valid_markers,
            score_name = f"score_{cell_type}",
        )
        log.info(f"  {cell_type}: {len(valid_markers)}/{len(markers)} markers found")

    # Assign cell type to each cluster based on highest mean score
    score_cols = [f"score_{ct}" for ct in MARKER_GENES.keys()]
    cluster_scores = adata.obs.groupby(cluster_key)[score_cols].mean()
    cluster_scores.columns = list(MARKER_GENES.keys())

    cluster_annotation = cluster_scores.idxmax(axis=1)
    log.info("Cluster annotations:")
    for cluster, cell_type in cluster_annotation.items():
        n_cells = (adata.obs[cluster_key] == cluster).sum()
        log.info(f"  Cluster {cluster:>3} → {cell_type:<15} ({n_cells:>6} cells)")

    # Map back to cells
    adata.obs["cell_type"] = adata.obs[cluster_key].map(cluster_annotation)
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

    return cluster_annotation, cluster_scores


def plot_umap_cell_type(adata, figures_dir):
    """Reproduce Fig. 1b — UMAP colored by cell type."""
    log.info("Plotting UMAP by cell type (Fig. 1b)...")

    # Build ordered palette matching paper
    cell_types_present = [ct for ct in CELL_TYPE_COLORS.keys()
                          if ct in adata.obs["cell_type"].cat.categories]
    palette = {ct: CELL_TYPE_COLORS[ct] for ct in cell_types_present}

    fig, ax = plt.subplots(figsize=(8, 7))

    for ct in cell_types_present:
        mask = adata.obs["cell_type"] == ct
        ax.scatter(
            adata.obsm["X_umap"][mask, 0],
            adata.obsm["X_umap"][mask, 1],
            c     = palette[ct],
            s     = 0.5,
            alpha = 0.7,
            rasterized = True,
            label = ct,
        )

    # Legend
    legend_handles = [
        Patch(facecolor=palette[ct], label=ct)
        for ct in cell_types_present
    ]
    ax.legend(
        handles       = legend_handles,
        loc           = "upper left",
        bbox_to_anchor = (1.01, 1),
        fontsize      = 9,
        frameon       = True,
        title         = "Cell type",
        title_fontsize = 10,
    )
    ax.set_xlabel("UMAP1", fontsize=11)
    ax.set_ylabel("UMAP2", fontsize=11)
    ax.set_title(f"UMAP — {adata.n_obs:,} cells, {len(cell_types_present)} cell types",
                 fontsize=12, fontweight="bold")
    ax.set_aspect("equal")

    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "umap_cell_type.pdf"),
                dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(figures_dir, "umap_cell_type.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    log.info("  Saved umap_cell_type.pdf/png")


def plot_umap_patient(adata, figures_dir):
    """Reproduce Fig. 1d — UMAP colored by patient."""
    log.info("Plotting UMAP by patient (Fig. 1d)...")

    patients = sorted(adata.obs["patient_id"].unique(),
                      key=lambda x: int(x.replace("P", "")))
    n_patients = len(patients)

    # Generate 42 visually distinct colors
    cmap = plt.cm.get_cmap("tab20", 20)
    cmap2 = plt.cm.get_cmap("tab20b", 20)
    cmap3 = plt.cm.get_cmap("tab20c", 20)
    colors_pool = (
        [cmap(i) for i in range(20)] +
        [cmap2(i) for i in range(20)] +
        [cmap3(i) for i in range(4)]
    )
    patient_colors = {p: colors_pool[i] for i, p in enumerate(patients)}

    fig, ax = plt.subplots(figsize=(10, 8))

    for patient in patients:
        mask = adata.obs["patient_id"] == patient
        ax.scatter(
            adata.obsm["X_umap"][mask, 0],
            adata.obsm["X_umap"][mask, 1],
            c          = [patient_colors[patient]],
            s          = 0.3,
            alpha      = 0.6,
            rasterized = True,
        )

    # Legend in 3 columns
    legend_handles = [
        Patch(facecolor=patient_colors[p], label=p)
        for p in patients
    ]
    ax.legend(
        handles        = legend_handles,
        loc            = "upper left",
        bbox_to_anchor = (1.01, 1),
        fontsize       = 7,
        ncol           = 2,
        frameon        = True,
        title          = "Patient",
        title_fontsize = 9,
        markerscale    = 1.5,
    )
    ax.set_xlabel("UMAP1", fontsize=11)
    ax.set_ylabel("UMAP2", fontsize=11)
    ax.set_title("UMAP — colored by patient", fontsize=12, fontweight="bold")
    ax.set_aspect("equal")

    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "umap_by_patient_annotated.pdf"),
                dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(figures_dir, "umap_by_patient_annotated.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    log.info("  Saved umap_by_patient_annotated.pdf/png")


def plot_marker_heatmap(adata, figures_dir):
    """Reproduce Fig. 1c — marker gene heatmap."""
    log.info("Plotting marker gene heatmap (Fig. 1c)...")

    # Subsample for heatmap speed — 200 cells per cell type max
    cell_types = [ct for ct in CELL_TYPE_COLORS.keys()
                  if ct in adata.obs["cell_type"].cat.categories
                  and ct != "Unknown"]

    sampled_indices = []
    for ct in cell_types:
        idx = np.where(adata.obs["cell_type"] == ct)[0]
        n_sample = min(200, len(idx))
        sampled_indices.extend(np.random.choice(idx, n_sample, replace=False))

    sampled_indices = sorted(sampled_indices)
    adata_sub = adata[sampled_indices].copy()

    # Get log1p normalized expression for heatmap
    if "log1p_norm" in adata_sub.layers:
        expr = adata_sub.layers["log1p_norm"]
    else:
        expr = adata_sub.X

    # Build expression dataframe for markers
    valid_markers = [g for g in ALL_MARKERS if g in adata_sub.var_names]
    marker_idx = [list(adata_sub.var_names).index(g) for g in valid_markers]

    if hasattr(expr, "toarray"):
        expr_dense = expr[:, marker_idx].toarray()
    else:
        expr_dense = np.array(expr[:, marker_idx])

    df_expr = pd.DataFrame(
        expr_dense,
        index   = adata_sub.obs.index,
        columns = valid_markers,
    )
    df_expr["cell_type"] = adata_sub.obs["cell_type"].values

    # Sort by cell type order
    df_expr["cell_type"] = pd.Categorical(
        df_expr["cell_type"], categories=cell_types, ordered=True
    )
    df_expr = df_expr.sort_values("cell_type")

    # Row colors
    row_colors = df_expr["cell_type"].map(CELL_TYPE_COLORS)
    expr_matrix = df_expr.drop(columns=["cell_type"]).T

    # Clustermap
    g = sns.clustermap(
        expr_matrix,
        col_colors      = row_colors.values,
        col_cluster     = False,
        row_cluster     = False,
        cmap            = "RdBu_r",
        vmin            = 0,
        vmax            = 3,
        figsize         = (14, 10),
        xticklabels     = False,
        yticklabels     = True,
        linewidths      = 0,
        rasterized      = True,
        cbar_kws        = {"label": "log1p expression", "shrink": 0.3},
    )

    g.ax_heatmap.set_ylabel("Marker genes", fontsize=10)
    g.ax_heatmap.set_xlabel("Cells", fontsize=10)
    g.ax_heatmap.yaxis.set_tick_params(labelsize=7)
    g.fig.suptitle("Canonical cell type markers", fontsize=13,
                   fontweight="bold", y=1.01)

    # Cell type color legend
    legend_handles = [
        Patch(facecolor=CELL_TYPE_COLORS[ct], label=ct)
        for ct in cell_types
    ]
    g.ax_heatmap.legend(
        handles        = legend_handles,
        loc            = "upper left",
        bbox_to_anchor = (1.15, 1),
        fontsize       = 8,
        title          = "Cell type",
        title_fontsize = 9,
        frameon        = True,
    )

    plt.savefig(os.path.join(figures_dir, "marker_heatmap.pdf"),
                dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(figures_dir, "marker_heatmap.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    log.info("  Saved marker_heatmap.pdf/png")


def plot_cell_type_proportions(adata, figures_dir):
    """Reproduce Fig. 1e — stacked bar chart of cell type proportions per patient."""
    log.info("Plotting cell type proportions per patient (Fig. 1e)...")

    patients = sorted(adata.obs["patient_id"].unique(),
                      key=lambda x: int(x.replace("P", "")))
    cell_types = [ct for ct in CELL_TYPE_COLORS.keys()
                  if ct in adata.obs["cell_type"].cat.categories
                  and ct != "Unknown"]

    # Compute proportions
    prop_df = pd.crosstab(
        adata.obs["patient_id"],
        adata.obs["cell_type"],
        normalize = "index",
    )
    prop_df = prop_df.reindex(patients)
    prop_df = prop_df.reindex(columns=cell_types, fill_value=0)

    fig, ax = plt.subplots(figsize=(18, 5))

    bottom = np.zeros(len(patients))
    for ct in cell_types:
        values = prop_df[ct].values
        ax.bar(
            range(len(patients)),
            values,
            bottom     = bottom,
            color      = CELL_TYPE_COLORS[ct],
            label      = ct,
            width      = 0.85,
            edgecolor  = "none",
        )
        bottom += values

    ax.set_xticks(range(len(patients)))
    ax.set_xticklabels(patients, rotation=90, fontsize=7)
    ax.set_xlabel("Patient", fontsize=11)
    ax.set_ylabel("Proportion", fontsize=11)
    ax.set_ylim(0, 1)
    ax.set_title("Cell type composition per patient", fontsize=12, fontweight="bold")

    handles = [Patch(facecolor=CELL_TYPE_COLORS[ct], label=ct) for ct in cell_types]
    ax.legend(
        handles        = handles,
        loc            = "upper left",
        bbox_to_anchor = (1.01, 1),
        fontsize       = 8,
        title          = "Cell type",
        title_fontsize = 9,
        frameon        = True,
    )

    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "cell_type_proportions.pdf"),
                dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(figures_dir, "cell_type_proportions.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    log.info("  Saved cell_type_proportions.pdf/png")


def save_tables(adata, cluster_annotation, tables_dir):
    """Save cluster marker and cell type assignment tables."""
    os.makedirs(tables_dir, exist_ok=True)

    # Cell type assignments
    assignments = adata.obs[["patient_id", "cell_type"]].copy()
    assignments.to_csv(os.path.join(tables_dir, "cell_type_assignments.csv"))
    log.info(f"  Saved cell_type_assignments.csv")

    # Cell type counts
    counts = adata.obs["cell_type"].value_counts().reset_index()
    counts.columns = ["cell_type", "n_cells"]
    counts["pct"] = (counts["n_cells"] / adata.n_obs * 100).round(2)
    counts.to_csv(os.path.join(tables_dir, "cell_type_counts.csv"), index=False)
    log.info(f"  Saved cell_type_counts.csv")
    log.info("Cell type summary:")
    for _, row in counts.iterrows():
        log.info(f"  {row['cell_type']:<15} {row['n_cells']:>7,} cells ({row['pct']:.1f}%)")


def main(input_path, output_path, figures_dir, tables_dir, resolution):

    log.info(f"Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    log.info(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)

    # ── Leiden clustering ─────────────────────────────────────────────────────
    run_leiden(adata, resolution)

    # ── Cell type annotation ──────────────────────────────────────────────────
    cluster_annotation, cluster_scores = score_and_annotate(adata, resolution)

    # ── Figures ───────────────────────────────────────────────────────────────
    plot_umap_cell_type(adata, figures_dir)
    plot_umap_patient(adata, figures_dir)
    plot_marker_heatmap(adata, figures_dir)
    plot_cell_type_proportions(adata, figures_dir)

    # ── Tables ────────────────────────────────────────────────────────────────
    save_tables(adata, cluster_annotation, tables_dir)

    # ── Save annotated AnnData ────────────────────────────────────────────────
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    log.info(f"Saved annotated AnnData to {output_path}")
    log.info("Stage 5 complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",       default="data/processed/embedded.h5ad")
    parser.add_argument("--output",      default="data/processed/annotated.h5ad")
    parser.add_argument("--figures_dir", default="results/figures/annotation")
    parser.add_argument("--tables_dir",  default="results/tables")
    parser.add_argument("--resolution",  type=float, default=1.0)
    args = parser.parse_args()
    main(args.input, args.output, args.figures_dir,
         args.tables_dir, args.resolution)
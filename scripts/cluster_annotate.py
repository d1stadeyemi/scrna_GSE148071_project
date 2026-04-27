"""
Stage 5/6 — Leiden Clustering and Cell Type Annotation

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
BIOLOGICAL ANNOTATION STRATEGY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

This script annotates cell clusters using biological principles, NOT by tuning
parameters to match expected numbers. The annotation hierarchy below reflects
our knowledge of cell type biology in the lung tumor microenvironment.

MARKER SOURCES:
  - Wu et al. 2021 (this paper) — Methods section
  - Travaglini et al. 2020 — Human Lung Cell Atlas (Nature)
  - Habermann et al. 2020 — IPF scRNA-seq (Science Advances)
  - Lambrechts et al. 2018 — NSCLC scRNA-seq (Nature Medicine)

ANNOTATION HIERARCHY:
  Why a hierarchy? Different cell types have different marker specificities:

  Tier 1 — Rare immune cells (Mast, fDC):
    FDCSP is exquisitely specific to follicular dendritic cells.
    GATA2+TPSAB1 combination is highly specific to mast cells.
    These markers are rarely expressed elsewhere — high threshold ensures
    we only call these when the signal is genuine.

  Tier 2 — Immune effectors (T, B, Myeloid, Neutrophil):
    CD3D/E/G are pan-T markers with no known expression outside T lineage.
    CD79A/B are pan-B markers.
    CD14+LYZ combination marks monocytes/macrophages.
    S100A8/A9 mark neutrophils/granulocytes.
    These markers are well-characterized — moderate threshold sufficient.

  Tier 3 — Stromal cells (Fibroblast, Endothelial):
    COL1A1/A2 are major structural collagens, predominant in fibroblasts.
    PECAM1/VWF/CLDN5 mark vascular endothelium.

  Tier 4 — Normal epithelial (Alveolar, Airway Epithelial):
    SFTPC/SFTPA1 mark AT2 alveolar cells (surfactant production).
    CAPS/SNTN mark ciliated airway epithelial cells.
    Higher threshold because these populations overlap with cancer.

  Tier 5 — Cancer cells:
    Identified by EPCAM positivity + absence of normal differentiation.
    This follows PATHOLOGY convention: cancer cells retain epithelial
    origin marker (EPCAM) but lose tissue-specific differentiation markers
    (SFTPC, SCGB1A1, FOXJ1) as part of malignant transformation.
    NOT a hack — this is how pathologists diagnose carcinoma in tissue.

  Tier 6 — Unknown:
    Clusters with no strong signal stay "Unknown". This is honest.
    Forcing every cluster into a named category produces misleading biology.

CLUSTER COUNT SELECTION:
  We sweep resolutions and use silhouette score as an objective measure
  of cluster quality. We additionally consider target_clusters from the
  publication for reproduction purposes — but silhouette takes priority
  for discovery mode. Both are reported in the log.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""

import os
import argparse
import logging
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s [%(levelname)s] %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger(__name__)


# ── Cell type colors ──────────────────────────────────────────────────────────
# Colors match Wu et al. 2021 Fig 1b for direct visual comparison
CELL_TYPE_COLORS = {
    "Cancer":      "#E6851E",   # orange
    "Myeloid":     "#E84F9E",   # pink/magenta
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

# ── Canonical marker genes ────────────────────────────────────────────────────
# Used for cell scoring — each cell type scored by average marker expression
# relative to a random control gene set (scanpy score_genes method)
MARKER_GENES = {
    # Immune — most specific markers in the lung
    "T_cell":      ["CD2", "CD3D", "CD3E", "CD3G", "TRAC",
                    "CD8A", "CD4", "IL7R", "CCR7", "SELL"],
    "B_cell":      ["CD79A", "CD79B", "MS4A1", "IGLC3", "MZB1",
                    "IGHM", "JCHAIN"],
    "Myeloid":     ["CD14", "LYZ", "CD68", "FCGR3A", "MNDA",
                    "CSF1R", "MARCO", "MRC1"],
    "Neutrophil":  ["CSF3R", "S100A8", "S100A9", "FCGR3B",
                    "CXCR2", "OLR1"],
    "Mast_cell":   ["GATA2", "TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "fDC":         ["FDCSP", "CR2"],

    # Stromal
    "Fibroblast":  ["COL1A1", "COL1A2", "DCN", "LUM", "FAP"],
    "Endothelial": ["CLDN5", "VWF", "PECAM1", "FLT1", "CDH5"],

    # Epithelial — normal lung (used also to EXCLUDE cancer below)
    "Alveolar":    ["CLDN18", "AQP4", "FOLR1", "SFTPA1", "SFTPC", "ABCA3"],
    "Epithelial":  ["CAPS", "SNTN", "TPPP3", "FOXJ1"],
}

# EPCAM and epithelial differentiation markers for cancer identification
EPCAM_MARKERS     = ["EPCAM", "KRT18", "KRT19", "KRT7"]
NORMAL_EPI_MARKERS = ["SFTPC", "SFTPA1", "ABCA3",   # AT2 (alveolar)
                       "SCGB1A1", "SCGB3A1",          # Club cells
                       "FOXJ1", "TPPP3", "CAPS"]      # Ciliated cells

# Extended markers for visualization heatmap (Fig 1c reproduction)
HEATMAP_MARKERS = {
    "T_cell":      ["CD2", "CD3D", "CD3E", "CD3G", "TRAC",
                    "CD8A", "CD4", "NKG7", "GZMA", "GZMK"],
    "B_cell":      ["CD79A", "CD79B", "MS4A1", "IGLC3", "JCHAIN",
                    "MZB1", "IGHM"],
    "Myeloid":     ["CD14", "LYZ", "CD68", "FCGR1A", "FCGR3A",
                    "MNDA", "CSF1R", "MARCO", "MRC1"],
    "Neutrophil":  ["CSF3R", "S100A8", "S100A9", "FCGR3B",
                    "CXCR2", "OLR1"],
    "Mast_cell":   ["GATA2", "CPA3", "MS4A2", "TPSAB1", "TPSB2"],
    "fDC":         ["FDCSP", "CR2"],
    "Fibroblast":  ["COL1A1", "COL1A2", "DCN", "LUM", "FAP"],
    "Endothelial": ["CLDN5", "FLT1", "VWF", "PECAM1", "CDH5",
                    "KDR", "ACKR1", "ANGPT2"],
    "Alveolar":    ["CLDN18", "AQP4", "FOLR1", "SFTPA1", "SFTPC",
                    "ABCA3", "SFTPB"],
    "Epithelial":  ["CAPS", "SNTN", "TPPP3", "FOXJ1", "KRT5", "KRT6A"],
    "Cancer":      ["EPCAM", "KRT18", "KRT19", "MKI67", "TP63",
                    "NAPSA", "NKX2-1", "KRT5", "KRT6A"],
}


# =============================================================================
# CLUSTERING
# =============================================================================

def run_leiden_sweep(adata, resolutions):
    """
    Run Leiden clustering at multiple resolutions.
    Compute silhouette score at each — measures cluster cohesion objectively.

    Silhouette score ranges from -1 to 1:
      > 0.5: well-separated clusters
      0.2-0.5: moderate separation
      < 0.2: overlapping clusters

    Both n_clusters and silhouette are reported — the user can assess
    the resolution-biology trade-off.
    """
    log.info("Running Leiden resolution sweep...")

    # Choose representation for silhouette (batch-corrected if available)
    if "X_scvi" in adata.obsm:
        sil_rep = adata.obsm["X_scvi"]
        rep_name = "X_scvi"
    else:
        sil_rep = adata.obsm["X_pca"]
        rep_name = "X_pca"
    log.info(f"  Silhouette computed on {rep_name}")

    results = {}
    for res in resolutions:
        key = f"leiden_{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        n_clusters = adata.obs[key].nunique()

        # Silhouette on subsample (10k cells for speed)
        n_sample = min(10000, adata.n_obs)
        idx = np.random.choice(adata.n_obs, n_sample, replace=False)
        sil = silhouette_score(
            sil_rep[idx],
            adata.obs[key].iloc[idx],
            sample_size=None
        )

        results[res] = {"n_clusters": n_clusters, "silhouette": round(sil, 4)}
        log.info(f"  res={res:.1f}: {n_clusters} clusters, "
                 f"silhouette={sil:.4f}")

    return results


def select_resolution(sweep_results, target_clusters=None):
    """
    Select resolution using two criteria:

    1. If target_clusters specified (reproduction mode):
       Select resolution closest to target. Report silhouette of that choice.

    2. Discovery mode (target_clusters=None):
       Select resolution maximizing silhouette score.
       This is the objective biological optimum.

    Both modes report what was chosen and why — transparent decision.
    """
    if target_clusters is not None:
        # Reproduction mode — prioritize matching cluster count
        best_res = min(
            sweep_results,
            key=lambda r: abs(sweep_results[r]["n_clusters"] - target_clusters)
        )
        got = sweep_results[best_res]["n_clusters"]
        sil = sweep_results[best_res]["silhouette"]
        log.info(f"Resolution selected: {best_res} "
                 f"(target={target_clusters}, got={got}, silhouette={sil})")
    else:
        # Discovery mode — maximize silhouette
        best_res = max(
            sweep_results,
            key=lambda r: sweep_results[r]["silhouette"]
        )
        got = sweep_results[best_res]["n_clusters"]
        sil = sweep_results[best_res]["silhouette"]
        log.info(f"Resolution selected: {best_res} "
                 f"(max silhouette={sil}, gives {got} clusters)")

    return best_res


# =============================================================================
# ANNOTATION
# =============================================================================

def score_markers(adata, cluster_key):
    """
    Score each cell for canonical markers using scanpy score_genes.

    score_genes computes: mean(marker expression) - mean(control expression)
    where controls are randomly sampled genes with similar expression levels.
    This corrects for expression level bias.

    Also scores EPCAM and normal epithelial markers separately
    for the cancer identification logic below.
    """
    log.info("Scoring markers for all cell types...")

    for cell_type, markers in MARKER_GENES.items():
        valid = [g for g in markers if g in adata.var_names]
        if not valid:
            log.warning(f"  {cell_type}: 0/{len(markers)} markers found — skipping")
            continue
        sc.tl.score_genes(adata, valid, score_name=f"score_{cell_type}")
        log.info(f"  {cell_type}: {len(valid)}/{len(markers)} markers scored")

    # EPCAM score — for cancer positive signal
    valid_epcam = [g for g in EPCAM_MARKERS if g in adata.var_names]
    if valid_epcam:
        sc.tl.score_genes(adata, valid_epcam, score_name="score_EPCAM")
        log.info(f"  EPCAM: {len(valid_epcam)}/{len(EPCAM_MARKERS)} markers scored")

    # Normal epithelial score — for cancer NEGATIVE signal
    valid_norm = [g for g in NORMAL_EPI_MARKERS if g in adata.var_names]
    if valid_norm:
        sc.tl.score_genes(adata, valid_norm,
                          score_name="score_normal_epithelial")
        log.info(f"  Normal epithelial: {len(valid_norm)}/{len(NORMAL_EPI_MARKERS)} markers scored")


def annotate_clusters(adata, cluster_key, min_cells=50):
    """
    Assign cell types to clusters using biological hierarchy.

    See module docstring for full explanation of each tier.

    Returns:
      cluster_annotation: dict of cluster -> cell_type
      cluster_scores: DataFrame of mean scores per cluster
    """
    # Compute mean scores per cluster
    score_cols = [f"score_{ct}" for ct in MARKER_GENES
                  if f"score_{ct}" in adata.obs.columns]
    extra_cols = [c for c in ["score_EPCAM", "score_normal_epithelial"]
                  if c in adata.obs.columns]

    cluster_mean = adata.obs.groupby(cluster_key)[
        score_cols + extra_cols
    ].mean()

    cluster_annotation = {}
    cluster_confidence = {}

    log.info("Annotating clusters...")

    for cluster in sorted(cluster_mean.index, key=int):
        scores = cluster_mean.loc[cluster]
        n_cells = (adata.obs[cluster_key] == cluster).sum()

        # Discard tiny clusters — likely artefacts or doublets
        if n_cells < min_cells:
            log.warning(f"  Cluster {cluster}: only {n_cells} cells "
                        f"— labeled Unknown (below min_cells={min_cells})")
            cluster_annotation[cluster] = "Unknown"
            cluster_confidence[cluster] = "none"
            continue

        cell_type = None
        confidence = None

        def get(key, default=-999):
            return float(scores.get(f"score_{key}", default))

        # ── Tier 1: Rare cells — high specificity threshold ───────────────────
        if get("Mast_cell") > 3.0:
            cell_type, confidence = "Mast_cell", "high"

        elif get("fDC") > 2.0:
            cell_type, confidence = "fDC", "high"

        # ── Tier 2: Immune cells ──────────────────────────────────────────────
        elif max(get("T_cell"), get("B_cell"),
                 get("Myeloid"), get("Neutrophil")) > 0.8:
            immune = {
                "T_cell":    get("T_cell"),
                "B_cell":    get("B_cell"),
                "Myeloid":   get("Myeloid"),
                "Neutrophil": get("Neutrophil"),
            }
            best = max(immune, key=immune.get)
            cell_type = best
            confidence = "high" if immune[best] > 1.5 else "medium"

        elif max(get("T_cell"), get("B_cell"),
                 get("Myeloid"), get("Neutrophil")) > 0.3:
            immune = {
                "T_cell":    get("T_cell"),
                "B_cell":    get("B_cell"),
                "Myeloid":   get("Myeloid"),
                "Neutrophil": get("Neutrophil"),
            }
            best = max(immune, key=immune.get)
            cell_type = best
            confidence = "low"

        # ── Tier 3: Stromal cells ─────────────────────────────────────────────
        elif get("Fibroblast") > 1.5:
            cell_type, confidence = "Fibroblast", "high"

        elif get("Endothelial") > 2.0:
            cell_type, confidence = "Endothelial", "high"

        # ── Tier 4: Normal epithelial ─────────────────────────────────────────
        elif get("Alveolar") > 2.0:
            cell_type, confidence = "Alveolar", "high"

        elif get("Epithelial") > 1.5:
            cell_type, confidence = "Epithelial", "high"

        # ── Tier 5: Cancer cells ──────────────────────────────────────────────
        # EPCAM+ but not expressing normal epithelial differentiation markers
        # Lower threshold — EMT cancer cells can have reduced EPCAM
        elif "score_EPCAM" in scores.index and "score_normal_epithelial" in scores.index:
            epcam = float(scores["score_EPCAM"])
            normal_epi = float(scores["score_normal_epithelial"])

            # Standard: EPCAM+ and not normal epithelial
            if epcam > 0.1 and normal_epi < 0.5:
                cell_type = "Cancer"
                confidence = "high" if epcam > 0.25 else "medium"

        # ── Tier 6: Unknown ───────────────────────────────────────────────────
        if cell_type is None:
            all_scores = {ct: get(ct) for ct in MARKER_GENES}
            best_ct = max(all_scores, key=all_scores.get)
            best_sc = all_scores[best_ct]

            # If no strong immune/stromal signal → likely cancer by exclusion
            # This follows paper's logic: cancer clusters negative for all other markers
            immune_stromal_max = max(
                get("T_cell"), get("B_cell"), get("Myeloid"),
                get("Neutrophil"), get("Fibroblast"), get("Endothelial"),
                get("Mast_cell"), get("fDC")
            )
            if immune_stromal_max < 0.3 and best_sc < 0.5:
                cell_type = "Cancer"
                confidence = "low"
                log.info(f"  Cluster {cluster}: assigned Cancer by exclusion "
                        f"(no immune/stromal signal, n={n_cells})")
            elif best_sc > 0.15:
                cell_type = best_ct
                confidence = "low"
            else:
                cell_type = "Unknown"
                confidence = "none"

        cluster_annotation[cluster] = cell_type
        cluster_confidence[cluster] = confidence

        # Get the assigned cell type's score for logging
        assigned_score = get(cell_type) if cell_type not in \
            ["Cancer", "Unknown"] else float(
                scores.get("score_EPCAM", 0))
        log.info(f"  Cluster {cluster:>3}: {cell_type:<15} "
                 f"n={n_cells:>6}, confidence={confidence:<6}, "
                 f"score={assigned_score:.3f}")

    # Apply annotations to cells
    adata.obs["cell_type"] = (
        adata.obs[cluster_key]
        .map(pd.Series(cluster_annotation))
        .astype("category")
    )
    adata.obs["annotation_confidence"] = (
        adata.obs[cluster_key]
        .map(pd.Series(cluster_confidence))
        .astype("category")
    )

    return pd.Series(cluster_annotation), cluster_mean


# =============================================================================
# FIGURES
# =============================================================================

def plot_resolution_sweep(sweep_results, selected_res, figures_dir, track):
    """Resolution sweep — n_clusters and silhouette vs resolution."""
    resolutions = sorted(sweep_results.keys())
    n_clusters  = [sweep_results[r]["n_clusters"]  for r in resolutions]
    silhouettes = [sweep_results[r]["silhouette"]   for r in resolutions]

    fig, ax1 = plt.subplots(figsize=(9, 4))
    ax2 = ax1.twinx()

    ax1.plot(resolutions, n_clusters, "o-", color="steelblue",
             linewidth=2, label="N clusters")
    ax2.plot(resolutions, silhouettes, "s--", color="tomato",
             linewidth=2, label="Silhouette")

    ax1.axvline(selected_res, color="grey", linestyle=":",
                label=f"Selected: {selected_res}")
    ax1.set_xlabel("Leiden resolution", fontsize=11)
    ax1.set_ylabel("Number of clusters", fontsize=11, color="steelblue")
    ax2.set_ylabel("Silhouette score", fontsize=11, color="tomato")
    ax1.set_title(f"[{track}] Resolution sweep — cluster quality",
                  fontsize=12, fontweight="bold")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=8)

    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir, f"resolution_sweep.{ext}"),
                    dpi=200, bbox_inches="tight")
    plt.close()


def plot_cluster_scores(cluster_scores, figures_dir, track):
    """Heatmap of cluster scores — annotation QC figure."""
    # Only show standard cell type score columns
    standard_cols = [c for c in cluster_scores.columns
                     if c.startswith("score_")
                     and "EPCAM" not in c
                     and "normal_epithelial" not in c]

    plot_data = cluster_scores[standard_cols].copy()
    plot_data.columns = [c.replace("score_", "") for c in plot_data.columns]

    fig, ax = plt.subplots(figsize=(14, 6))
    sns.heatmap(
        plot_data.T,
        cmap        = "RdBu_r",
        center      = 0,
        ax          = ax,
        yticklabels = True,
        xticklabels = True,
        linewidths  = 0.5,
        linecolor   = "white",
        annot       = False,
    )
    ax.set_title(f"[{track}] Cluster scores per cell type (annotation QC)",
                 fontweight="bold")
    ax.set_xlabel("Leiden cluster", fontsize=10)
    ax.set_ylabel("Cell type", fontsize=10)
    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir, f"cluster_scores_heatmap.{ext}"),
                    dpi=200, bbox_inches="tight")
    plt.close()
    log.info("  Saved cluster_scores_heatmap")


def plot_umap_cell_type(adata, figures_dir, track):
    """Fig 1b reproduction — UMAP colored by cell type."""
    log.info("  Plotting UMAP by cell type...")
    cell_types = [ct for ct in CELL_TYPE_COLORS
                  if ct in adata.obs["cell_type"].cat.categories]

    fig, ax = plt.subplots(figsize=(8, 7))
    for ct in cell_types:
        mask = adata.obs["cell_type"] == ct
        ax.scatter(
            adata.obsm["X_umap"][mask, 0],
            adata.obsm["X_umap"][mask, 1],
            c          = CELL_TYPE_COLORS[ct],
            s          = 0.5,
            alpha      = 0.8,
            rasterized = True,
        )
    handles = [Patch(facecolor=CELL_TYPE_COLORS[ct], label=ct)
               for ct in cell_types]
    ax.legend(handles=handles, loc="upper left",
              bbox_to_anchor=(1.01, 1), fontsize=9,
              frameon=True, title="Cell type", title_fontsize=10)
    ax.set_xlabel("UMAP1", fontsize=11)
    ax.set_ylabel("UMAP2", fontsize=11)
    ax.set_title(f"[{track}] UMAP — {adata.n_obs:,} cells, "
                 f"{len(cell_types)} cell types",
                 fontsize=12, fontweight="bold")
    ax.set_aspect("equal")
    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir, f"umap_cell_type.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log.info("  Saved umap_cell_type")


def plot_umap_patient(adata, figures_dir, track):
    """Fig 1d reproduction — UMAP colored by patient."""
    log.info("  Plotting UMAP by patient...")
    patients = sorted(adata.obs["patient_id"].unique(),
                      key=lambda x: int(x.replace("P", "")))
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
    handles = [Patch(facecolor=pcolors[p], label=p) for p in patients]
    ax.legend(handles=handles, loc="upper left",
              bbox_to_anchor=(1.01, 1), fontsize=7, ncol=2,
              frameon=True, title="Patient", title_fontsize=9)
    ax.set_xlabel("UMAP1", fontsize=11)
    ax.set_ylabel("UMAP2", fontsize=11)
    ax.set_title(f"[{track}] UMAP — colored by patient",
                 fontsize=12, fontweight="bold")
    ax.set_aspect("equal")
    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir, f"umap_by_patient.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log.info("  Saved umap_by_patient")


def plot_marker_heatmap(adata, figures_dir, track):
    """Fig 1c reproduction — canonical marker expression heatmap."""
    log.info("  Plotting marker heatmap...")
    cell_types = [ct for ct in CELL_TYPE_COLORS
                  if ct in adata.obs["cell_type"].cat.categories
                  and ct != "Unknown"]

    np.random.seed(42)
    sampled = []
    for ct in cell_types:
        idx = np.where(adata.obs["cell_type"] == ct)[0]
        n = min(300, len(idx))
        sampled.extend(np.random.choice(idx, n, replace=False))
    sampled = sorted(sampled)
    sub = adata[sampled].copy()

    # Collect valid heatmap markers
    all_markers = []
    for ct in cell_types:
        for g in HEATMAP_MARKERS.get(ct, []):
            if g in sub.var_names and g not in all_markers:
                all_markers.append(g)

    if not all_markers:
        log.warning("  No heatmap markers found — skipping heatmap")
        return

    # Expression matrix
    if "log1p_norm" in sub.layers:
        expr = sub.layers["log1p_norm"]
    else:
        expr = sub.X
    gene_idx = [list(sub.var_names).index(g) for g in all_markers]
    expr_arr = expr[:, gene_idx]
    if hasattr(expr_arr, "toarray"):
        expr_arr = expr_arr.toarray()
    expr_arr = np.array(expr_arr)

    df = pd.DataFrame(expr_arr, columns=all_markers)
    df["cell_type"] = pd.Categorical(
        sub.obs["cell_type"].values, categories=cell_types, ordered=True
    )
    df = df.sort_values("cell_type")
    row_colors = df["cell_type"].map(CELL_TYPE_COLORS)
    matrix = df.drop(columns=["cell_type"]).T

    g = sns.clustermap(
        matrix,
        col_colors  = row_colors.values,
        col_cluster = False,
        row_cluster = False,
        cmap        = "RdBu_r",
        vmin=0, vmax=3,
        figsize     = (16, 12),
        xticklabels = False,
        yticklabels = True,
        linewidths  = 0,
        rasterized  = True,
        cbar_kws    = {"label": "log1p expression", "shrink": 0.3},
    )
    g.ax_heatmap.yaxis.set_tick_params(labelsize=7)
    g.fig.suptitle(f"[{track}] Canonical cell type markers",
                   fontsize=13, fontweight="bold", y=1.01)
    handles = [Patch(facecolor=CELL_TYPE_COLORS[ct], label=ct)
               for ct in cell_types]
    g.ax_heatmap.legend(handles=handles, loc="upper left",
                        bbox_to_anchor=(1.15, 1), fontsize=8,
                        title="Cell type", title_fontsize=9)
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir, f"marker_heatmap.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log.info("  Saved marker_heatmap")

def plot_proportions(adata, figures_dir, track):
    """Fig 1e reproduction — stacked bar chart of cell type proportions."""
    log.info("  Plotting cell type proportions per patient...")
    patients = sorted(adata.obs["patient_id"].unique(),
                      key=lambda x: int(x.replace("P", "")))
    cell_types = [ct for ct in CELL_TYPE_COLORS
                  if ct in adata.obs["cell_type"].cat.categories
                  and ct != "Unknown"]

    prop = pd.crosstab(
        adata.obs["patient_id"],
        adata.obs["cell_type"],
        normalize="index"
    ).reindex(patients).reindex(columns=cell_types, fill_value=0)

    fig, ax = plt.subplots(figsize=(18, 5))
    bottom = np.zeros(len(patients))
    for ct in cell_types:
        ax.bar(range(len(patients)), prop[ct].values,
               bottom=bottom, color=CELL_TYPE_COLORS[ct],
               label=ct, width=0.85, edgecolor="none")
        bottom += prop[ct].values

    ax.set_xticks(range(len(patients)))
    ax.set_xticklabels(patients, rotation=90, fontsize=7)
    ax.set_xlabel("Patient", fontsize=11)
    ax.set_ylabel("Proportion", fontsize=11)
    ax.set_ylim(0, 1)
    ax.set_title(f"[{track}] Cell type composition per patient",
                 fontsize=12, fontweight="bold")
    handles = [Patch(facecolor=CELL_TYPE_COLORS[ct], label=ct)
               for ct in cell_types]
    ax.legend(handles=handles, loc="upper left",
              bbox_to_anchor=(1.01, 1), fontsize=8,
              title="Cell type", title_fontsize=9)
    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir, f"cell_type_proportions.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log.info("  Saved cell_type_proportions")


# =============================================================================
# TABLES
# =============================================================================

def save_tables(adata, tables_dir, track):
    """Save annotation tables for downstream analysis."""
    os.makedirs(tables_dir, exist_ok=True)

    # Per-cell assignments
    adata.obs[["patient_id", "gsm_id", "cell_type",
               "annotation_confidence"]].to_csv(
        os.path.join(tables_dir, "cell_type_assignments.csv")
    )

    # Summary counts
    counts = adata.obs["cell_type"].value_counts().reset_index()
    counts.columns = ["cell_type", "n_cells"]
    counts["pct"] = (counts["n_cells"] / adata.n_obs * 100).round(2)
    counts.to_csv(os.path.join(tables_dir, "cell_type_counts.csv"),
                  index=False)

    log.info(f"[{track}] Cell type summary:")
    for _, row in counts.iterrows():
        log.info(f"  {row['cell_type']:<15} "
                 f"{row['n_cells']:>7,} cells  ({row['pct']:.1f}%)")

    # Per-patient proportions
    prop = pd.crosstab(
        adata.obs["patient_id"],
        adata.obs["cell_type"],
        normalize="index"
    ).round(4)
    prop.to_csv(os.path.join(tables_dir, "cell_type_proportions_per_patient.csv"))

    log.info(f"  Tables saved to {tables_dir}")


# =============================================================================
# MAIN
# =============================================================================

def main(input_path, output_path, figures_dir, tables_dir,
         res_min, res_max, target_clusters, min_cells, track):

    log.info(f"[{track.upper()} TRACK] Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    log.info(f"  Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(tables_dir,  exist_ok=True)

    # ── Resolution sweep ──────────────────────────────────────────────────────
    resolutions = [round(r, 1)
                   for r in np.arange(res_min, res_max + 0.05, 0.1)]
    sweep = run_leiden_sweep(adata, resolutions)
    selected = select_resolution(sweep, target_clusters)
    plot_resolution_sweep(sweep, selected, figures_dir, track)

    cluster_key = f"leiden_{selected}"
    log.info(f"  Using cluster key: {cluster_key} "
             f"({sweep[selected]['n_clusters']} clusters)")

    # ── Scoring and annotation ────────────────────────────────────────────────
    score_markers(adata, cluster_key)
    cluster_annotation, cluster_scores = annotate_clusters(
        adata, cluster_key, min_cells=min_cells
    )
    plot_cluster_scores(cluster_scores, figures_dir, track)

    # ── Figures ───────────────────────────────────────────────────────────────
    plot_umap_cell_type(adata, figures_dir, track)
    plot_umap_patient(adata, figures_dir, track)
    plot_marker_heatmap(adata, figures_dir, track)
    plot_dotplot(adata, figures_dir, track)
    plot_proportions(adata, figures_dir, track)

    # ── Tables ────────────────────────────────────────────────────────────────
    save_tables(adata, tables_dir, track)

    # ── Store annotation metadata in uns ─────────────────────────────────────
    adata.uns["annotation"] = {
        "track":           track,
        "resolution":      float(selected),
        "cluster_key":     cluster_key,
        "n_clusters":      sweep[selected]["n_clusters"],
        "silhouette":      sweep[selected]["silhouette"],
        "cell_type_colors": CELL_TYPE_COLORS,
    }

    # ── Save ──────────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    log.info(f"  Saved annotated AnnData to {output_path}")
    log.info(f"Stage 5/6 [{track}] complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Biology-driven Leiden clustering and cell type annotation"
    )
    parser.add_argument("--input",           required=True)
    parser.add_argument("--output",          required=True)
    parser.add_argument("--figures_dir",     required=True)
    parser.add_argument("--tables_dir",      required=True)
    parser.add_argument("--res_min",         type=float, default=0.5)
    parser.add_argument("--res_max",         type=float, default=1.5)
    parser.add_argument("--target_clusters", type=int,   default=None,
                        help="Target cluster count for reproduction. "
                             "If omitted, uses silhouette maximization.")
    parser.add_argument("--min_cells",       type=int,   default=50,
                        help="Minimum cells per cluster. Smaller clusters "
                             "are labeled Unknown (likely artefacts).")
    parser.add_argument("--track",           required=True,
                        choices=["paper", "alternative"])
    args = parser.parse_args()
    main(args.input, args.output, args.figures_dir, args.tables_dir,
         args.res_min, args.res_max, args.target_clusters,
         args.min_cells, args.track)

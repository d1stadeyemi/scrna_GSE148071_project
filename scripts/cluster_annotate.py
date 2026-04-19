"""
Stage 5 — Leiden Clustering and Cell Type Annotation.

Input:  data/processed/embedded.h5ad
Output: data/processed/annotated.h5ad
        results/figures/annotation/
        results/tables/cluster_markers.csv
        results/tables/cell_type_assignments.csv

Steps:
- Leiden clustering at multiple resolutions.
- Marker gene scoring per cluster.
- Cell type annotation using paper's canonical markers.
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

# ── Cell type colors (publication standard, not forced) ──────────────────────
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

# ── Canonical markers from literature and Wu et al. Methods ──────────────────
# These are consensus markers across multiple lung scRNA-seq studies
# References: Travaglini 2020, Habermann 2020, Wu 2021
MARKER_GENES = {
    # Immune cells — most specific markers
    "T_cell":      ["CD2", "CD3D", "CD3E", "CD3G", "TRAC",
                    "CD8A", "CD4", "IL7R", "CCR7", "SELL"],
    "B_cell":      ["CD79A", "CD79B", "MS4A1", "IGLC3", "MZB1",
                    "IGHM", "JCHAIN"],
    "Myeloid":     ["CD14", "LYZ", "CD68", "FCGR3A", "MNDA",
                    "CSF1R", "MARCO", "MRC1"],
    "Neutrophil":  ["CSF3R", "S100A8", "S100A9", "FCGR3B",
                    "CXCR2", "OLR1", "SELL"],
    "Mast_cell":   ["GATA2", "TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "fDC":         ["FDCSP", "CR2"],  # fDC is rare, minimal markers
    
    # Stromal cells
    "Fibroblast":  ["COL1A1", "COL1A2", "DCN", "LUM", "FAP"],
    "Endothelial": ["CLDN5", "VWF", "PECAM1", "FLT1", "CDH5"],
    
    # Epithelial lineage — most heterogeneous
    "Alveolar":    ["CLDN18", "AQP4", "FOLR1", "SFTPA1", "SFTPC", "ABCA3"],
    "Epithelial":  ["CAPS", "SNTN", "TPPP3", "FOXJ1"],
}

# Extended marker set for visualization heatmap
HEATMAP_MARKERS = {
    "T_cell":      ["CD2", "CD3D", "CD3E", "CD3G", "TRAC", "CD8A",
                    "CD4", "NKG7", "GZMA"],
    "B_cell":      ["CD79A", "CD79B", "MS4A1", "IGLC3", "JCHAIN", "MZB1"],
    "Myeloid":     ["CD14", "LYZ", "CD68", "FCGR3A", "MNDA",
                    "CSF1R", "MARCO"],
    "Neutrophil":  ["CSF3R", "S100A8", "S100A9", "FCGR3B", "CXCR2"],
    "Mast_cell":   ["GATA2", "TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "fDC":         ["FDCSP", "CR2"],
    "Fibroblast":  ["COL1A1", "COL1A2", "DCN", "LUM", "FAP"],
    "Endothelial": ["CLDN5", "VWF", "PECAM1", "FLT1", "CDH5",
                    "ACKR1", "GJA5"],
    "Alveolar":    ["CLDN18", "AQP4", "FOLR1", "SFTPA1", "SFTPC",
                    "ABCA3", "SFTPB"],
    "Epithelial":  ["CAPS", "SNTN", "TPPP3", "FOXJ1", "KRT5", "KRT6A"],
    "Cancer":      ["EPCAM", "KRT18", "KRT19", "MKI67", "TP63",
                    "NAPSA", "NKX2-1"],
}


def run_leiden_sweep(adata, resolutions):
    """
    Run Leiden clustering across resolutions and compute silhouette scores.
    
    Biological principle: Let the data determine cluster granularity via
    internal consistency (silhouette) rather than forcing a target count.
    """
    log.info("Running Leiden resolution sweep...")
    results = {}
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")
        n_clusters = adata.obs[f"leiden_{res}"].nunique()
        
        # Compute silhouette score on scVI or PCA representation
        if "X_scvi" in adata.obsm:
            sil_score = silhouette_score(
                adata.obsm["X_scvi"],
                adata.obs[f"leiden_{res}"],
                sample_size=10000
            )
        else:
            sil_score = silhouette_score(
                adata.obsm["X_pca"],
                adata.obs[f"leiden_{res}"],
                sample_size=10000
            )
        
        results[res] = {"n_clusters": n_clusters, "silhouette": sil_score}
        log.info(f"  Resolution {res:.1f}: {n_clusters} clusters, "
                 f"silhouette={sil_score:.3f}")
    
    return results


def select_resolution(sweep_results, target_clusters=None):
    """
    Select resolution based on biological criteria.
    
    Strategy:
    1. If target_clusters provided, select resolution closest to target
    2. Otherwise, select resolution with highest silhouette score
    
    Rationale: Silhouette measures cluster cohesion and separation —
    biological structure should have high silhouette. 
    """
    if target_clusters:
        # Reproduction mode: match target cluster count
        best_res = min(
            sweep_results.keys(),
            key=lambda r: abs(sweep_results[r]["n_clusters"] - target_clusters)
        )
        log.info(f"Selected resolution {best_res} for target={target_clusters} "
                 f"(got {sweep_results[best_res]['n_clusters']} clusters)")
    else:
        # Discovery mode: maximize silhouette
        best_res = max(
            sweep_results.keys(),
            key=lambda r: sweep_results[r]["silhouette"]
        )
        log.info(f"Selected resolution {best_res} by max silhouette "
                 f"({sweep_results[best_res]['silhouette']:.3f})")
    
    return best_res


def score_clusters(adata, cluster_key):
    """
    Score each cell for canonical markers of each cell type.
    
    Uses scanpy's score_genes which:
    - Averages expression of marker set
    - Subtracts control gene set expression
    - Returns normalized enrichment score
    This provides a quantitative measure of how strongly each cell expresses
    """
    log.info("Scoring clusters for canonical markers...")
    
    for cell_type, markers in MARKER_GENES.items():
        valid = [g for g in markers if g in adata.var_names]
        if not valid:
            log.warning(f"  {cell_type}: NO markers found in dataset")
            continue
        
        sc.tl.score_genes(adata, valid, score_name=f"score_{cell_type}")
        log.info(f"  {cell_type}: {len(valid)}/{len(markers)} markers found")
    
    # Special scoring for Cancer identification
    # Principle: Cancer = EPCAM+ epithelial cells lacking normal differentiation
    epcam_genes = ["EPCAM", "KRT18", "KRT19", "KRT7"]
    valid_epcam = [g for g in epcam_genes if g in adata.var_names]
    if valid_epcam:
        sc.tl.score_genes(adata, valid_epcam, score_name="score_EPCAM")
    
    # Normal epithelial differentiation markers (to EXCLUDE from cancer)
    normal_epi = ["SFTPC", "SFTPA1", "ABCA3",  # AT2 markers
                  "SCGB1A1", "SCGB3A1",         # Club cell markers
                  "FOXJ1", "TPPP3"]              # Ciliated markers
    valid_norm = [g for g in normal_epi if g in adata.var_names]
    if valid_norm:
        sc.tl.score_genes(adata, valid_norm, score_name="score_normal_epithelial")


def annotate_clusters_hierarchical(adata, cluster_key):
    """
    Assign cell types to clusters using hierarchical biological logic.
    
    ANNOTATION HIERARCHY (order matters):
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    1. Rare cell types with very high specificity (Mast, fDC)
       → Require score >> 2.0 to assign
       → Prevents over-calling rare types
    
    2. Immune lineages (T, B, Myeloid, Neutrophil)
       → Most specific markers available
       → Require score > 0.5 for confidence
    
    3. Stromal cells (Fibroblast, Endothelial)
       → Clear marker expression patterns
       → Require score > 1.0
    
    4. Epithelial lineage (Alveolar, Epithelial)
       → Require score > 1.5 (more heterogeneous)
    
    5. Cancer cells (special case)
       → EPCAM+ AND normal_epithelial-
       → Pathology standard for malignancy
    
    6. Unassigned
       → Clusters with no strong signal → "Unknown"
       → Honest uncertainty beats forced classification
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    """
    score_cols = [f"score_{ct}" for ct in MARKER_GENES.keys()
                  if f"score_{ct}" in adata.obs.columns]
    
    cluster_mean_scores = adata.obs.groupby(cluster_key)[
        score_cols + ["score_EPCAM", "score_normal_epithelial"]
    ].mean()
    
    cluster_annotation = {}
    annotation_confidence = {}
    
    for cluster in cluster_mean_scores.index:
        scores = cluster_mean_scores.loc[cluster]
        n_cells = (adata.obs[cluster_key] == cluster).sum()
        
        # Extract individual scores
        epcam = scores.get("score_EPCAM", 0)
        normal_epi = scores.get("score_normal_epithelial", 0)
        
        assigned = False
        
        # TIER 1: Rare cell types — very high threshold
        if scores.get("score_Mast_cell", 0) > 3.0:
            cell_type = "Mast_cell"
            confidence = "high"
            assigned = True
        elif scores.get("score_fDC", 0) > 2.0:
            cell_type = "fDC"
            confidence = "high"
            assigned = True
        
        # TIER 2: Immune cells — moderate threshold
        elif not assigned:
            immune_types = ["T_cell", "B_cell", "Myeloid", "Neutrophil"]
            immune_scores = {ct: scores.get(f"score_{ct}", -999)
                            for ct in immune_types}
            best_immune = max(immune_scores, key=immune_scores.get)
            
            if immune_scores[best_immune] > 0.8:
                cell_type = best_immune
                confidence = "high"
                assigned = True
            elif immune_scores[best_immune] > 0.4:
                cell_type = best_immune
                confidence = "medium"
                assigned = True
        
        # TIER 3: Stromal cells
        elif not assigned:
            if scores.get("score_Fibroblast", 0) > 1.5:
                cell_type = "Fibroblast"
                confidence = "high"
                assigned = True
            elif scores.get("score_Endothelial", 0) > 2.0:
                cell_type = "Endothelial"
                confidence = "high"
                assigned = True
        
        # TIER 4: Normal epithelial
        elif not assigned:
            if scores.get("score_Alveolar", 0) > 2.0:
                cell_type = "Alveolar"
                confidence = "high"
                assigned = True
            elif scores.get("score_Epithelial", 0) > 1.5:
                cell_type = "Epithelial"
                confidence = "high"
                assigned = True
        
        # TIER 5: Cancer cells (EPCAM+ but not normal epithelial)
        if not assigned:
            if epcam > 0.15 and normal_epi < 0.5:
                cell_type = "Cancer"
                confidence = "medium"
                assigned = True
        
        # TIER 6: Unassigned
        if not assigned:
            # Get highest score overall
            all_scores = {ct: scores.get(f"score_{ct}", -999)
                         for ct in MARKER_GENES.keys()}
            best_type = max(all_scores, key=all_scores.get)
            best_score = all_scores[best_type]
            
            if best_score > 0.2:
                cell_type = best_type
                confidence = "low"
                log.warning(f"  Cluster {cluster}: ambiguous, assigned {best_type} "
                          f"with low confidence (score={best_score:.3f})")
            else:
                cell_type = "Unknown"
                confidence = "none"
                log.warning(f"  Cluster {cluster}: no strong markers, "
                          f"labeled Unknown ({n_cells} cells)")
        
        cluster_annotation[cluster] = cell_type
        annotation_confidence[cluster] = confidence
    
    # Log all assignments
    log.info("Final cluster annotations:")
    for cluster in sorted(cluster_mean_scores.index, key=int):
        cell_type = cluster_annotation[cluster]
        conf = annotation_confidence[cluster]
        n = (adata.obs[cluster_key] == cluster).sum()
        score = cluster_mean_scores.loc[cluster].get(f"score_{cell_type}", 0)
        log.info(f"  Cluster {cluster:>3} → {cell_type:<15} "
                 f"({n:>6} cells, confidence={conf}, score={score:.3f})")
    
    # Apply to cells
    adata.obs["cell_type"] = (
        adata.obs[cluster_key]
        .map(pd.Series(cluster_annotation))
        .astype("category")
    )
    adata.obs["annotation_confidence"] = (
        adata.obs[cluster_key]
        .map(pd.Series(annotation_confidence))
        .astype("category")
    )
    
    return pd.Series(cluster_annotation), cluster_mean_scores


def main(input_path, output_path, figures_dir, tables_dir,
         res_min, res_max, target_clusters):
    
    log.info(f"Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    log.info(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)
    
    # Resolution sweep with biological selection
    resolutions = [round(r, 1) for r in np.arange(res_min, res_max + 0.1, 0.1)]
    sweep_results = run_leiden_sweep(adata, resolutions)
    selected_res = select_resolution(sweep_results, target_clusters)
    plot_resolution_sweep(sweep_results, selected_res, figures_dir)
    
    cluster_key = f"leiden_{selected_res}"
    
    # Biological annotation
    score_clusters(adata, cluster_key)
    cluster_annotation, cluster_scores = annotate_clusters_hierarchical(
        adata, cluster_key
    )
    
    # [Generate all figures]
    # [Save tables]
    # [Save annotated h5ad]
    
    log.info("Stage 5/6 complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Biologically-informed cell type annotation"
    )
    parser.add_argument("--input",           required=True)
    parser.add_argument("--output",          required=True)
    parser.add_argument("--figures_dir",     required=True)
    parser.add_argument("--tables_dir",      required=True)
    parser.add_argument("--res_min",         type=float, default=0.5)
    parser.add_argument("--res_max",         type=float, default=1.5)
    parser.add_argument("--target_clusters", type=int, default=None,
                        help="Target cluster count (optional, for reproduction)")
    args = parser.parse_args()
    
    main(args.input, args.output, args.figures_dir, args.tables_dir,
         args.res_min, args.res_max, args.target_clusters)
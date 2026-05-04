"""
Stage 7 — CNV Inference and Malignant Cell Identification

Reproduces Wu et al. 2021 paper Fig. 2a (CNV heatmap) and Fig. 2b-2e
(intratumoral heterogeneity scores).

Biological rationale:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Cancer cells accumulate copy number variations (CNVs) — large-scale
chromosomal gains and losses driven by genomic instability. Normal cells
maintain diploid genomes. This difference can be detected from scRNA-seq
data because genes on gained chromosomal regions are over-expressed, and
genes on deleted regions are under-expressed, relative to normal cells.

infercnvpy implements the InferCNV algorithm:
1. Sort genes by genomic position
2. Smooth expression in sliding windows (101 genes default)
3. Subtract reference (immune/stromal) cell mean expression
4. Residuals = CNV signal: positive = gain, negative = loss

Reference cells (diploid ground truth):
  - T cells, B cells, Myeloid, Neutrophil, Fibroblast, Endothelial
  - These are not epithelial-derived and have normal diploid genomes
  - Mast cells and fDC excluded — too rare, unstable reference

Malignant classification:
  High CNV score → malignant (confirmed cancer)
  Low CNV score + EPCAM+ → normal epithelial contamination

Intratumoral heterogeneity (ITH):
  Paper defines two ITH metrics:
  - ITHCNA: IQR of pairwise Pearson correlation of CNV profiles
  - ITHGEX: IQR of pairwise Pearson correlation of gene expression profiles
  Higher IQR = more heterogeneous tumor

Reference:
  Wu et al. 2021, Nature Communications
  infercnvpy: Sturm et al., Bioinformatics 2023
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
import seaborn as sns

try:
    import infercnvpy as cnv
    HAS_INFERCNVPY = True
except ImportError:
    HAS_INFERCNVPY = False

logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s [%(levelname)s] %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger(__name__)


# ── Reference cell types (known diploid) ─────────────────────────────────────
# These will be used to establish the baseline expression profile.
# ALL must be non-epithelial (immune and stromal only).
REFERENCE_CELL_TYPES = [
    "T_cell", "B_cell", "Myeloid", "Neutrophil",
    "Fibroblast", "Endothelial"
]

# Cell types to analyze for CNV (epithelial origin)
EPITHELIAL_CELL_TYPES = ["Cancer", "Alveolar", "Epithelial"]

# CNV score threshold for malignancy classification
# Cells above this threshold are classified as confirmed malignant
# Based on InferCNV standard practice and paper's approach
CNV_SCORE_THRESHOLD = 0.03  # tunable parameter


def load_gene_positions(gtf_path=None):
    """
    Load gene chromosomal positions for CNV inference.

    infercnvpy requires gene positions to sort by genomic location.
    We use the GENCODE/Ensembl annotation that matches the paper's GRCh38 v92.

    If GTF not available, infercnvpy can use its built-in gene position database.
    """
    if gtf_path and os.path.exists(gtf_path):
        log.info(f"  Loading gene positions from {gtf_path}")
        # infercnvpy reads GTF directly
        return gtf_path
    else:
        log.info("  Using infercnvpy built-in gene positions (GRCh38)")
        return None


def prepare_adata_for_cnv(adata):
    """
    Prepare AnnData for infercnvpy:
    1. Subset to epithelial + reference cells
    2. Use raw counts (infercnvpy normalizes internally)
    3. Add cnv_label obs column
    """
    log.info("Preparing AnnData for CNV inference...")

    # Check which cell types are present
    present = adata.obs["cell_type"].unique()
    ref_present = [ct for ct in REFERENCE_CELL_TYPES if ct in present]
    epi_present = [ct for ct in EPITHELIAL_CELL_TYPES if ct in present]

    log.info(f"  Reference cell types: {ref_present}")
    log.info(f"  Epithelial cell types: {epi_present}")

    # Subset to relevant cells
    mask = adata.obs["cell_type"].isin(ref_present + epi_present)
    adata_sub = adata[mask].copy()
    log.info(f"  Subset: {adata_sub.n_obs:,} cells x {adata_sub.n_vars:,} genes")

    # Add cnv_label: reference vs epithelial/cancer
    adata_sub.obs["cnv_label"] = adata_sub.obs["cell_type"].apply(
        lambda ct: "reference" if ct in REFERENCE_CELL_TYPES else ct
    )

    # Use raw counts for CNV inference
    if "counts" in adata_sub.layers:
        adata_sub.X = adata_sub.layers["counts"].copy()
        log.info("  Using raw counts from layers['counts']")
    else:
        log.warning("  No raw counts layer found — using adata.X")

    n_ref = (adata_sub.obs["cnv_label"] == "reference").sum()
    n_epi = (adata_sub.obs["cnv_label"] != "reference").sum()
    log.info(f"  Reference cells: {n_ref:,}")
    log.info(f"  Epithelial/Cancer cells: {n_epi:,}")

    return adata_sub


def run_infercnv(adata_sub, gtf_path=None):
    """
    Run infercnvpy CNV inference.

    Parameters matching paper's InferCNV settings (Methods):
    - window_size=101 genes (paper: "101 genes as slide window")
    - step=1
    - reference: immune + stromal cells
    - Outputs: adata.obsm['X_cnv']
    """
    log.info("Running infercnvpy CNV inference...")

    if not HAS_INFERCNVPY:
        log.error("infercnvpy not installed — cannot run CNV inference")
        log.error("Add 'infercnvpy' to workflow/envs/scvi.yaml and rebuild env")
        raise ImportError("infercnvpy required for Stage 7")

    # Annotate gene positions
    # infercnvpy uses gtf or its built-in database
    if gtf_path:
        cnv.io.genomic_position_from_gtf(gtf_path, adata_sub, gtf_gene_id="gene_name")
    else:
        # Use built-in — matches GRCh38
        try:
            cnv.io.genomic_position_from_gtf(
                "grch38", adata_sub, gtf_gene_id="gene_name"
            )
        except Exception:
            log.warning("  Built-in GTF failed — trying gene symbol lookup")
            cnv.io.genomic_position_from_gtf(
                "grch38_110", adata_sub, gtf_gene_id="gene_name"
            )

    # Filter to genes with position info
    adata_sub = adata_sub[:, ~adata_sub.var["chromosome"].isnull()].copy()
    log.info(f"  Genes with position info: {adata_sub.n_vars:,}")

    # Run CNV inference
    # Paper Methods: "101 genes as a slide window to smooth relative expression"
    log.info("  Running infercnvpy (window_size=101, step=1)...")
    cnv.tl.infercnv(
        adata_sub,
        reference_key    = "cnv_label",
        reference_cat    = "reference",
        window_size      = 101,
        step             = 1,
        n_jobs           = 8,
        inplace          = True,
    )

    log.info(f"  CNV matrix shape: {adata_sub.obsm['X_cnv'].shape}")
    return adata_sub


def compute_cnv_scores(adata_cnv):
    """
    Compute per-cell CNV scores for malignancy classification.

    CNV score = mean of squared CNV values across all genomic positions.
    High score → large deviations from diploid → malignant cell.

    Paper: "1.5 standard deviation of the residual normalized expression values
    as the ceiling and floor for visualization"
    """
    log.info("Computing CNV scores...")

    X_cnv = adata_cnv.obsm["X_cnv"]
    if hasattr(X_cnv, "toarray"):
        X_cnv = X_cnv.toarray()

    # CNV score = mean squared deviation from 0 (diploid)
    cnv_score = np.mean(X_cnv ** 2, axis=1)
    adata_cnv.obs["cnv_score"] = cnv_score

    log.info(f"  CNV score range: {cnv_score.min():.4f} — {cnv_score.max():.4f}")
    log.info(f"  CNV score mean: {cnv_score.mean():.4f}")
    log.info(f"  CNV score median: {np.median(cnv_score):.4f}")

    # Reference cells should have low CNV scores — validate
    ref_scores = cnv_score[adata_cnv.obs["cnv_label"] == "reference"]
    epi_scores = cnv_score[adata_cnv.obs["cnv_label"] != "reference"]
    log.info(f"  Reference CNV mean: {ref_scores.mean():.4f} (should be low)")
    log.info(f"  Epithelial/Cancer CNV mean: {epi_scores.mean():.4f} (should be high)")

    return adata_cnv


def classify_malignant(adata_cnv, threshold=CNV_SCORE_THRESHOLD):
    """
    Classify epithelial cells as malignant vs non-malignant based on CNV score.

    This implements the paper's cancer cell validation:
    "We used non-malignant cells including immune cells and stromal cells
    as baselines to estimate the CNA of malignant cells"

    Cells originally labeled Cancer with LOW CNV → likely normal epithelium
    Cells originally labeled Cancer with HIGH CNV → confirmed malignant
    """
    log.info(f"Classifying malignant cells (threshold={threshold})...")

    adata_cnv.obs["malignant"] = "reference"

    epi_mask = adata_cnv.obs["cnv_label"] != "reference"
    high_cnv = adata_cnv.obs["cnv_score"] > threshold

    adata_cnv.obs.loc[epi_mask & high_cnv, "malignant"] = "malignant"
    adata_cnv.obs.loc[epi_mask & ~high_cnv, "malignant"] = "non_malignant_epithelial"

    # Summary
    n_malignant = (adata_cnv.obs["malignant"] == "malignant").sum()
    n_normal    = (adata_cnv.obs["malignant"] == "non_malignant_epithelial").sum()
    n_ref       = (adata_cnv.obs["malignant"] == "reference").sum()

    log.info(f"  Malignant: {n_malignant:,} cells")
    log.info(f"  Non-malignant epithelial: {n_normal:,} cells")
    log.info(f"  Reference (immune/stromal): {n_ref:,} cells")

    return adata_cnv


def compute_ith_scores(adata_cnv):
    """
    Compute intratumoral heterogeneity (ITH) scores per patient.

    Paper Methods:
    "ITHCNA: IQR of the distribution for all malignant cell pairs'
     Pearson correlation coefficients"
    "ITHGEX: IQR of the distribution using gene expression profiles"

    Higher IQR = more heterogeneous tumor (cells are more different from each other).
    """
    log.info("Computing intratumoral heterogeneity (ITH) scores...")

    malignant_mask = adata_cnv.obs["malignant"] == "malignant"
    adata_mal = adata_cnv[malignant_mask]

    X_cnv = adata_mal.obsm["X_cnv"]
    if hasattr(X_cnv, "toarray"):
        X_cnv = X_cnv.toarray()

    patients = sorted(
        adata_mal.obs["patient_id"].unique(),
        key=lambda x: int(x.replace("P", ""))
    )

    ith_results = []

    for patient in patients:
        pmask = adata_mal.obs["patient_id"] == patient
        n_cells = pmask.sum()

        if n_cells < 10:
            log.warning(f"  {patient}: only {n_cells} malignant cells — skipping ITH")
            continue

        # Subsample to max 100 cells per patient (paper: "100 malignant cells")
        idx = np.where(pmask)[0]
        n_sample = min(100, len(idx))
        np.random.seed(42)
        sampled = np.random.choice(idx, n_sample, replace=False)

        # ITHCNA: IQR of pairwise Pearson correlation on CNV profiles
        cnv_mat = X_cnv[sampled]
        corr_cnv = np.corrcoef(cnv_mat)
        upper_tri = corr_cnv[np.triu_indices_from(corr_cnv, k=1)]
        ith_cna = np.percentile(upper_tri, 75) - np.percentile(upper_tri, 25)

        # ITHGEX: IQR on log-normalized expression profiles
        if "log1p_norm" in adata_mal.layers:
            gex_mat = adata_mal.layers["log1p_norm"][sampled]
        else:
            gex_mat = adata_mal.X[sampled]
        if hasattr(gex_mat, "toarray"):
            gex_mat = gex_mat.toarray()
        corr_gex = np.corrcoef(gex_mat)
        upper_gex = corr_gex[np.triu_indices_from(corr_gex, k=1)]
        ith_gex = np.percentile(upper_gex, 75) - np.percentile(upper_gex, 25)

        ith_results.append({
            "patient_id": patient,
            "n_malignant": n_cells,
            "n_sampled":   n_sample,
            "ITHCNA":      round(ith_cna, 4),
            "ITHGEX":      round(ith_gex, 4),
        })
        log.info(f"  {patient}: n={n_cells:>4}, "
                 f"ITHCNA={ith_cna:.4f}, ITHGEX={ith_gex:.4f}")

    ith_df = pd.DataFrame(ith_results)
    log.info(f"\n  ITH summary:")
    log.info(f"  Mean ITHCNA: {ith_df['ITHCNA'].mean():.4f}")
    log.info(f"  Mean ITHGEX: {ith_df['ITHGEX'].mean():.4f}")

    return ith_df


# =============================================================================
# FIGURES
# =============================================================================

def plot_cnv_heatmap(adata_cnv, figures_dir, n_cells_per_patient=100):
    """
    Reproduce paper Fig. 2a — CNV heatmap.

    Paper: "randomly sampled 100 malignant cells of each patient
    were shown as their representative CNA profiles"
    """
    log.info("Plotting CNV heatmap (Fig. 2a)...")

    malignant_mask = adata_cnv.obs["malignant"] == "malignant"
    adata_mal = adata_cnv[malignant_mask]

    X_cnv = adata_mal.obsm["X_cnv"]
    if hasattr(X_cnv, "toarray"):
        X_cnv = X_cnv.toarray()

    patients = sorted(
        adata_mal.obs["patient_id"].unique(),
        key=lambda x: int(x.replace("P", ""))
    )

    # Sample cells per patient
    np.random.seed(42)
    sampled_idx = []
    sampled_patients = []

    for patient in patients:
        pmask = np.where(adata_mal.obs["patient_id"] == patient)[0]
        n = min(n_cells_per_patient, len(pmask))
        if n < 5:
            continue
        chosen = np.random.choice(pmask, n, replace=False)
        sampled_idx.extend(chosen)
        sampled_patients.extend([patient] * n)

    cnv_matrix = X_cnv[sampled_idx]

    # Clip at 1.5 SD as in paper
    # "1.5 standard deviation of the residual normalized expression values
    # as the ceiling and floor for visualization"
    sd = np.std(cnv_matrix)
    vmax = 1.5 * sd
    cnv_clipped = np.clip(cnv_matrix, -vmax, vmax)

    # Patient color bar
    unique_patients = sorted(set(sampled_patients),
                             key=lambda x: int(x.replace("P", "")))
    cmap_p = plt.cm.get_cmap("tab20", 20)
    cmap_p2 = plt.cm.get_cmap("tab20b", 20)
    pool = [cmap_p(i) for i in range(20)] + [cmap_p2(i) for i in range(20)]
    patient_colors = {p: pool[i] for i, p in enumerate(unique_patients)}
    row_colors = [patient_colors[p] for p in sampled_patients]

    fig, axes = plt.subplots(
        2, 1,
        figsize=(16, 10),
        gridspec_kw={"height_ratios": [0.5, 10]},
    )

    # Patient color bar
    patient_bar = np.array([
        list(patient_colors[p]) for p in sampled_patients
    ])
    axes[0].imshow(
        patient_bar.reshape(1, len(sampled_patients), 4),
        aspect="auto", interpolation="nearest"
    )
    axes[0].set_yticks([0])
    axes[0].set_yticklabels(["Patient"], fontsize=8)
    axes[0].set_xticks([])

    # CNV heatmap
    im = axes[1].imshow(
        cnv_clipped.T,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vmax, vmax=vmax,
        interpolation="nearest",
        rasterized=True,
    )
    axes[1].set_xlabel("Cells (grouped by patient)", fontsize=10)
    axes[1].set_ylabel("Chromosomal position →", fontsize=10)
    axes[1].set_yticks([])

    plt.colorbar(im, ax=axes[1], label="CNV signal", shrink=0.3,
                 orientation="vertical")

    # Patient legend
    handles = [
        Patch(facecolor=patient_colors[p], label=p)
        for p in unique_patients[:21]  # max 21 in legend
    ]
    axes[0].legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        fontsize=6,
        ncol=2,
        title="Patient",
        title_fontsize=7,
    )

    plt.suptitle("CNV profiles of malignant cells per patient (Fig. 2a)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()

    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir, f"cnv_heatmap.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log.info("  Saved cnv_heatmap")


def plot_cnv_scores_distribution(adata_cnv, figures_dir):
    """
    CNV score distribution — shows separation between malignant and normal cells.
    """
    log.info("Plotting CNV score distribution...")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Distribution by malignancy status
    colors = {
        "malignant":               "#E6851E",
        "non_malignant_epithelial": "#9B59B6",
        "reference":               "#4A90C4",
    }
    labels = {
        "malignant":               "Malignant",
        "non_malignant_epithelial": "Non-malignant epithelial",
        "reference":               "Reference (immune/stromal)",
    }

    ax = axes[0]
    for status, color in colors.items():
        mask = adata_cnv.obs["malignant"] == status
        if mask.sum() == 0:
            continue
        scores = adata_cnv.obs.loc[mask, "cnv_score"]
        ax.hist(scores, bins=50, alpha=0.6, color=color,
                label=f"{labels[status]} (n={mask.sum():,})",
                density=True)
    ax.axvline(CNV_SCORE_THRESHOLD, color="red", linestyle="--",
               label=f"Threshold={CNV_SCORE_THRESHOLD}")
    ax.set_xlabel("CNV score (mean squared deviation)", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_title("CNV score distribution by cell status", fontsize=11,
                 fontweight="bold")
    ax.legend(fontsize=8)

    # CNV score on UMAP (if UMAP is available)
    ax2 = axes[1]
    if "X_umap" in adata_cnv.obsm:
        sc_obj = ax2.scatter(
            adata_cnv.obsm["X_umap"][:, 0],
            adata_cnv.obsm["X_umap"][:, 1],
            c=adata_cnv.obs["cnv_score"],
            s=0.5, alpha=0.6,
            cmap="Reds",
            rasterized=True,
        )
        plt.colorbar(sc_obj, ax=ax2, label="CNV score", shrink=0.7)
        ax2.set_xlabel("UMAP1", fontsize=10)
        ax2.set_ylabel("UMAP2", fontsize=10)
        ax2.set_title("CNV score on UMAP", fontsize=11, fontweight="bold")
        ax2.set_aspect("equal")
    else:
        ax2.text(0.5, 0.5, "UMAP not available\nin CNV subset",
                 ha="center", va="center", transform=ax2.transAxes)

    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir, f"cnv_score_distribution.{ext}"),
                    dpi=200, bbox_inches="tight")
    plt.close()
    log.info("  Saved cnv_score_distribution")


def plot_ith_scores(ith_df, figures_dir):
    """
    Plot ITHCNA and ITHGEX per patient — reproduces paper Fig. 2b/2c.
    """
    log.info("Plotting ITH scores per patient...")

    fig, axes = plt.subplots(1, 2, figsize=(16, 4))

    for ax, metric, color in zip(
        axes,
        ["ITHCNA", "ITHGEX"],
        ["#E6851E", "#4A90C4"]
    ):
        patients = ith_df["patient_id"].tolist()
        values   = ith_df[metric].tolist()

        ax.bar(range(len(patients)), values, color=color, alpha=0.8,
               edgecolor="none")
        ax.set_xticks(range(len(patients)))
        ax.set_xticklabels(patients, rotation=90, fontsize=7)
        ax.set_xlabel("Patient", fontsize=10)
        ax.set_ylabel(metric, fontsize=10)
        ax.set_title(f"{metric} per patient\n"
                     f"(IQR of pairwise Pearson correlation)",
                     fontsize=11, fontweight="bold")
        ax.axhline(ith_df[metric].mean(), color="red", linestyle="--",
                   alpha=0.7, label=f"Mean={ith_df[metric].mean():.3f}")
        ax.legend(fontsize=8)

    plt.suptitle("Intratumoral heterogeneity scores (ITH) per patient",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir, f"ith_scores.{ext}"),
                    dpi=200, bbox_inches="tight")
    plt.close()
    log.info("  Saved ith_scores")


def plot_malignant_umap(adata_cnv, adata_full, figures_dir):
    """
    UMAP showing malignancy status — confirms CNV-based classification.
    """
    log.info("Plotting malignancy UMAP...")

    if "X_umap" not in adata_cnv.obsm:
        log.warning("  No UMAP in CNV subset — skipping malignancy UMAP")
        return

    colors = {
        "malignant":               "#E6851E",
        "non_malignant_epithelial": "#9B59B6",
        "reference":               "#AAAAAA",
    }
    labels = {
        "malignant":               "Malignant cancer",
        "non_malignant_epithelial": "Non-malignant epithelial",
        "reference":               "Reference (immune/stromal)",
    }

    fig, ax = plt.subplots(figsize=(8, 7))

    for status in ["reference", "non_malignant_epithelial", "malignant"]:
        mask = adata_cnv.obs["malignant"] == status
        if mask.sum() == 0:
            continue
        ax.scatter(
            adata_cnv.obsm["X_umap"][mask, 0],
            adata_cnv.obsm["X_umap"][mask, 1],
            c=colors[status], s=0.5, alpha=0.7,
            rasterized=True, label=labels[status],
        )

    handles = [Patch(facecolor=colors[s], label=labels[s])
               for s in colors]
    ax.legend(handles=handles, loc="upper left",
              bbox_to_anchor=(1.01, 1), fontsize=9)
    ax.set_xlabel("UMAP1", fontsize=11)
    ax.set_ylabel("UMAP2", fontsize=11)
    ax.set_title("Malignancy classification by CNV score",
                 fontsize=12, fontweight="bold")
    ax.set_aspect("equal")

    plt.tight_layout()
    for ext in ["pdf", "png"]:
        plt.savefig(os.path.join(figures_dir, f"malignancy_umap.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()
    log.info("  Saved malignancy_umap")


# =============================================================================
# MAIN
# =============================================================================

def main(input_path, output_path, figures_dir, tables_dir,
         gtf_path, cnv_threshold, n_cells_heatmap):

    log.info(f"Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    log.info(f"  Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    log.info(f"  Cell types: {dict(adata.obs['cell_type'].value_counts())}")

    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(tables_dir,  exist_ok=True)

    # ── Prepare CNV subset ────────────────────────────────────────────────────
    adata_cnv = prepare_adata_for_cnv(adata)

    # ── Run infercnvpy ────────────────────────────────────────────────────────
    adata_cnv = run_infercnv(adata_cnv, gtf_path)

    # ── CNV scores and classification ─────────────────────────────────────────
    adata_cnv = compute_cnv_scores(adata_cnv)
    adata_cnv = classify_malignant(adata_cnv, threshold=cnv_threshold)

    # ── ITH scores ────────────────────────────────────────────────────────────
    ith_df = compute_ith_scores(adata_cnv)

    # ── Figures ───────────────────────────────────────────────────────────────
    plot_cnv_heatmap(adata_cnv, figures_dir, n_cells_heatmap)
    plot_cnv_scores_distribution(adata_cnv, figures_dir)
    plot_ith_scores(ith_df, figures_dir)
    plot_malignant_umap(adata_cnv, adata, figures_dir)

    # ── Save tables ───────────────────────────────────────────────────────────
    # Per-cell malignancy status
    malignancy_table = adata_cnv.obs[[
        "patient_id", "cell_type", "cnv_label",
        "cnv_score", "malignant"
    ]].copy()
    malignancy_table.to_csv(
        os.path.join(tables_dir, "malignancy_classification.csv")
    )
    log.info("  Saved malignancy_classification.csv")

    # ITH scores
    ith_df.to_csv(
        os.path.join(tables_dir, "ith_scores.csv"), index=False
    )
    log.info("  Saved ith_scores.csv")

    # ── Propagate malignancy back to full AnnData ─────────────────────────────
    # Add CNV score and malignancy status to the full annotated object
    adata.obs["cnv_score"]  = np.nan
    adata.obs["malignant"]  = "not_analyzed"

    shared_idx = adata_cnv.obs_names.intersection(adata.obs_names)
    adata.obs.loc[shared_idx, "cnv_score"] = \
        adata_cnv.obs.loc[shared_idx, "cnv_score"]
    adata.obs.loc[shared_idx, "malignant"] = \
        adata_cnv.obs.loc[shared_idx, "malignant"]

    # Refine cell_type: confirmed_malignant replaces Cancer where CNV confirms
    refine_mask = (
        (adata.obs["cell_type"] == "Cancer") &
        (adata.obs["malignant"] == "malignant")
    )
    log.info(f"  Confirmed malignant cells: {refine_mask.sum():,}")

    # ── Save updated AnnData ──────────────────────────────────────────────────
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    log.info(f"  Saved updated AnnData to {output_path}")
    log.info("Stage 7 complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="CNV inference and malignant cell identification"
    )
    parser.add_argument("--input",           required=True,
                        help="Annotated AnnData from Stage 5/6")
    parser.add_argument("--output",          required=True,
                        help="Updated AnnData with CNV annotations")
    parser.add_argument("--figures_dir",     required=True)
    parser.add_argument("--tables_dir",      required=True)
    parser.add_argument("--gtf_path",        default=None,
                        help="Path to GTF file for gene positions. "
                             "If not provided, uses infercnvpy built-in GRCh38")
    parser.add_argument("--cnv_threshold",   type=float, default=0.03,
                        help="CNV score threshold for malignancy. "
                             "Cells above threshold = malignant")
    parser.add_argument("--n_cells_heatmap", type=int,   default=100,
                        help="Cells per patient in CNV heatmap (paper: 100)")
    args = parser.parse_args()
    main(args.input, args.output, args.figures_dir, args.tables_dir,
         args.gtf_path, args.cnv_threshold, args.n_cells_heatmap)
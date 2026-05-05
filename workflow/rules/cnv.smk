# =============================================================================
# Stage 7 — CNV Inference and Malignant Cell Identification
# =============================================================================
# Reproduces Wu et al. 2021 Fig. 2a (CNV heatmap) and Fig. 2b-2e (ITH scores)
#
# Biological rationale:
#   Cancer cells accumulate chromosomal copy number alterations (CNAs).
#   Normal (immune, stromal) cells are diploid — used as reference baseline.
#   CNV inference detects gene expression deviations consistent with
#   chromosomal gains/losses, distinguishing truly malignant cells from
#   normal epithelial contamination in our Cancer annotation.
#
# Method: infercnvpy (Python reimplementation of InferCNV)
#   - Window: 101 genes (paper exact match)
#   - Reference: T_cell, B_cell, Myeloid, Neutrophil, Fibroblast, Endothelial
#   - Gene positions: GRCh38 (matches paper's ensemble v92)
#
# Input:  paper_track/annotated.h5ad
# Output: paper_track/cnv_annotated.h5ad
#         results/figures/paper_track/cnv/
#         results/tables/paper_track/cnv/
# =============================================================================

rule infer_cnv:
    """
    CNV inference using infercnvpy.

    Runtime estimate: 60-120 min for 89,819 cells on campus-bigmem.
    Memory: ~64GB sufficient for this dataset size.

    Note: Uses paper track only — batch correction would confound CNV inference
    by harmonizing the very expression differences that indicate CNVs.
    """
    input:
        h5ad = "data/processed/paper_track/annotated.h5ad"
    output:
        h5ad    = "data/processed/paper_track/cnv_annotated.h5ad",
        figures = directory("results/figures/paper_track/cnv"),
        tables  = directory("results/tables/paper_track/cnv")
    params:
        gtf_path        = config["cnv"].get("gtf_path", ""),
        cnv_threshold   = config["cnv"]["cnv_threshold"],
        n_cells_heatmap = config["cnv"]["n_cells_heatmap"]
    threads: 8
    resources:
        mem_mb          = 128000,        # 128GB — CNV matrix is large
        runtime         = 300,           # 6 hours to be safe
        slurm_partition = "campus"
    conda:
        "../envs/cnv.yaml"               # separate env with infercnvpy
    log:
        "logs/paper_track/infer_cnv.log"
    shell:
        """
        python scripts/infer_cnv.py \
            --input           {input.h5ad} \
            --output          {output.h5ad} \
            --figures_dir     {output.figures} \
            --tables_dir      {output.tables} \
            --cnv_threshold   {params.cnv_threshold} \
            --n_cells_heatmap {params.n_cells_heatmap} \
            2>&1 | tee {log}
        """

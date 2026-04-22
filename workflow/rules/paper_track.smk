# =============================================================================
# Paper Track — Faithful reproduction of Wu et al. 2021
# =============================================================================
# No batch correction, 600 HVGs, 20 PCs, PCA-based UMAP
#
# Biological rationale:
#   The paper explicitly chose no batch correction to preserve cancer cell
#   heterogeneity. Advanced NSCLC tumors evolve independently in each patient.
#   Somatic mutations and clonal evolution create real transcriptional
#   differences that batch correction would erroneously remove.
# =============================================================================


rule paper_normalize_hvg:
    """
    Normalization and HVG selection for paper-faithful track.
    
    Follows Wu et al. Methods:
    - NormalizeData (target sum = 10,000)
    - Log1p transformation
    - 600 HVGs using Seurat 2.3 FindVariable equivalent (seurat flavor)
    - 20 PCA components (paper uses first 20 PCs)
    
    Note: seurat flavor (not seurat_v3) does NOT require raw counts layer.
    It uses the Fano factor-based dispersion method of Seurat 2.3.
    """
    input:
        h5ad = "data/processed/qc_filtered.h5ad"
    output:
        h5ad    = "data/processed/paper_track/normalized.h5ad",
        figures = directory("results/figures/paper_track/hvg_pca")
    params:
        n_hvgs      = config["paper_track"]["n_hvgs"],
        n_pcs       = config["paper_track"]["n_pcs"],
        hvg_flavor  = config["paper_track"]["hvg_flavor"]
    threads: 8
    resources:
        mem_mb          = 64000,
        runtime         = 30,
        slurm_partition = "campus-bigmem"
    conda:
        "../envs/scvi.yaml"
    log:
        "logs/paper_track/normalize_hvg.log"
    shell:
        """
        mkdir -p logs/paper_track
        python scripts/normalize_hvg.py \
            --input       {input.h5ad} \
            --output      {output.h5ad} \
            --figures_dir {output.figures} \
            --n_hvgs      {params.n_hvgs} \
            --n_pcs       {params.n_pcs} \
            --hvg_flavor  {params.hvg_flavor} \
            --track       paper \
            2>&1 | tee {log}
        """


rule paper_umap_embed:
    """
    UMAP embedding for paper track — uses PCA directly, no batch correction.
    
    Biological rationale:
    Using X_pca (no batch correction) means:
    - Cancer cells will cluster by patient (patient-specific biology preserved)
    - Immune/stromal cells will mix across patients (shared biology)
    This matches Fig 1d of Wu et al. and validates their batch correction decision.
    """
    input:
        h5ad = "data/processed/paper_track/normalized.h5ad"
    output:
        h5ad    = "data/processed/paper_track/embedded.h5ad",
        figures = directory("results/figures/paper_track/umap")
    params:
        n_neighbors = config["paper_track"]["n_neighbors"],
        n_pcs       = config["paper_track"]["n_pcs"],
        use_rep     = config["paper_track"]["use_rep"]
    threads: 8
    resources:
        mem_mb          = 64000,
        runtime         = 30,
        slurm_partition = "campus-bigmem"
    conda:
        "../envs/scvi.yaml"
    log:
        "logs/paper_track/umap_embed.log"
    shell:
        """
        python scripts/umap_embed.py \
            --input       {input.h5ad} \
            --output      {output.h5ad} \
            --figures_dir {output.figures} \
            --n_neighbors {params.n_neighbors} \
            --n_pcs       {params.n_pcs} \
            --use_rep     {params.use_rep} \
            --track       paper \
            2>&1 | tee {log}
        """


rule paper_cluster_annotate:
    """
    Leiden clustering and biologically-informed cell type annotation
    for the paper track.
    
    Annotation follows biological hierarchy — see scripts/cluster_annotate.py
    for full documentation of the annotation strategy.
    """
    input:
        h5ad = "data/processed/paper_track/embedded.h5ad"
    output:
        h5ad    = "data/processed/paper_track/annotated.h5ad",
        figures = directory("results/figures/paper_track/annotation")
    params:
        tables_dir      = "results/tables/paper_track",
        res_min         = config["annotation"]["res_min"],
        res_max         = config["annotation"]["res_max"],
        target_clusters = config["annotation"]["target_clusters"],
        min_cells       = config["annotation"]["min_cells_per_cluster"]
    threads: 8
    resources:
        mem_mb          = 64000,
        runtime         = 90,
        slurm_partition = "campus-bigmem"
    conda:
        "../envs/scvi.yaml"
    log:
        "logs/paper_track/cluster_annotate.log"
    shell:
        """
        python scripts/cluster_annotate.py \
            --input           {input.h5ad} \
            --output          {output.h5ad} \
            --figures_dir     {output.figures} \
            --tables_dir      {params.tables_dir} \
            --res_min         {params.res_min} \
            --res_max         {params.res_max} \
            --target_clusters {params.target_clusters} \
            --min_cells       {params.min_cells} \
            --track           paper \
            2>&1 | tee {log}
        """

# ── Track A: According to the paper (no batch correction) ───────────────────
# Follows Wu et al. 2021 methods exactly:
# - 600 HVGs (Seurat 2.3 FindVariable equivalent)
# - 20 PCs
# - No batch correction
# - Harmony used only if needed (paper mentions it as available fallback)

rule paper_normalize_hvg:
    input:
        h5ad = "data/processed/qc_filtered.h5ad"
    output:
        h5ad    = "data/processed/paper_track/normalized.h5ad",
        figures = directory("results/figures/paper_track/hvg_pca")
    params:
        n_hvgs  = config["paper_track"]["n_hvgs"],
        n_pcs   = config["paper_track"]["n_pcs"],
        flavor  = config["paper_track"]["hvg_flavor"],
    threads: 8
    resources:
        mem_mb          = 64000,
        runtime         = 30,
        slurm_partition = "campus"
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
            --hvg_flavor  {params.flavor} \
            2>&1 | tee {log}
        """


rule paper_umap_embed:
    input:
        h5ad = "data/processed/paper_track/normalized.h5ad"
    output:
        h5ad    = "data/processed/paper_track/embedded.h5ad",
        figures = directory("results/figures/paper_track/umap")
    params:
        n_neighbors = config["paper_track"]["n_neighbors"],
        n_pcs       = config["paper_track"]["n_pcs"],
        use_rep     = "X_pca"   # no batch correction — use PCA directly
    threads: 8
    resources:
        mem_mb          = 64000,
        runtime         = 30,
        slurm_partition = "campus"
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
            2>&1 | tee {log}
        """
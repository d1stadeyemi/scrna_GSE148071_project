# ── Track B: Alternative implimenting modern best-practice pipeline (scVI batch correction) ────────────
# Methodological upgrade over paper:
# - 3000 HVGs (seurat_v3)
# - scVI latent space (30 dims)
# - Batch correction across 42 patients

rule alternative_normalize_hvg:
    input:
        h5ad = "data/processed/qc_filtered.h5ad"
    output:
        h5ad    = "data/processed/alternative_track/normalized.h5ad",
        figures = directory("results/figures/alternative_track/hvg_pca")
    params:
        n_hvgs  = config["alternative_track"]["n_hvgs"],
        n_pcs   = 50,
        flavor  = config["alternative_track"]["hvg_flavor"],
    threads: 8
    resources:
        mem_mb          = 64000,
        runtime         = 30,
        slurm_partition = "campus"
    conda:
        "../envs/scvi.yaml"
    log:
        "logs/alternative_track/normalize_hvg.log"
    shell:
        """
        mkdir -p logs/alternative_track
        python scripts/normalize_hvg.py \
            --input       {input.h5ad} \
            --output      {output.h5ad} \
            --figures_dir {output.figures} \
            --n_hvgs      {params.n_hvgs} \
            --n_pcs       {params.n_pcs} \
            --hvg_flavor  {params.flavor} \
            2>&1 | tee {log}
        """


rule alternative_scvi_integrate:
    input:
        h5ad = "data/processed/alternative_track/normalized.h5ad"
    output:
        h5ad  = "data/processed/alternative_track/integrated.h5ad",
        model = directory("models/scvi_model")
    params:
        n_latent   = config["alternative_track"]["n_latent"],
        n_epochs   = config["alternative_track"]["n_epochs"],
        batch_size = config["alternative_track"]["batch_size"],
    threads: 16
    resources:
        mem_mb          = 64000,
        runtime         = 480,
        slurm_partition = "campus"
    conda:
        "../envs/scvi.yaml"
    log:
        "logs/alternative_track/scvi_integrate.log"
    shell:
        """
        python scripts/scvi_integrate.py \
            --input      {input.h5ad} \
            --output     {output.h5ad} \
            --model_dir  {output.model} \
            --n_latent   {params.n_latent} \
            --n_epochs   {params.n_epochs} \
            --batch_size {params.batch_size} \
            2>&1 | tee {log}
        """


rule alternative_umap_embed:
    input:
        h5ad = "data/processed/alternative_track/integrated.h5ad"
    output:
        h5ad    = "data/processed/alternative_track/embedded.h5ad",
        figures = directory("results/figures/alternative_track/umap")
    params:
        n_neighbors = config["alternative_track"]["n_neighbors"],
        n_pcs       = config["alternative_track"]["n_latent"],
        use_rep     = "X_scvi"
    threads: 8
    resources:
        mem_mb          = 64000,
        runtime         = 60,
        slurm_partition = "campus"
    conda:
        "../envs/scvi.yaml"
    log:
        "logs/alternative_track/umap_embed.log"
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
# =============================================================================
# Alternative Track — Modern best practice with scVI batch correction
# =============================================================================
# scVI integration, 3000 HVGs, 50 PCs, scVI latent UMAP
#
# Biological rationale:
#   scVI uses a VAE (Variational Autoencoder) to model gene expression
#   as a function of both biology (latent variables) and batch (patient ID).
#   This removes patient-level technical variation while theoretically
#   preserving biological variation.
#
#   Trade-off: May also remove BIOLOGICAL patient variation (cancer heterogeneity).
#   This track exists to quantify that trade-off empirically.
# =============================================================================


rule alternative_normalize_hvg:
    """
    Normalization and HVG selection for alternative track.
    
    Uses seurat_v3 flavor (modern standard):
    - Stabilizes variance across expression levels
    - Accounts for batch effects during HVG selection (batch_key=patient_id)
    - 3000 HVGs captures more subtle within-cell-type variation
    - Requires raw counts in layer (stored before normalization)
    """
    input:
        h5ad = "data/processed/qc_filtered.h5ad"
    output:
        h5ad    = "data/processed/alternative_track/normalized.h5ad",
        figures = directory("results/figures/alternative_track/hvg_pca")
    params:
        n_hvgs      = config["alternative_track"]["n_hvgs"],
        n_pcs       = config["alternative_track"]["n_pcs"],
        hvg_flavor  = config["alternative_track"]["hvg_flavor"]
    threads: 8
    resources:
        mem_mb          = 64000,
        runtime         = 30,
        slurm_partition = "campus-bigmem"
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
            --hvg_flavor  {params.hvg_flavor} \
            --track       alternative \
            2>&1 | tee {log}
        """


rule alternative_scvi_integrate:
    """
    scVI integration for alternative track.
    
    Model: VAE with:
    - n_latent=30 dimensions (standard for scRNA-seq)
    - n_layers=2, n_hidden=128 (paper-matching architecture)
    - batch_key=patient_id (42 patients as batches)
    - Early stopping (patience=10) to prevent over-training
    
    If model already exists at model_dir, loads instead of retraining.
    This allows the post-training steps to rerun without repeating GPU training.
    """
    input:
        h5ad = "data/processed/alternative_track/normalized.h5ad"
    output:
        h5ad  = "data/processed/alternative_track/integrated.h5ad",
        model = directory("models/scvi_model")
    params:
        n_latent   = config["alternative_track"]["n_latent"],
        n_epochs   = config["alternative_track"]["n_epochs"],
        batch_size = config["alternative_track"]["batch_size"]
    threads: 16
    resources:
        mem_mb          = 64000,
        runtime         = 480,
        slurm_partition = "campus-bigmem"
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
    """
    UMAP embedding for alternative track — uses scVI latent space.
    
    Using X_scvi (batch corrected) means:
    - All cell types harmonized across patients
    - Patient-specific cancer signals may be reduced
    - Better for identifying shared immune cell states
    Comparing with paper track UMAP will quantify batch correction impact.
    """
    input:
        h5ad = "data/processed/alternative_track/integrated.h5ad"
    output:
        h5ad    = "data/processed/alternative_track/embedded.h5ad",
        figures = directory("results/figures/alternative_track/umap")
    params:
        n_neighbors = config["alternative_track"]["n_neighbors"],
        n_pcs       = config["alternative_track"]["n_latent"],
        use_rep     = config["alternative_track"]["use_rep"]
    threads: 8
    resources:
        mem_mb          = 64000,
        runtime         = 60,
        slurm_partition = "campus-bigmem"
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
            --track       alternative \
            2>&1 | tee {log}
        """


rule alternative_cluster_annotate:
    """
    Leiden clustering and cell type annotation for alternative track.
    Same biological annotation logic as paper track.
    """
    input:
        h5ad = "data/processed/alternative_track/embedded.h5ad"
    output:
        h5ad    = "data/processed/alternative_track/annotated.h5ad",
        figures = directory("results/figures/alternative_track/annotation")
    params:
        tables_dir      = "results/tables/alternative_track",
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
        "logs/alternative_track/cluster_annotate.log"
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
            --track           alternative \
            2>&1 | tee {log}
        """

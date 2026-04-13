rule normalize_hvg:
    input:
        h5ad = "data/processed/qc_filtered.h5ad"
    output:
        h5ad    = "data/processed/normalized.h5ad",
        figures = directory("results/figures/hvg_pca")
    params:
        n_hvgs = config["dimred"]["n_hvgs"],
        n_pcs  = config["dimred"]["n_pcs"],
    threads: 8
    resources:
        mem_mb  = 64000,
        runtime = 30
    conda:
        "../envs/scvi.yaml"
    log:
        "logs/normalize_hvg.log"
    shell:
        """
        python scripts/normalize_hvg.py \
            --input       {input.h5ad} \
            --output      {output.h5ad} \
            --figures_dir {output.figures} \
            --n_hvgs      {params.n_hvgs} \
            --n_pcs       {params.n_pcs} \
            2>&1 | tee {log}
        """


rule scvi_integrate:
    input:
        h5ad = "data/processed/normalized.h5ad"
    output:
        h5ad  = "data/processed/scvi_integrated.h5ad",
        model = directory("models/scvi_model")
    params:
        n_latent   = config["dimred"]["n_latent"],
        n_epochs   = config["dimred"]["n_epochs"],
        batch_size = config["dimred"]["batch_size"],
    threads: 8
    resources:
        mem_mb         = 64000,
        runtime        = 360,
    conda:
        "../envs/scvi.yaml"
    log:
        "logs/scvi_integrate.log"
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


rule umap_embed:
    input:
        h5ad = "data/processed/scvi_integrated.h5ad"
    output:
        h5ad    = "data/processed/embedded.h5ad",
        figures = directory("results/figures/umap")
    params:
        n_neighbors = config["dimred"]["n_neighbors"],
        n_pcs       = config["dimred"]["n_latent"],
    threads: 8
    resources:
        mem_mb  = 64000,
        runtime = 60
    conda:
        "../envs/scvi.yaml"
    log:
        "logs/umap_embed.log"
    shell:
        """
        python scripts/umap_embed.py \
            --input       {input.h5ad} \
            --output      {output.h5ad} \
            --figures_dir {output.figures} \
            --n_neighbors {params.n_neighbors} \
            --n_pcs       {params.n_pcs} \
            2>&1 | tee {log}
        """

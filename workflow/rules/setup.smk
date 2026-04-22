# ── Setup: Create all required directories ───────────────────────────────────

rule setup_dirs:
    output:
        touch("logs/.setup_done")
    shell:
        """
        mkdir -p data/raw data/processed/paper_track data/processed/alternative_track
        mkdir -p data/metadata
        mkdir -p logs/paper_track logs/alternative_track
        mkdir -p models/scvi_model
        mkdir -p results/figures/qc
        mkdir -p results/figures/paper_track/{{hvg_pca,umap,annotation}}
        mkdir -p results/figures/alternative_track/{{hvg_pca,umap,annotation}}
        mkdir -p results/figures/comparison
        mkdir -p results/tables/paper_track results/tables/alternative_track
        mkdir -p notebooks
        touch {output}
        """

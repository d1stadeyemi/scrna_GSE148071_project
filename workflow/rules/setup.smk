# ── Setup: Create all required directories ───────────────────────────────────

rule setup_dirs:
    output:
        touch("logs/.setup_done")
    localrule: True         
    shell:
        """
        mkdir -p data/raw \
                 data/processed/paper_track \
                 data/processed/alternative_track \
                 data/metadata \
                 logs/paper_track \
                 logs/alternative_track \
                 models/scvi_model \
                 results/figures/qc \
                 results/figures/paper_track/hvg_pca \
                 results/figures/paper_track/umap \
                 results/figures/paper_track/annotation \
                 results/figures/alternative_track/hvg_pca \
                 results/figures/alternative_track/umap \
                 results/figures/alternative_track/annotation \
                 results/figures/comparison \
                 results/tables/paper_track \
                 results/tables/alternative_track \
                 notebooks
        """
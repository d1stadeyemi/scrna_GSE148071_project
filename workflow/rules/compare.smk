# ── Compare: Generate side-by-side track comparison report ───────────────────

rule compare_tracks:
    """
    Generate a comparison report between paper track and alternative track.
    
    Compares:
    - UMAP topology (patient mixing vs separation)
    - Cell type proportions per patient
    - Cluster count and resolution
    - Cancer cell patient specificity
    
    This comparison quantifies the biological impact of batch correction
    and justifies the use of paper track as primary analysis.
    """
    input:
        paper_h5ad = "data/processed/paper_track/annotated.h5ad",
        alt_h5ad   = "data/processed/alternative_track/annotated.h5ad"
    output:
        report = "results/figures/comparison/track_comparison.pdf"
    params:
        figures_dir = "results/figures/comparison"
    threads: 4
    resources:
        mem_mb          = 64000,
        runtime         = 30,
        slurm_partition = "campus-bigmem"
    conda:
        "../envs/scvi.yaml"
    log:
        "logs/compare_tracks.log"
    shell:
        """
        python scripts/compare_tracks.py \
            --paper_h5ad  {input.paper_h5ad} \
            --alt_h5ad    {input.alt_h5ad} \
            --figures_dir {params.figures_dir} \
            --output      {output.report} \
            2>&1 | tee {log}
        """

# ── QC: Cell and gene quality filtering ──────────────────────────────────────

rule qc_filter:
    """
    Apply QC filters following Wu et al. 2021 Methods:
      - min_genes >= 200   (remove empty droplets and dead cells)
      - max_genes <= 5000  (remove likely multiplets/doublets)
      - max_pct_mt <= 30%  (remove apoptotic cells)
      - min_cells >= 3     (remove very rare genes)
    
    Generates QC violin plots and scatter plots before and after filtering.
    Expected output: ~89,819 cells (paper: 90,406 — 0.6% difference due to
    minor GEO data updates since original publication).
    """
    input:
        h5ad = "data/processed/merged_raw.h5ad"
    output:
        h5ad    = "data/processed/qc_filtered.h5ad",
        figures = directory("results/figures/qc")
    params:
        min_genes  = config["qc"]["min_genes"],
        max_genes  = config["qc"]["max_genes"],
        max_pct_mt = config["qc"]["max_pct_mt"],
        min_cells  = config["qc"]["min_cells"]
    threads: 4
    resources:
        mem_mb          = 64000,
        runtime         = 30,
        slurm_partition = "campus"
    conda:
        "../envs/scvi.yaml"
    log:
        "logs/qc_filter.log"
    shell:
        """
        python scripts/qc_filter.py \
            --input      {input.h5ad} \
            --output     {output.h5ad} \
            --figures_dir {output.figures} \
            --min_genes  {params.min_genes} \
            --max_genes  {params.max_genes} \
            --max_pct_mt {params.max_pct_mt} \
            --min_cells  {params.min_cells} \
            2>&1 | tee {log}
        """

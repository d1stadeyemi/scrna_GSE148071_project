rule qc_filter:
    input:
        h5ad = "data/processed/merged_raw.h5ad"
    output:
        h5ad    = "data/processed/qc_filtered.h5ad",
        figures = directory("results/figures/qc")
    params:
        min_genes  = config["qc"]["min_genes"],
        max_genes  = config["qc"]["max_genes"],
        max_pct_mt = config["qc"]["max_pct_mt"],
        min_cells  = config["qc"]["min_cells"],
    threads: 4
    resources:
        mem_mb  = 32000,
        runtime = 30
    conda:
        "../envs/scanpy.yaml"
    log:
        "logs/qc_filter.log"
    shell:
        """
        python scripts/qc_filter.py \
            --input       {input.h5ad} \
            --output      {output.h5ad} \
            --figures_dir {output.figures} \
            --min_genes   {params.min_genes} \
            --max_genes   {params.max_genes} \
            --max_pct_mt  {params.max_pct_mt} \
            --min_cells   {params.min_cells} \
            2>&1 | tee {log}
        """
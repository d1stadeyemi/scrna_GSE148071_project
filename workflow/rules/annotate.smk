rule cluster_annotate:
    input:
        h5ad = "data/processed/{track}/embedded.h5ad"
    output:
        h5ad    = "data/processed/{track}/annotated.h5ad",
        figures = directory("results/figures/{track}/annotation")
    params:
        resolution = config["annotation"]["resolution"],
        tables_dir = "results/tables"
    threads: 8
    resources:
        mem_mb          = 64000,
        runtime         = 120,
    conda:
        "../envs/scvi.yaml"
    log:
        "logs/{track}/cluster_annotate.log"
    shell:
        """
        python scripts/cluster_annotate.py \
            --input       {input.h5ad} \
            --output      {output.h5ad} \
            --figures_dir {output.figures} \
            --tables_dir  {params.tables_dir} \
            --resolution  {params.resolution} \
            2>&1 | tee {log}
        """
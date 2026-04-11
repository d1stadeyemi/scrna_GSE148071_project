rule parse_to_anndata:
    input:
        inventory = "data/raw/inventory.txt"
    output:
        h5ad = "data/processed/merged_raw.h5ad"
    threads: 4
    resources:
        mem_mb  = 32000,    
        runtime = 60        
    conda:
        "../envs/scanpy.yaml"
    log:
        "logs/parse_to_anndata.log"
    shell:
        """
        python scripts/parse_to_anndata.py \
            --extracted_dir data/raw/extracted \
            --output {output.h5ad} \
            2>&1 | tee {log}
        """
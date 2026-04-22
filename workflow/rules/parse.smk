# ── Parse: Convert per-sample matrices to unified AnnData ────────────────────

rule parse_to_anndata:
    """
    Parse 42 per-sample expression matrices from GEO into a single AnnData.
    
    Each file: GSM4453576_P1_exp.txt.gz → P1 patient expression matrix
    Format: genes × cells (transposed from standard)
    Barcodes prefixed with sample-specific ID for global uniqueness.
    
    Expected output: ~90,406 cells × 29,528 genes (paper reports 90,406)
    """
    input:
        matrices = "data/raw/matrices",
        metadata = "data/metadata/geo_metadata.tsv"
    output:
        h5ad = "data/processed/merged_raw.h5ad"
    resources:
        mem_mb          = 64000,
        runtime         = 60,
        slurm_partition = "campus"
    conda:
        "../envs/scvi.yaml"
    log:
        "logs/parse_to_anndata.log"
    shell:
        """
        python scripts/parse_to_anndata.py \
            --input_dir  {input.matrices} \
            --metadata   {input.metadata} \
            --output     {output.h5ad} \
            2>&1 | tee {log}
        """

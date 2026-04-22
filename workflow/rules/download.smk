# ── Download: Fetch raw data from GEO ────────────────────────────────────────
# GSE148071 — Wu et al. 2021 NSCLC scRNA-seq
# 42 samples, HiSeq X Ten, 10x-like Singleron GEXSCOPE protocol

rule download_raw_tar:
    input:
        "logs/.setup_done"
    output:
        temp("data/raw/GSE148071_RAW.tar")
    params:
        url = config["geo"]["base_url"] + config["files"]["raw_tar"]
    resources:
        mem_mb          = 4000,
        runtime         = 120,
        slurm_partition = "campus"
    log:
        "logs/download_raw_tar.log"
    shell:
        """
        wget -q --show-progress \
            -O {output} \
            {params.url} \
            2>&1 | tee {log}
        echo "Download complete: $(du -sh {output} | cut -f1)" >> {log}
        """


rule extract_raw_tar:
    input:
        "data/raw/GSE148071_RAW.tar"
    output:
        directory("data/raw/matrices")
    resources:
        mem_mb          = 4000,
        runtime         = 30,
        slurm_partition = "campus"
    log:
        "logs/extract_raw_tar.log"
    shell:
        """
        mkdir -p {output}
        tar -xf {input} -C {output} 2>&1 | tee {log}
        echo "Extracted $(ls {output}/*.gz | wc -l) files" >> {log}
        """


rule download_series_matrix:
    input:
        "logs/.setup_done"
    output:
        "data/raw/GSE148071_series_matrix.txt.gz"
    params:
        url = config["geo"]["matrix_url"]
    resources:
        mem_mb          = 2000,
        runtime         = 20,
        slurm_partition = "campus"
    log:
        "logs/download_series_matrix.log"
    shell:
        """
        wget -q --show-progress \
            -O {output} \
            {params.url} \
            2>&1 | tee {log}
        """


rule parse_series_matrix:
    """
    Extract patient metadata from GEO series matrix.
    GEO only provides age and gender for this dataset.
    Clinical metadata (histology, mutation, smoking) is in Supplementary Table 1
    of the paper (PDF only, not machine-readable from GEO).
    """
    input:
        "data/raw/GSE148071_series_matrix.txt.gz"
    output:
        "data/metadata/geo_metadata.tsv"
    resources:
        mem_mb          = 4000,
        runtime         = 10,
        slurm_partition = "campus"
    log:
        "logs/parse_series_matrix.log"
    shell:
        """
        python scripts/parse_metadata.py \
            --input  {input} \
            --output {output} \
            2>&1 | tee {log}
        """

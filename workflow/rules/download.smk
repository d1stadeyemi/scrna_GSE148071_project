rule download_raw_tar:
    input:
        "logs/.setup_done"         
    output:
        f"{config['paths']['raw_dir']}/GSE148071_RAW.tar"
    params:
        url = config["geo"]["base_url"] + config["files"]["raw_tar"]
    threads: 2
    resources: 
        mem_mb = 4000,
        runtime = 60
    log:
        "logs/download_raw_tar.log"
    shell:
        """
        wget -q --show-progress -O {output} {params.url} > {log} 2>&1
        """

rule extract_raw_tar:
    input:
        "data/raw/GSE148071_RAW.tar"
    output:
        directory("data/raw/extracted/")
    log:
        "logs/extract_raw_tar.log"
    shell:
        """
        mkdir -p data/raw/extracted
        tar -xvf {input} -C data/raw/extracted 2> {log}
        """

rule inventory_files:
    input:
        "data/raw/extracted/"
    output:
        "data/raw/inventory.txt"
    log:
        "logs/inventory_files.log"
    shell:
        """
        find data/raw/extracted -type f | sort > {output} 2> {log}
        """
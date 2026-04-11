rule setup:
    output:
        touch("logs/.setup_done")  
    threads: 1
    resources:
        mem_mb = 1000,
        runtime = 10
    shell:
        """
        mkdir -p config \
                 workflow/rules workflow/envs \
                 scripts resources notebooks \
                 data/raw data/processed data/metadata \
                 results/figures results/tables \
                 logs
        touch data/raw/.gitkeep
        touch data/processed/.gitkeep
        touch results/figures/.gitkeep
        touch results/tables/.gitkeep
        echo "Directory scaffold complete"
        """
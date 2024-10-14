rule seurat:
    input: 
        mtx=rules.counts.output.results,
    output:
        rds="results/seurat/{sample}/{sample}.rds",
    params:
        url=config["get_cellranger"]["url"],
    log:
        "results/logs/seurat/{sample}.log",
    benchmark:
        "results/benchmarks/seurat/{sample}.txt",
    conda:
        "/public/home/weiyifan/miniforge3/envs/seurat4",
    script:
        "scripts/basic_seurat_new.R"

rule integration:
    input: 
        rds=[f"results/seurat/{sample}/{sample}.rds" for sample in SAMPLES],
    output:
        rds="results/integration/integrated.rds",
    params:
        url=config["get_cellranger"]["url"],
    log:
        "results/logs/integration/integration.log",
    benchmark:
        "results/benchmarks/integration/integration.txt",
    conda:
        "/public/home/weiyifan/miniforge3/envs/seurat4",
    script:
        "scripts/multi_Seurat.R"
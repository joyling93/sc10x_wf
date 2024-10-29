rule seurat:
    input: 
        mtx=rules.counts.output.mtx,
    output:
        rds="results/seurat/{sample}/{sample}_seurat.rds",
    log:
        "results/logs/seurat/{sample}.log",
    benchmark:
        "results/benchmarks/seurat/{sample}.txt",
    conda:
        "/public/home/weiyifan/miniforge3/envs/seurat4",
    script:
        "../scripts/basic_seurat_new.R"

rule integration:
    input: 
        rds=[f"results/seurat/{sample}/{sample}_seurat.rds" for sample in SAMPLES],
    output:
        rds="results/integration/integrated.rds",
    log:
        "results/logs/integration/integration.log",
    benchmark:
        "results/benchmarks/integration/integration.txt",
    conda:
        "/public/home/weiyifan/miniforge3/envs/sc_integration",
    script:
        "../scripts/multi_Seurat.R"
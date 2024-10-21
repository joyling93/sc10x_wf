rule seurat:
    input: 
        mtx="results/counts/{sample}_cr/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
    output:
        rds="results/seurat/{sample}/{sample}_seurat.rds",
    params:
        url=config["get_cellranger"]["url"],
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
    params:
        url=config["get_cellranger"]["url"],
    log:
        "results/logs/integration/integration.log",
    benchmark:
        "results/benchmarks/integration/integration.txt",
    conda:
        "/public/home/weiyifan/miniforge3/envs/sc_integration",
    script:
        "../scripts/multi_Seurat.R"
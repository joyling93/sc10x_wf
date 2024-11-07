rule clone_typing:
    input: 
        [f"results/seurat/{sample}/{sample}_seurat.rds" for sample in SAMPLES],
    output:
        directory("results/ct"),
    log:
        "results/logs/ct.log",
    benchmark:
        "results/benchmarks/ct.txt",
    params:
        sample=SAMPLES,
        vdj_type=config[["vdj_type"]],
    conda:
        "scRepertoire",
    script:
        "../scripts/scRepertoire.R"
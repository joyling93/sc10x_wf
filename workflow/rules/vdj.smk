rule clone_typing:
    input: 
        [f"results/counts/{sample}_cr/per_sample_outs/{sample}/vdj_b/filtered_contig_annotations.csv" for sample in SAMPLES],
    output:
        directory("results/ct"),
    log:
        "results/logs/ct.log",
    benchmark:
        "results/benchmarks/ct.txt",
    params:
        sample=SAMPLES,
        vdj_type=config["vdj_type"],
    conda:
        "scRepertoire",
    script:
        "../scripts/scRepertoire.R"
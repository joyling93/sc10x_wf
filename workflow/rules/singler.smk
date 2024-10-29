rule singler:
    input: 
        rds="results/integration/integrated.rds",
    output:
        outdir=directory("results/singler/"),
        anf="results/singler/singler_results.csv",
        xls="results/singler/FindAllMarkers.xls"
    log:
        "results/logs/singler.log",
    benchmark:
        "results/benchmarks/singler.txt",
    conda:
        "/public/home/weiyifan/miniforge3/envs/singler",
    script:
        "../scripts/singler.R"

rule enrichment:
    input:
        "results/singler/FindAllMarkers.xls",
    output:
        directory("results/enrichment/singler/"),
    threads: 10,
    params:
        db=config["species"],
    conda:
        "/public/home/weiyifan/miniforge3/envs/enrichment"
    log:
        "logs/enrich/singler.log",
    script:
        "../scripts/enrich.R"
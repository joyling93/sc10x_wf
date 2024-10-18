rule singler:
    input: 
        rds="results/integration/integrated.rds",
    output:
        outdir="results/singler/",
        anf="results/singler/singler_results.csv",
        "results/singler/FindAllMarkers.xls"
    params:
        url=config["get_cellranger"]["url"],
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
        directory("results/enrichment/{contrast}"),
    params:
        contrast="singler",
    conda:
        "/public/home/weiyifan/miniforge3/envs/enrichment"
    log:
        "logs/enrich/{contrast}.log",
    script:
        "../scripts/enrich.R"
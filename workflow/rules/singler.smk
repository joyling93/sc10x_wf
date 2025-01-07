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
        "singler",
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
        "enrichment"
    log:
        "results/logs/enrich/enrich.log",
    script:
        "../scripts/enrich.R"
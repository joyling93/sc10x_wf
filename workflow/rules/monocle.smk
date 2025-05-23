rule monocle:
    input: 
        rds="results/integration/integrated.rds",
        meta=rules.cpdb.output[4],
    output:
        "results/monocle/monocle.rds",
        "results/monocle/cell_trajectory.png",
        "results/monocle/Pseudotime_cell_trajectory.png",
        "results/monocle/ordering_genes.png",
        "results/monocle/cell_Pseudotime.txt",
        "results/monocle/cell_trajectory_stat.png",
    log:
        "results/logs/monocle.log",
    benchmark:
        "results/benchmarks/monocle.txt",
    priority:
        1,
    conda:
        "monocle2",
    threads:
        11,
    script:
        "../scripts/monocle2.R"
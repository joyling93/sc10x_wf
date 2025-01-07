rule rmd_report:
    input: 
        directory("resources/template"),
        #integration=rules.integration.output,
        #singler=rules.singler.output,
        rules.enrichment.output,
        #rules.cpdb.output[0],
        rules.monocle.output,
    output:
        "results/report.html",
    benchmark:
        "results/benchmarks/report.txt",
    conda:
        "seurat4",
    shell:
        """
            Rscript --vanilla -e 'rmarkdown::render("resources/template/report.Rmd", output_file="../../results/report.html", quiet=FALSE)'
        """
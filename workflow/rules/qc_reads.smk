rule fastqc:
    input:
        "data/{sample}_S1_L00{lane}_{read}_001.fastq.gz",
    output:
        html="results/fastqc/{lane}_{sample}_{read}.html",
        zip="results/fastqc/{lane}_{sample}_{read}_fastqc.zip",
    log:
        "results/logs/fastqc/{lane}_{sample}_{read}.log",
    group:
        "qc",
    benchmark:
        "results/benchmarks/fastqc/{lane}_{sample}_{read}.txt"
    threads: 5
    wrapper:
        "0.79.0/bio/fastqc"


rule multiqc:
    input:
        **agg_fastqc(),
    output:
        html=report(
            "results/multiqc/multiqc.html",
            caption="../reports/multiqc.rst",
            category="QC Metrics for Sequencing Reads",
        ),
    log:
        "results/logs/multiqc/multiqc.log",
    benchmark:
        "results/benchmarks/fastqc/multiqc.txt"
    wrapper:
        "v3.5.3/bio/multiqc"
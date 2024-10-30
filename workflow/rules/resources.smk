# Downloads and extracts cell ranger
# The naming and rules ensure that any verion can be used
rule get_cellranger:
    output:
        cr="resources/cellranger",
    params:
        url=config["get_cellranger"][config["pipeline"]]["bin"],
    log:
        "results/logs/get_cellranger/get_cellranger.log",
    shell:
        """
           ln -s "{params.url}"  resources/cellranger
        """


rule get_reference:
    output:
        dir=directory("resources/genome"),
    params:
        url=config["get_reference"][config["species"]],
    log:
        "results/logs/get_reference/get_reference.log",
    shell:
        """
            ln -s "{params.url}" resources/genome
        """


# This currently expects names as I receive them from the core
rule clean_names:
    input:
        get_fastqs,
    output:
        "data/{sample}_S1_L00{lane}_{read}_001.fastq.gz",
    log:
        "results/logs/clean_names/{lane}_{sample}_{read}.log",
    benchmark:
        "results/benchmarks/clean_names/{lane}_{sample}_{read}.txt"
    shell:
        """
        mv {input} {output} &> {log}
        """

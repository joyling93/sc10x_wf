# No conda env used as cellranger cannot be installed in one
rule counts:
    input:
        unpack(get_sample_reads),
        bin=rules.get_cellranger.output.cr,
        genome=rules.get_reference.output.dir,
    output:
        mtx=temp("results/counts/{sample}_cr/"),
        html=report(
            "results/counts/{sample}_cr/outs/web_summary.html",
            caption="../reports/counts.rst",
            category="Cellranger Counts",
            subcategory="{sample}",
        ),
    params:
        introns=convert_introns(),
        mem=config["counts"]["mem"],
        sp_extra=convert_sp_extra,
    log:
        "results/logs/counts/{sample}.log",
    benchmark:
        "results/benchmarks/counts/{sample}.txt"
    threads: 16
    shell:
        """
        {input.bin} \
        count \
        --nosecondary \
        {params.introns} \
        --id {wildcards.sample} \
        --transcriptome {input.genome} \
        --fastqs data \
        --sample {wildcards.sample} \
        --localcores {threads} \
        --localmem {params.mem} \
        {params.sp_extra} \
        --output-dir results/counts/{wildcards.sample} &> {log} ; \
        cp -r results/counts/{wildcards.sample}/outs/ results/counts/{wildcards.sample}_cr/
        """
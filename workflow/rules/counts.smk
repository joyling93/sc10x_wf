# No conda env used as cellranger cannot be installed in one
if config["pipeline"]=="multi" :
    rule counts:
        input:
            #unpack(get_sample_reads),
            rules.multiqc.output,
            csv=get_csv,
            bin=rules.get_cellranger.output.cr,
            genome=rules.get_reference.output.dir,
        wildcard_constraints:
            sample=".+-gex",
        output:
            mtx=directory("results/counts/{sample}_cr/"),
        params:
            sp_extra=convert_sp_extra,
        log:
            "results/logs/counts/{sample}.log",
        benchmark:
            "results/benchmarks/counts/{sample}.txt"
        threads: 16
        shell:
            """
            {input.bin} \
            multi \
            --id={wildcards.sample} \
            --csv={input.csv}\
            --output-dir results/counts/{wildcards.sample} &> {log} ; \
            cp -r results/counts/{wildcards.sample}/outs/ results/counts/{wildcards.sample}_cr/
            """
else :
    rule counts:
        input:
            unpack(get_sample_reads),
            bin=rules.get_cellranger.output.cr,
            genome=rules.get_reference.output.dir,
        output:
            mtx=temp(directory("results/counts/{sample}_cr/")),
        params:
            introns=convert_introns(),
            mem=config["counts"]["mem"],
            sp_extra=convert_sp_extra,
        log:
            "results/logs/counts/{sample}.log",
        resources:
            mem_mb=100000
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
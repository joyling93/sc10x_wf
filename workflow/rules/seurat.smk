rule seurat:
    input: 
        mtx=rules.counts.output.results,
    output:
        rds="results/seurat/{sample}/{sample}.rds",
    params:
        url=config["get_cellranger"]["url"],
    log:
        "results/logs/seurat/{sample}.log"",
    benchmark:
        "results/benchmarks/seurat/seurat.txt",
    conda:
        "/public/home/weiyifan/miniforge3/envs/seurat4",
    shell:
        """
           Rscript /public/home/weiyifan/xzm/workshop/basic_qc/basic_seurat_new.R --indir {input.mtx} --name {wildcards.sample} --outdir results/seurat/{wildcards.sample}
        """

rule integration:
    input: 
        rds=";".join([f"results/seurat/{sample}/{sample}.rds" for sample in SAMPLES]),
    output:
        rds="results/integration/integrated.rds,
    params:
        url=config["get_cellranger"]["url"],
    log:
        "results/logs/integration/integration.log"",
    benchmark:
        "results/benchmarks/integration/integration.txt",
    conda:
        "/public/home/weiyifan/miniforge3/envs/seurat4",
    shell:
        """
           Rscript /public/home/weiyifan/xzm/workshop/integration/multi_Seurat.R --compare {input.rds} --species GRCh38 --outdir results/integration/
        """
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(reticulate)
obj <- readRDS(snakemake@input[[1]])
outdir <- dirname(snakemake@output[[1]])


if (!file.exists(file.path(outdir, "matrix.mtx"))) {
        Matrix::writeMM(obj@assays$RNA@counts, 
                file = file.path(outdir, "matrix.mtx"))
        write(x = rownames(obj@assays$RNA@counts), 
                file = file.path(outdir, "features.tsv"))
        write(x = colnames(obj@assays$RNA@counts), 
                file = file.path(outdir, "barcodes.tsv"))
}
obj@meta.data$Cell = rownames(obj@meta.data)
df = obj@meta.data[, c("Cell", "seurat_clusters")]
write.table(df, file = file.path(outdir, "meta.tsv"), sep = "\t", quote = F, row.names = F)

#use_condaenv("/public/home/xiezhuoming/miniforge3/envs/cpdb")
use_python('/public/home/xiezhuoming/miniforge3/envs/seurat4/bin/python')
cpdb<-import("cellphonedb.src.core.methods")
if(snakemake@config[["species"]]=="homo_sapiens"){
        res<-cpdb$cpdb_statistical_analysis_method$call(
        cpdb_file_path = "/public/home/xiezhuoming/xzm/code/workshop/cellphonedb/v2/v5.0.0/cellphonedb.zip",
        meta_file_path = "./results/cpdb/meta.tsv",
        counts_file_path = "./results/cpdb/",
        counts_data = 'hgnc_symbol',
        score_interactions = TRUE,
        threshold = 0.1,
        threads = as.integer(snakemake@threads),
        output_path = outdir
)
}else{
        res<-NULL
}

saveRDS(res, file = snakemake@output[[1]])
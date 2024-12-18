library(SingleR)
library(Seurat)
library(tidyverse)
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
        

use_condaenv("/public/home/weiyifan/miniforge3/envs/cpdb")
cpdb<-import("cellphonedb.src.core.methods")

res<-cpdb$cpdb_statistical_analysis_method$call(
        cpdb_file_path = db,
        meta_file_path = "./results/cpdb/meta.tsv",
        counts_file_path = "./results/cpdb/",
        counts_data = 'hgnc_symbol',
        score_interactions = TRUE,
        threshold = 0.1,
        output_path = outdir
)

saveRDS(res, file = snakemake@output[[1]])
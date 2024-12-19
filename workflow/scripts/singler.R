log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(SingleR)
library(Seurat)
library(tidyverse)
obj <- readRDS(snakemake@input[[1]])
outdir <- snakemake@output[[1]]
rd<-switch(snakemake@config[["species"]],
       homo_sapiens = "/repository/Test/scRNA_mj/SingleR_database/HumanPrimaryCellAtlas_hpca.se_human.RData",
       mus_musculus = "/repository/Test/scRNA_mj/SingleR_database/MouseRNAseqData.Rdata"
)
hpca.se <- get(load(rd))

#转化SingleCellExperiment对象
Seurat_Object_Diet <- DietSeurat(obj, graphs = "pca")
SCE <- as.SingleCellExperiment(Seurat_Object_Diet)

pred.hesc <- SingleR(SCE, hpca.se,  labels=hpca.se$label.main, clusters=obj$seurat_clusters)

#写出注释结果
write.csv(pred.hesc,file.path(outdir,'singler_results.csv'))

#写入注释结果
correct_celltype <- tibble(
        seurat_clusters=levels(obj$seurat_clusters),
        singler=pred.hesc$labels
)

obj@meta.data['singler'] <- left_join(obj[['seurat_clusters']],correct_celltype)$singler

DimPlot(obj, label=TRUE, group.by = 'singler',reduction = 'umap')


Idents(obj) <- 'singler'
obj.markers <- FindAllMarkers(obj, logfc.threshold=0.25, min.pct=0.1) 

write_tsv(obj.markers,file.path(outdir,'FindAllMarkers.xls'))

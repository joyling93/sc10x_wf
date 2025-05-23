log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

scdir<-'/public/home/xiezhuoming/xzm/code/workshop'
source('/public/home/xiezhuoming/xzm/code/workshop/utilis.R')
dims <- 15
minpct <- 0.1
species <- snakemake@params[["db"]]
logfc<-'0.25'
minumi <- 3
mito <- NULL
resolution <- NULL
suppressMessages({
    library(Seurat)
    library(tidyverse)
})


input_dir <- paste0(snakemake@input[[1]],"/filtered_feature_bc_matrix/")
if(snakemake@config[["pipeline"]]=="multi"){
    input_dir <- paste0(snakemake@input[[1]],
        "/outs/per_sample_outs/",snakemake@wildcards[['sample']],"/count/sample_filtered_feature_bc_matrix/")
}

id_map <- read_delim(paste0( input_dir,'features.tsv.gz'), delim="\t",
                    col_names = c('Ensembl', 'Symbol', 'Type')) %>% dplyr::select(-Type) %>%
    mutate(Symbol_uniq=make.unique(Symbol))
# }
outdir <- dirname(snakemake@output[[1]])
sp <- snakemake@wildcards[['sample']]
print(c(outdir, sp, input_dir))

if(snakemake@config[["pipeline"]]=="spaceranger"){
    #library(patchwork)
    input_dir <- snakemake@input[[1]]
    brain <- Load10X_Spatial(input_dir)
    plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
    plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
    setwd(outdir)
    ggsave(paste0(sp, '_features_VlnPlot.png'), plot1 + plot2, dpi=300,width=14)
    brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
    brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
    brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
    brain <- FindClusters(brain, verbose = FALSE)
    brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
    p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
    p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
    p1 + p2
    ggsave(paste0(sp, '_DimPlot.png'), p1 + p2, dpi=300,width=14)
    saveRDS(brain, file = paste0(sp, "_seurat.rds"))
    q(save = "no")
}


data.count <- Read10X(input_dir)


obj <- CreateSeuratObject(counts = data.count,min.cells = 3,min.features = 200)
setwd(outdir)
human_hemo_gene <- unlist(strsplit('HBA1 HBA2 HBB HBD HBE1 HBG1 HBG2 HBM HBQ1 HBZ', ' '))
mouse_hemo_gene <- unlist(strsplit('Hbb-bt Hbb-bs Hbb-bh2 Hbb-bh1 Hbb-y Hba-x Hba-a1 Hbq1b Hba-a2 Hbq1a', ' '))
hemo_gene <- switch(species, homo_sapiens=human_hemo_gene, mus_musculus=mouse_hemo_gene, NA)

obj[['percent.mito']] <- PercentageFeatureSet(object = obj, pattern = '^(MT|mt|Mt)-')
ggsave(paste0(sp, '_features_VlnPlot.png'), VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size=0), dpi=300)
# if (!is.na(hemo_gene[1])){
#     obj[['percent.hemo']] <- PercentageFeatureSet(object = obj, features=hemo_gene)
# }

write.table(round(do.call("cbind", tapply(obj$percent.mito, Idents(obj), quantile, probs=seq(0,1,0.05))), digits = 3),
            file=paste0(sp, '_mito_quantile.xls'), sep="\t", col.names=FALSE, quote=FALSE)

# if (!is.na(hemo_gene[1])){
#     write.table(round(do.call("cbind", tapply(obj$percent.hemo, Idents(obj), quantile, probs=seq(0,1,0.05))), digits = 3),
#             file=paste0(sp, '_hemo_quantile.xls'), sep="\t", col.names=FALSE, quote=FALSE)
# }

filter_cells_tsv <- tibble()
# UMI过滤
umi_threshold <- ifelse(is.null(minumi), 0, minumi)
cellNum <- dim(obj)[2]
filter_cells_tsv <-  filter_cells_tsv %>% bind_rows(obj[[]] %>% filter(nCount_RNA < umi_threshold))
obj <- subset(obj, subset = nCount_RNA >= umi_threshold)
umiFilter <- dim(obj)[2]

# mito阈值判断
mt.p <- pnorm(obj$percent.mito, mean = median(obj$percent.mito), sd = mad(obj$percent.mito), lower.tail = FALSE)
mt.lim <- min(obj$percent.mito[which(p.adjust(mt.p, method = "fdr") < 0.05)])
mito_threshold <- ifelse(is.null(mito), mt.lim, mito) 

# mito过滤
filter_cells_tsv <- filter_cells_tsv %>% bind_rows(obj[[]] %>% filter(percent.mito > mito_threshold))
obj <- subset(obj, subset = `percent.mito` <= mito_threshold)
mitoFilter <- dim(obj)[2]

obj <- subset(obj, subset = nFeature_RNA<5000)

filter_stat_tsv <- tibble(`Sample name`=sp,
                           cellNum=format(cellNum,  big.mark=','),
                           umiThreshold=format(umi_threshold, big.mark=','),
                           umiFilter=format(umiFilter, big.mark=','),
                           mitoThreshold=format(round(mito_threshold,2)),
                           mitoFilter=format(mitoFilter, big.mark=','))

write_tsv(filter_stat_tsv, paste0(sp, '_filter_stat.tsv'))
write_tsv(filter_cells_tsv, paste0(sp, '_filter_cells.tsv'))

ggsave(paste0(sp, '_features_VlnPlot_filter.png'), VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size=0), dpi=300)

if (!is.na(hemo_gene[1])){
    ggsave(paste0(sp, '_features_VlnPlot2.png'), VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3), dpi=300)
}

plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mito")
# if (!is.na(hemo_gene[1])){
#     plot3 <- FeatureScatter(obj, feature1 = "nCount_RNA")
#     ggsave(paste0(sp, '_FeatureScatter.png'), CombinePlots(plots = list(plot1, plot2, plot3)), dpi=300)
# }else{
#     ggsave(paste0(sp, '_FeatureScatter.png'), CombinePlots(plots = list(plot1, plot2)), dpi=300)
# }

write.table(do.call("cbind", tapply(obj$nFeature_RNA, Idents(obj), quantile, probs=seq(0,1,0.05))),
            file=paste0(sp, '_nFeature_quantile.xls'), sep="\t", col.names=FALSE, quote=FALSE)
write.table(do.call("cbind", tapply(obj$nCount_RNA, Idents(obj), quantile,probs=seq(0,1,0.05))),
            file=paste0(sp, '_nCount_quantile.xls'), sep="\t", col.names=FALSE, quote=FALSE)

dims <- 1:dims

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)

obj <- JackStraw(obj, num.replicate = 100)
obj <- ScoreJackStraw(obj, dims= dims)
p0 <- JackStrawPlot(obj, dims = dims)
ggsave(paste0(sp, '_JackStrawPlot.png'), p0, dpi=300)

p0 <- ElbowPlot(obj)
ggsave(paste0(sp, '_ElbowPlot.png'), p0, dpi=300)

obj <- FindNeighbors(obj, dims = dims)

if (is.null(resolution)){
    resolution <- seq(0.2, 1.4, 0.3)
    resolution_idents <- 'RNA_snn_res.0.8'
}else{
    resolution <- resolution
    resolution_idents <- paste0('RNA_snn_res.', resolution)
}

obj <- FindClusters(obj, resolution = resolution)

res.table  <- obj[[]] %>% summarize(across(starts_with("RNA_snn_res"), n_distinct)) %>% t() %>%
                  as.data.frame() %>% dplyr::rename(cluster_num=V1) %>%rownames_to_column('resolution')
write_tsv(res.table, paste0(sp, '_resolution.xls'))


cluster_num <- res.table$cluster_num
names(cluster_num) <- res.table$resolution
if (cluster_num[resolution_idents]>length(MYCOLOR)) {
    printf('too many cluster with %s, %s>%s', resolution_idents, cluster_num[resolution_idents], length(MYCOLOR))
    for (i in (match(resolution_idents,res.table$resolution)-1):1) {
        if (cluster_num[i] <=42) {
            resolution_idents <- names(cluster_num[i])
            printf('ok with %s, %s<%s', resolution_idents, cluster_num[resolution_idents], length(MYCOLOR))
            break
        }
    }
}


Idents(obj) <- resolution_idents
obj[['seurat_clusters']] <- obj[[resolution_idents]]

cellnum_df <- obj[[]] %>% count(!!sym(resolution_idents)) %>%
              dplyr::rename(cluster=!!sym(resolution_idents), cellnum=n) %>%
              mutate(percent=cellnum/sum(cellnum)*100)
write_tsv(cellnum_df, paste0(sp, '_cellnum.xls'))
#ggplot(cellnum_df, aes(x='sample', y=percent, fill=cluster)) + geom_bar(stat='identity') +  coord_polar(theta = 'y')
p0 <- ggplot(cellnum_df, aes(x=cluster, y=percent, fill=cluster)) + 
          geom_bar(stat='identity') + 
          xlab('cluster') + ylab('percent (%)') + 
          scale_fill_manual(values=MYCOLOR) +
          theme_classic() + theme(legend.title = element_blank())
ggsave(paste0(sp, '_cellnum.png'), p0, dpi=300)

obj <- RunUMAP(obj, reduction = "pca", dims = dims)
obj <- RunTSNE(obj, reduction = "pca", dims = dims, check_duplicates = FALSE)

p0 <- DimPlot(obj, reduction = 'umap', cols=MYCOLOR)
ggsave(paste0(sp, '_umap.png'), p0, dpi=300)

p0 <- DimPlot(obj, reduction = 'tsne', cols=MYCOLOR)
ggsave(paste0(sp, '_tsne.png'), p0, dpi=300)

pal <-colorRampPalette(c("blue","cyan", "yellow","red"))
reduction_loci <- as.data.frame(Embeddings(obj, reduction="tsne"))
reduction_loci <- cbind(reduction_loci, obj[[]])

p <- ggplot(reduction_loci, aes_string(colnames(reduction_loci)[1], colnames(reduction_loci)[2]))
p1 <- p + geom_point(aes(color=nCount_RNA)) + scale_colour_gradientn(colors=pal(500)) + theme_classic()
ggsave(paste0(sp, '_reduction_umi.png'), plot=p1, dpi=300)
write_tsv(reduction_loci %>% rownames_to_column('cellID'),
              paste0(sp, '_reduction_umi.xls'))

obj.markers <- FindAllMarkers(obj, logfc.threshold=logfc, min.pct=minpct) %>% 
                    left_join(id_map, by=c('gene'='Symbol_uniq')) %>% 
                    dplyr::select(-Symbol) %>% dplyr::relocate(gene, Ensembl)
#write.table(obj.markers, file=paste0(sp, '_FindAllMarkers.xls'), row.names=FALSE, sep="\t", quote=FALSE)
write_tsv(obj.markers, paste0(sp, '_FindAllMarkers.xls'))

top9 <- obj.markers %>% group_by(cluster) %>% top_n(n = 9, wt = avg_log2FC)
obj <- ScaleData(obj, features=top9$gene)
p2 <- refineDoHeatmap(obj, features = unique(top9$gene), label=FALSE)
ggsave(paste0(sp, '_top9_heatmap.png'), plot=p2, dpi=300)

l0 <- group_split(top9)
names(l0) <- group_keys(top9)[['cluster']]
# for (x in names(l0)){
#     tmp0 <- arrangeTop9(
#                 FeaturePlot(obj, features=l0[[x]][['gene']], reduction=tsne, pt.size=0.8, combine=FALSE)
#             )
#     ggsave(paste0(sp, '_cluster', x, '_top9_FeaturePlot.png'), plot=tmp0)
#     tmp1 <- arrangeTop9(
#                 VlnPlot(obj, features=l0[[x]][['gene']], pt.size=0, cols=MYCOLOR[1:length(levels(Idents(obj)))], combine=FALSE),
#                 legend.position='none'
#             )
#     ggsave(paste0(sp, '_cluster', x, '_top9_VlnPlot.png'), plot=tmp1, dpi=300)
# }

#预测多倍体
#obj <- doublets_pred(obj)

saveRDS(obj, file = paste0(sp, "_seurat.rds"))


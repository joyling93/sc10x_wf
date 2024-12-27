log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(Seurat)
library(dplyr)
library(stringr)
library(argparse)
library(reshape2)
library(tidyr)
library(ggplot2)
library(scales)
library(cowplot)
library(purrr)

source('/public/home/weiyifan/xzm/code/workshop/integration/utils_std.R')


parser = ArgumentParser()
parser$add_argument("--Nfeatures",help="nfeatures of FindVariableFeatures()", default='2000' )
parser$add_argument("--prefix", help="prefix of results",default="")
parser$add_argument("--resolution", help="resolution for cluster",default='0.6')
parser$add_argument("--dims", help="dims_use for cluster",default='30')
parser$add_argument("--min_cells", help="genes expressed in min cells",default='3')
parser$add_argument("--x_low_cutoff", help="x_low_cutoff for search high variable genes",default='0.0125')
parser$add_argument("--x_high_cutoff", help="x_high_cutoff for search high variable genes",default='3')
parser$add_argument("--y_cutoff", help="y_cutoff for search high variable genes",default='0.5')
args <- parser$parse_args()


str(args)

#path=args$path
#compare=args$compare
####gene_path=args$gene_path
species = "GRCh38"
outdir = "results/integration/"
prefix=args$prefix
resolution=args$resolution
dims=args$dims
min_cells=args$min_cells
x_low_cutoff=args$x_low_cutoff
x_high_cutoff=args$x_high_cutoff
y_cutoff=args$y_cutoff
Nfeatures = args$Nfeatures

resolution <- as.numeric(resolution)
min_cells <- as.numeric(min_cells)
Nfeatures<-as.numeric(Nfeatures)
dims<- as.numeric(dims)
x_low_cutoff <- as.numeric(x_low_cutoff)
x_high_cutoff <- as.numeric(x_high_cutoff)
y_cutoff <- as.numeric(y_cutoff)

if (is.null(prefix)){
	prefix="inte"
}

if (!dir.exists(paste0(outdir,'/Anchors'))){
  dir.create(paste0(outdir,'/Anchors'))
}


if (!dir.exists(paste0(outdir,'/DIFF'))){
  dir.create(paste0(outdir,'/DIFF'))
}

if (!dir.exists(paste0(outdir,'/Marker'))){
  dir.create(paste0(outdir,'/Marker'))
}

if (!dir.exists(paste0(outdir,'/CellsRatio'))){
  dir.create(paste0(outdir,'/CellsRatio'))
}

ob.list <- list()
samples<- list()

numsap=1

for (fp in snakemake@input[["rds"]]){
  	pbmc <- readRDS(fp)
    each <- str_split(basename(fp),'.rds')[[1]][1]
        if(length(grep('-1',colnames(pbmc@assays$RNA@counts)[1]))){
	colnames(pbmc@assays$RNA@counts) <- str_replace_all(colnames(pbmc@assays$RNA@counts), '-1',paste0('-',numsap))
        }else{
        colnames(pbmc@assays$RNA@counts) <- paste0(colnames(pbmc@assays$RNA@counts),'-',numsap)
        }
	ob <- CreateSeuratObject(counts =pbmc@assays$RNA@counts,project =each,min.cells = min_cells)
	ob <- subset(ob,subset = nCount_RNA>500 & nFeature_RNA > 200 & nFeature_RNA <6000)
	ob$stim <-each
	#对于不平行的数据容易报错
	#if(c('scrublet.pred')%in%colnames(pbmc@meta.data)){
	#        ob$scrublet.pred <- pbmc$scrublet.pred
	#        ob$scrublet.score <- pbmc$scrublet.score   
	#}
	ob <- NormalizeData(ob)
	ob <- FindVariableFeatures(ob,  selection.method = "vst",nfeatures = Nfeatures)
	numsap=numsap+1
	ob.list[[each]] <- ob
  samples[[each]] <- each
}


obj.anchors <- FindIntegrationAnchors(object.list = ob.list, dims = 1:dims)
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:dims)

DefaultAssay(obj.integrated) <- "integrated"

obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = dims, verbose = FALSE)
obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:dims)
obj.integrated <- RunTSNE(obj.integrated, reduction = "pca", dims = 1:dims)
obj.integrated <- FindNeighbors(obj.integrated, dims = 1:dims)
obj.integrated <- FindClusters(obj.integrated, resolution = resolution)

saveRDS(obj.integrated, file =paste0(outdir,"/","integrated.rds"))
saveRDS(obj.integrated@meta.data, file =paste0(outdir,"/",prefix,"_integrated_medt.rds"))

p1 <- DimPlot(obj.integrated, reduction = "umap", group.by = "stim")
p2 <- DimPlot(obj.integrated, reduction = "umap", label = TRUE, repel = TRUE)

p3 <- DimPlot(obj.integrated, reduction = "tsne", group.by = "stim")
p4 <- DimPlot(obj.integrated, reduction = "tsne", label = TRUE, repel = TRUE)

ggsave(paste0(outdir,'/Anchors/',prefix,'_clusters_UMAP.pdf'),p2)
ggsave(paste0(outdir,'/Anchors/',prefix,'_clusters_UMAP.png'),p2)
ggsave(paste0(outdir,'/Anchors/',prefix,'_sample_UMAP.pdf'),p1)
ggsave(paste0(outdir,'/Anchors/',prefix,'_sample_UMAP.png'),p1)

ggsave(paste0(outdir,'/Anchors/',prefix,'_clusters_tSNE.pdf'),p4)
ggsave(paste0(outdir,'/Anchors/',prefix,'_clusters_tSNE.png'),p4)
ggsave(paste0(outdir,'/Anchors/',prefix,'_sample_tSNE.pdf'),p3)
ggsave(paste0(outdir,'/Anchors/',prefix,'_sample_tSNE.png'),p3)


## 写出细胞纬度坐标
umap<-obj.integrated@reductions$umap@cell.embeddings
umap<-cbind(rownames(umap),umap)
colnames(umap)[1]<-'Barcode'
write.table(umap,file=paste0(outdir,'/Anchors/',prefix,'_UMAP.csv'),quote=F,sep=',',row.names=F,col.names=T)

tsne<-obj.integrated@reductions$tsne@cell.embeddings
tsne<-cbind(rownames(tsne),tsne)
colnames(tsne)[1]<-'Barcode'
write.table(tsne,file=paste0(outdir,'/Anchors/',prefix,'_tSNE.csv'),quote=F,sep=',',row.names=F,col.names=T)

### 写出clusster 信息
clus<-as.data.frame(Idents(obj.integrated))
clus<-cbind(rownames(clus),clus)
colnames(clus)<-c('Barcode','Cluster')
write.table(clus,file=paste0(outdir,'/Anchors/',prefix,'_cluster.csv'),quote=F,sep=',',row.names=F,col.names=T)


sample_number <- length(samples)
num <- sample_number + 1

b <- mutate(clus,Sample=as.numeric(str_split(clus$Barcode,'-',simplify = TRUE)[,2]))

c <- b %>%
  group_by(Cluster) %>%
  count(Sample) %>%
  spread(key=Sample,value=n)
colnames(c)[2:num] <- samples

#calculate the relative percent of each sample in each cluster
df <- c[,-1]
rowsum <- rowSums(df,na.rm=T)
df_r <- df/rowsum

#print out the relative table with heads
df_r <- cbind(as.numeric(rownames(df_r))-1, df_r)
colnames(df_r)[1] <- "Cluster"
colnames(df_r)[2:num] <- samples

write.table(c,file=paste0(outdir,'/CellsRatio/',prefix,'_cluster_abundance.csv'),sep=',',row.names = F, quote=F)
write.table(df_r,file=paste0(outdir,'/CellsRatio/',prefix,'_cluster_persent.csv'),sep=',',row.names = F, quote=F)


colour1 <-MYCOLOR[1:sample_number]
##colour1 <- hue_pal()(sample_number)

td <- gather(df_r,key="Cluster Name",value="Cells Ratio",-Cluster)
td[,1] <- factor(td[,1], levels = sort(as.numeric(df_r$Cluster)))
td[,2] <- factor(td[,2], levels = samples)


plt<- ggplot(td,aes(x=td[,1],y=td[,3],fill=td[,2]))+
  geom_bar(position = 'stack',stat="identity")+
  labs(x="Cluster Name",y="Cells Ratio")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),axis.text = element_text(color="black"))+
  scale_y_continuous(expand=c(0.001,0.001))+
  scale_fill_manual(values=colour1)+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Sample'))

pdf(paste0(outdir,'/CellsRatio/',prefix,'_cluster_persent.pdf'))
plt
dev.off()

png(paste0(outdir,'/CellsRatio/',prefix,'_cluster_persent.png'),type="cairo-png")
plt
dev.off()

#clusters in sample percent
df_c <- t(t(df)/rowSums(t(df),na.rm=T))
row.names(df_c) <-c(1:dim(df_c)[1])
df_c <- as.data.frame(cbind((as.numeric(rownames(df_c))-1), df_c))
colnames(df_c)[1] <- "Cluster"
write.table(df_c,file=paste0(outdir,'/CellsRatio/',prefix,'_sample_persent.csv'),sep=',',row.names = F, quote=F)

cluster_number <- length(row.names(df_c))

#colour2 <- hue_pal()(cluster_number)
colour2 <-MYCOLOR[1:cluster_number]


td_c <- gather(df_c,key="Sample Name",value="Cells Ratio",-Cluster)
td_c[,1] <- factor(td_c[,1], levels = sort(as.numeric(df_c$Cluster)))
td_c[,2] <- factor(td_c[,2], levels = samples)

plt_c<- ggplot(td_c,aes(x=td_c[,2],y=td_c[,3],fill=td_c[,1]))+
  geom_bar(position = 'stack',stat="identity")+
  labs(x="Sample Name",y="Cells Ratio")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),axis.text = element_text(color="black"))+
  scale_y_continuous(expand=c(0.001,0.001))+
  scale_fill_manual(values=colour2)+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Cluster'))

pdf(paste0(outdir,'/CellsRatio/',prefix,'_sample_persent.pdf'))
plt_c
dev.off()

png(paste0(outdir,'/CellsRatio/',prefix,'_sample_persent.png'),type="cairo-png")
plt_c
dev.off()


####Find All Markers

DefaultAssay(obj.integrated) <- "RNA"

obj.markers <- FindAllMarkers(obj.integrated, only.pos = F, min.pct = 0.25,logfc.threshold=0.25)

obj.markers <- obj.markers[,c(7,6,2,1,5,3,4)]
write.table(obj.markers,file=paste0(outdir,'/DIFF/','Cluster_diff.xls'),quote=F,row.names=F,col.names=T,sep='\t')

cluster <- unique(obj.markers$cluster)
for (i in cluster) {
  data <- filter(obj.markers,cluster == i)
  write.table(data,paste0(outdir,'/DIFF/','Cluster_',i,'_diff.xls'),quote=F,row.names=F,col.names=T,sep='\t')
  data1 <- filter(obj.markers,cluster == i, p_val_adj < 0.05, avg_log2FC > 0)
  write.table(data1,paste0(outdir,'/DIFF/','Cluster_',i,'_diff_significant.xls'),quote=F,row.names=F,col.names=T,sep='\t')
}

top9 <- obj.markers %>% group_by(cluster) %>% top_n(n = 9, wt = avg_log2FC)
obj.integrated <- ScaleData(obj.integrated, features=top9$gene, assay='RNA')

p2 <- refineDoHeatmap(obj.integrated, features = unique(top9$gene), label=FALSE, assay='RNA')

ggsave(paste0(outdir,'/Marker/integrated_top9_Heatmap.png'), plot=p2, dpi=300)

# l0 <- group_split(top9)
# names(l0) <- group_keys(top9)[['cluster']]
# for (x in names(l0)){
#     tmp0 <- arrangeTop9(
#                 FeaturePlot(obj.integrated, features=l0[[x]][['gene']], reduction='umap', pt.size=0.8, combine=FALSE)
#             )
#     tmp0 <- tmp0+labs(title = x)
#     ggsave(paste0(outdir,'/Marker/integrated_cluster', x, '_top9_FeaturePlot.png'), plot=tmp0)
#     tmp1 <- arrangeTop9(
#                 VlnPlot(obj.integrated, features=l0[[x]][['gene']], pt.size=0, cols=MYCOLOR[1:length(levels(Idents(obj.integrated)))], combine=FALSE),
#                 legend.position='none')
#     tmp1 <- tmp1+labs(title = x)
#     ggsave(paste0(outdir,'/Marker/integrated_cluster', x, '_top9_VlnPlot.png'), plot=tmp1, dpi=300)
# }

#saveRDS(obj.integrated, file =paste0(outdir,"/",prefix,"_integrated_seurat.rds"))


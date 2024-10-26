suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(org.Mm.eg.db))
#suppressMessages(library(org.Rn.eg.db))
#suppressMessages(library(org.Ss.eg.db))
suppressMessages(library(KEGG.db))
library(STRINGdb)
#读取差异表达基因
library(tidyverse)
gene_list=read.table(snakemake@input[[1]],header = T)
db<-"homo_sapiens"
outdir<-snakemake@output[[1]]
dir.create(outdir,recursive = T)

if(!db%in%c("homo_sapiens","org.Mm.eg.db")){
    db<-"homo_sapiens"
}
db<-
        switch(db,
                homo_sapiens=c('org.Hs.eg.db','Homo sapiens','hsa',"9606"),
                mus_musculus=c('org.Mm.eg.db','Mus musculus','mmu',"10090"),
                rno=c('org.Rn.eg.db','Rattus norvegicus','rno'),
                ss=c('org.Ss.eg.db','Sus scrofa','susScr')
        )

print(head(gene_list$gene))
eg <- bitr(gene_list$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=db[1])
gene_list <- gene_list%>%left_join(eg,by=join_by(gene==SYMBOL))%>%drop_na()
#geneList<-eg$ENTREZID

enrich_ora<- function(gl,db,out_dir,use_internal_data=F){
                dir.create(out_dir,recursive = T)
                #gl<-gl[[1]]
                saveRDS(gl,'test.rds')
                ##ppi
                string_db <- STRINGdb$new(version="12.0",species=as.numeric(db[4]),score_threshold=700,
                                input_directory= "/public/home/weiyifan/database/stringdb/")
                print(head(gl$gene))

                deg_mapped <- string_db$map(as.data.frame(gl), "gene", removeUnmappedRows = TRUE )
                
                cat("Total String id mapped :", dim(deg_mapped)[1])
                
                info <- string_db$get_interactions(deg_mapped$STRING_id)
                
                pdb<-string_db$get_proteins()
                info <- left_join(info,pdb,by = c('from'='protein_external_id'))%>%
                            left_join(pdb,by = c('to'='protein_external_id'))%>%
                            dplyr::select(4,7,3)%>%
                            rename('from'=preferred_name.x,'to'=preferred_name.y)
                write.table(info, file = file.path(out_dir,"STRING_info.txt"),sep="\t", row.names =F, quote = F)
                
                pl <- list()
                geneList <- gl$ENTREZID
                safe_enrichGO<-possibly(enrichGO,'no results')
                ego1 <- enrichGO(gene         = geneList,
                                 OrgDb         = db[1],
                                 keyType       = 'ENTREZID',
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.01,
                                 qvalueCutoff  = 0.05,
                                 readable=T)
                write.table(ego1,file=file.path(out_dir,'GO_bp_enrich_result.tsv'),quote=F,sep='\t',row.names=F,col.names=T)
                
                if(!is.null(ego1)){
                        p<-dotplot(ego1, showCategory=30) + ggtitle("BP for ORA")
                        if(dim(p$data)[1]>0){
                                pl<-append(pl,list(p))
                                ggsave(file.path(out_dir,'bp_ora.png'),p,width=20,height = 12)
                                ggsave(file.path(out_dir,'bp_ora.pdf'),p,width=20,height = 12)
                        }
                }
                print('bp done')
                
                ego2 <- enrichGO(gene         = geneList,
                                 OrgDb         = db[1],
                                 keyType       = 'ENTREZID',
                                 ont           = "MF",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.01,
                                 qvalueCutoff  = 0.05,
                                 readable=T)
                write.table(ego2,file=file.path(out_dir,'GO_mf_enrich_result.tsv'),quote=F,sep='\t',row.names=F,col.names=T)
                
                if(!is.null(ego2)){
                        p<-dotplot(ego2, showCategory=30) + ggtitle("MF for ORA")
                        if(dim(p$data)[1]>0){
                                pl<-append(pl,list(p))
                                ggsave(file.path(out_dir,'mf_ora.png'),p,width=20,height = 12)
                                ggsave(file.path(out_dir,'mf_ora.pdf'),p,width=20,height = 12)
                        }
                }
                print('mf done')
                
                ego5 <- enrichGO(gene         = geneList,
                                 OrgDb         = db[1],
                                 keyType       = 'ENTREZID',
                                 ont           = "CC",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.01,
                                 qvalueCutoff  = 0.05,
                                 readable=T)
                write.table(ego5,file=file.path(out_dir,'GO_cc_enrich_result.tsv'),quote=F,sep='\t',row.names=F,col.names=T)
                
                if(!is.null(ego5)){
                        p<-dotplot(ego5, showCategory=30) + ggtitle("CC for ORA")
                        if(dim(p$data)[1]>0){
                                pl<-append(pl,list(p))
                                ggsave(file.path(out_dir,'CC_ora.png'),p,width=20,height = 12)
                                ggsave(file.path(out_dir,'CC_ora.pdf'),p,width=20,height = 12)
                        }
                }
                print('cc done')
                
                safe_enrichKEGG<-possibly(enrichKEGG,'no results')
                ego4 <- enrichKEGG(gene         = geneList,
                                   organism     = db[3],
                                   pvalueCutoff = 0.05,
                                   use_internal_data = use_internal_data
                                   )
                
                if(!is.null(ego4)){
                        p<-dotplot(ego4, showCategory=30) + ggtitle("kegg for ORA")
                        if(dim(p$data)[1]>0){
                                pl<-append(pl,list(p))
                                ggsave(file.path(out_dir,'kegg_ora.png'),p,width=20,height = 12)
                                ggsave(file.path(out_dir,'kegg_ora.pdf'),p,width=20,height = 12)
                                ego4<-setReadable(ego4,OrgDb = db[1], keyType="ENTREZID")
                                write.table(ego4,file=file.path(out_dir,'kegg_enrich_result.tsv'),quote=F,sep='\t',row.names=F,col.names=T)
                                
                        }
                }
                print('kegg done')
                
                dl <- 
                list('anno_index'=gl,
                     'bp_ora'=ego1,
                     'mf_ora'=ego2,
                     'cc_ora'=ego5,
                     'kegg_ora'=ego4)
                #saveRDS(dl,file.path(out_dir,'enrich_list.rds'))
                openxlsx::write.xlsx(dl,file.path(out_dir,'enrich_list.xlsx'),overwrite =T)
                
                #Sys.sleep(30)
}

gene_list <- gene_list%>%
    mutate(type=ifelse(p_val_adj>0.05,'not_significant',
        ifelse(avg_log2FC>0,'up','down')))

##ora
gl <- gene_list%>%dplyr::filter(type!='not_significant')%>%
        group_by(cluster,type)%>%nest()
pwalk(list(gl$cluster,gl$type,gl$data),~enrich_ora(gl=..3,db=db,out_dir=file.path(outdir,..1,..2)))
#walk2(gl,cluster,type,~enrich_ora(gl=data,db=db,out_dir=file.path(outdir,.x,.y)))
##gsea
library(BiocParallel)
options(MulticoreParam=MulticoreParam(workers=snakemake@threads[[1]]))

gl2 <- gene_list%>%group_by(cluster)%>%nest()

enrich_gsea<-function(gl,db,out_dir){
        geneList <- gl$avg_log2FC
        names(geneList) <- gl$ENTREZID
        geneList <- sort(geneList, decreasing = TRUE)

        ego <- gseGO(geneList     = geneList,
                OrgDb        = db[1],
                ont          = "all",
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = FALSE
                )
        write.csv(ego,file.path(outdir,'go_gsea.csv'))
        saveRDS(ego,file.path(outdir,'go_gsea.rds'))

        kk <- gseKEGG(geneList     = geneList,
                organism     = db[3],
                minGSSize    = 120,
                pvalueCutoff = 0.05,
                verbose      = FALSE,
                use_internal_data = F
                )

        write.csv(kk,file.path(outdir,'kegg_gsea.csv'))
        saveRDS(kk,file.path(outdir,'kegg_gsea.rds'))
}
walk2(gl2$cluster,gl2$data,~enrich_gsea(gl=..2,db=db,out_dir=file.path(outdir,..1)))
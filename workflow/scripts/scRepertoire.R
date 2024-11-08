log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

suppressMessages(library(scRepertoire))

#input <- paste0(snakemake@input[[1]],"per_sample_outs/",{snakemake@wildcards[["sample"]],"/vdj_b/filtered_contig_annotations.csv")
input <- snakemake@input
print(input)
contig_list <- 
purrr::map(input,function(sample){
        read.csv(sample, stringsAsFactors = F)
})
#分别有BCR和TCR，为避免cellbarcode重复加上 samples_ID 前缀，
#本质是一个以samples_ID分组的list,把同一细胞的等位基因整合到一起
# combined <- combineTCR(contig_list, 
#                        samples = sample, 
#                        ID = snakemake@params[["sample"]], cells =snakemake@params[["vdj_type"]])#T-GD for gamma-delta TCR
print(as.character(snakemake@params[["sample"]]))
#str(contig_list)
combined <- combineBCR(contig_list, 
                       samples = as.character(snakemake@params[["sample"]]),
                        threshold = 0.85)

print(snakemake@output[[1]])
dir.create(snakemake@output[[1]])
setwd(snakemake@output[[1]])
#分组计算 unique clonetype 
quantContig(combined, cloneCall="gene+nt", scale = TRUE)
ggsave('unique_clonetype_per_sample.png', width = 10, height = 10)

#计算clonetype数目分布
abundanceContig(combined, cloneCall = "gene", scale = FALSE)
ggsave('number_of_clonetypes_per_sample.png', width = 10, height = 10)

#计算cdr3长度分布
lengthContig(combined, cloneCall="aa") 
ggsave('cdr3_distribution.png', width = 10, height = 10)
#Clonal Space Homeostasis
clonalHomeostasis(combined, cloneCall = "gene")
ggsave('Homeostasis.png', width = 10, height = 10)
clonalProportion(combined, cloneCall = "gene") 
ggsave('Homeostasis_numbin.png', width = 10, height = 10)
#clonetype 种类相似性
clonalOverlap(combined, cloneCall = "gene+nt", method = "morisita")
ggsave("clonalOverlap.png", width = 10, height = 10)
#多样性分析
clonalDiversity(combined, cloneCall = "gene")
ggsave("Diversity.png", width = 10, height = 10)

#
saveRDS(combined, 'combined.rds')
rm(list = ls())
p_load(maftools, GenVisR, ggplot2, reshape2, data.table, ggpubr, export)


cln_survival <- read.csv("outdata/cln_xena.csv")
allid <- read.csv("outdata/allid_493.txt") %>% .[,2]
load("outdata/sample_group.Rdata")
load("outdata/4groups_mafs.Rdata")
maf2 <- read.maf(maf = "data/0458c57f-316c-4a7c-9294-ccd11c97c2f9/TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.maf")

for (i in levels(sample_group$group)) {
print(i)
a <- subsetMaf(TCGAdataMAF,genes = i,)
tsbi1 <- a@clinical.data$Tumor_Sample_Barcode
iid <- sample_group$sample[sample_group$group == i] %>% substring(1,12)
tsbi2 <- intersect(x = tsbi1,iid)
a1 <- subsetMaf(maf = TCGAdataMAF,tsb = tsbi1)
oncoplot(maf = a1)
graph2ppt(file = "output/plots/oncoplot_top10_mutatedgenes.pptx",app = T)
}

a <- subsetMaf(TCGAdataMAF,genes = c("KEAP1",'NFE2L2','CUL3'))
tsbi1 <- a@clinical.data$Tumor_Sample_Barcode
iid <- sample_group$sample[sample_group$group != 'Wild'] %>% substring(1,12)
tsbi2 <- intersect(tsbi1,iid)
a1 <- subsetMaf(maf = TCGAdataMAF,tsb = iid)
genes_show <- union(c('KEAP1','NFE2L2','CUL3'), getGeneSummary(a1)$Hugo_Symbol[1:20])
oncoplot(maf = a1,genes = genes_show,keepGeneOrder = T)
graph2ppt(file = "output/plots/oncoplot_top10_mutatedgenes.pptx",app = T)
graph2pdf(file = "output/plots/oncoplot_top10_mutatedgenes_Mutgroup")

a <- subsetMaf(TCGAdataMAF,genes = c("KEAP1",'NFE2L2','CUL3'))
alltsb <- TCGAdataMAF@clinical.data$Tumor_Sample_Barcode
tsbi1 <- alltsb[ !alltsb %in% a@clinical.data$Tumor_Sample_Barcode]
iid <- sample_group$sample[sample_group$group == 'Wild'] %>% substring(1,12)
tsbi2 <- intersect(x = tsbi1,iid)
a1 <- subsetMaf(maf = TCGAdataMAF,tsb = iid)
genes_show <- union(c('KEAP1','NFE2L2','CUL3'),getGeneSummary(a1)$Hugo_Symbol[1:20])
oncoplot(maf = a1,genes = genes_show,keepGeneOrder = T)
graph2pdf(file = "output/plots/oncoplot_top10_mutatedgenes_Wildgroup")
graph2ppt(file = "output/plots/oncoplot_top10_mutatedgenes.pptx",app = T)

# 这里有个isTCGA选项，选择T，它会自动精简我们的样本barcode成前12位
# subsetmaf 给maf文件按照我们自己的样本取子集的maf
{
#   # subsetmaf 把原来的TCGA-LUADmaf文件根据TSB取子集。
#   TCGAdataMAF<-read.maf("data/LUAD.maf",isTCGA = T)
#   tsb_keap1 <- sample_group$sample[sample_group$group == "KEAP1"] |> substring(1,12)
#   tsb_nfe2l2 <- sample_group$sample[sample_group$group == 'NFE2L2'] |> substring(1,12)
#   tsb_cul3 <- sample_group$sample[sample_group$group == "CUL3"] |> substring(1,12)
#   tsb_wt <- sample_group$sample[sample_group$group == "Wild"] |> substring(1,12)
# 
#   keap1_maf <- subsetMaf(maf = TCGAdataMAF,tsb = tsb_keap1)
#   nfe2l2_maf <- subsetMaf(maf = TCGAdataMAF,tsb = tsb_nfe2l2)
#   cul3_maf <- subsetMaf(maf = TCGAdataMAF,tsb = tsb_cul3)
#   wild_maf <- subsetMaf(maf = TCGAdataMAF,tsb = tsb_wt)
#   save(TCGAdataMAF,keap1_maf,nfe2l2_maf,cul3_maf,wild_maf,
#        file = "outdata/4groups_mafs.Rdata")
}



# oncoplot ----------------------------------------------------------------
oncoplot(maf = keap1_maf, top = 15)
graph2ppt(file = "output/plots/oncoplot_top10_mutatedgenes.pptx",append = T)
oncoplot(maf = nfe2l2_maf, top = 15)
graph2ppt(file = "output/plots/oncoplot_top10_mutatedgenes.pptx",append = T)
oncoplot(maf = cul3_maf, top = 10)
graph2ppt(file = "output/plots/oncoplot_top10_mutatedgenes.pptx",append = T)
oncoplot(maf = wild_maf, top = 10)
graph2ppt(file = "output/plots/oncoplot_top10_mutatedgenes.pptx",append = T)


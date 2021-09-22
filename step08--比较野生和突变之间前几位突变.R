rm(list = ls())
p_load(maftools, GenVisR, ggplot2, reshape2, data.table, ggpubr, export)
cln_survival <- read.csv("outdata/cln_xena.csv")
cln <- read.csv("outdata/phe_finish.csv") %>% select(-1) %>% 
  dplyr::rename(Tumor_Sample_Barcode = sample) %>% 
  select(1,6) %>% 
  mutate(group = as.factor(group))
cln_annotation <- cln %>% 
  mutate(Mut = na_if(group, "wild"),
         Wild = na_if(group, "mut"),) 
top21 <- read.csv(file = "outdata/top21genes.csv", header = F) 
top21[1,1] <- "TP53" 
TCGAdataMAFA.plus.cn<-read.maf(maf = 'outdata/TCGAdataMAFA(mut).plus.cn_maftools.maf',
                               clinicalData = cln_annotation)

TCGAdataMAFB.plus.cn<-read.maf(maf = 'outdata/TCGAdataMAFB.plus.cn_maftools.maf',
                               clinicalData = cln_annotation)

#通过写出data再读回去的方式分亚组
TCGAdataMAF<-read.maf("data/LUAD.maf")
data <- TCGAdataMAF@data
data <- as.data.frame(data)
data$Tumor_Sample_Barcode <- gsub('-','.',data$Tumor_Sample_Barcode)
data$Tumor_Sample_Barcode <- gsub('............$','',data$Tumor_Sample_Barcode)
data$Tumor_Sample_Barcode <- gsub( '\\.','-', x = data$Tumor_Sample_Barcode)
data <- data[data$Tumor_Sample_Barcode %in% cln_survival$sample, ]
TCGAdataMAF_493 <- read.maf(maf = data, clinicalData = cln_annotation)

## 先搞出前21
genes <- top21[,1] %>% as.character()

oncoplot(maf = TCGAdataMAF_493,top = 21 ,clinicalFeatures = "group",sortByAnnotation = F)
oncoplot(maf = TCGAdataMAF_493, clinicalFeatures = "Mut",sortByAnnotation = F)
oncoplot(maf = TCGAdataMAF_493, clinicalFeatures = "Wild",sortByAnnotation = F)
oncoplot(maf = TCGAdataMAFA.plus.cn, genes = genes)
oncoplot(maf = TCGAdataMAFB.plus.cn, genes = genes)


# 最后选了10个基因，是mut和wild看起来有区别的。然后也是比较突变高得基因
genes <- c("TTN", "RYR2", "CSMD3","USH2A", "SPTA1", "ZFHX4", "NAV3","SPEF2", "SNTG2", "EGFR")
oncoplot(maf = TCGAdataMAF_493,genes = genes, clinicalFeatures = "Mut",sortByAnnotation = F)
oncoplot(maf = TCGAdataMAF_493,genes = genes, clinicalFeatures = "Wild",sortByAnnotation = F)
oncoplot(maf = TCGAdataMAFA.plus.cn, genes = genes, keepGeneOrder = T,clinicalFeatures = "group",annotationColor = groupcolors)
graph2ppt(file = "output/plots/突变组10基因oncoplot.pptx")
graph2pdf(file = "output/plots/突变组10基因oncoplot")

{
# laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
# laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 
# laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
# fabcolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')
# names(fabcolors) = c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7")
# fabcolors = list(FAB_classification = fabcolors)
# oncoplot(maf = laml, clinicalFeatures = 'FAB_classification', sortByAnnotation = TRUE, annotationColor = fabcolors)
}
groupcolors <- RColorBrewer::brewer.pal(3,name = 'Spectral') %>% .[-3]
names(groupcolors) <- c("wild","mut")
groupcolors <- list(group = groupcolors)
wildcolors = list(Wild = c(wild = "BLUE"))
# groupcolors = list(group = c(Wild = "#0000ff", Mut = "#FF0000"))
# oncoplot(maf = TCGAdataMAFB.plus.cn, genes = genes, keepGeneOrder = T,clinicalFeatures = "Wild",annotationColor = wildcolors)
oncoplot(maf = TCGAdataMAFB.plus.cn, genes = genes, keepGeneOrder = T,clinicalFeatures = "group",annotationColor = groupcolors)
graph2ppt(file = "output/plots/野生组10基因oncoplot.pptx")
graph2pdf(file = "output/plots/野生组10基因oncoplot")

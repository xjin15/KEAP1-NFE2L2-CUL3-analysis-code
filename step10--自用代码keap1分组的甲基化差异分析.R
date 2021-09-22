rm(list = ls())
#加载数据
library(tidyverse)
data_beta=read_tsv(file = "data/TCGA-LUAD.methylation450.tsv")
data_samp=read_tsv(file = "data/TCGA-LUAD.GDC_phenotype.tsv.gz")
cln_data <- read_csv("outdata/cln_xena.csv") %>% 
  select(sample,group) %>% 
  dplyr::rename(Sample_Name = sample,
                Sample_Group = group)
dim(data_beta)
## [1] 485577    504

ID=intersect(cln_data$Sample_Name, colnames(data_beta))

#提取具有共同交集的sample信息
cln_data <- cln_data[cln_data$Sample_Name %in% ID,] 
#查看整理好的表型数据
dim(cln_data)
data_samp <- cln_data
#提取具有共同交集的beta矩阵
names(data_beta)[1]="CpG"
data_beta <- column_to_rownames(data_beta,"CpG") 
data_beta=data_beta[ , ID]
####清洁数据不易，记得保存
save(data_beta, data_samp, file = "methylation.rds") # RDS文件是二进制文件，保存起来很小,比Rdata还小

#加载数据和包
rm(list = ls())
library(ChAMP)
library(tidyverse)
load("methylation.rds")     ######加载速度快得飞起！
data_samp <- as.data.frame(data_samp) %>% arrange(Sample_Group)

###beta矩阵的处理
##1.保持beta矩阵的列名和表型数据的样本名一致，这很重要！
##2.beta矩阵里不能有NA
##3.beta矩阵中最小值如果为0，可以加上0.00001防止报错
##4.beta矩阵必须为矩阵数据格式
data_order=data_beta[ , data_samp$Sample_Name]
data_order=as.matrix(data_order)
sum(is.na(data_order))
## [1] 39199574
data_order <- na.omit(data_order)
data_order=data_order+0.000001
###CHAMP包导入
myLoad=champ.filter(beta = data_order,pd = data_samp)

## [===========================]
## [You may want to process champ.QC() next.]
##成功导入后及时保存
save(myLoad,file="meth_load.rds")
# load("meth_load.rds") 
champ.QC(beta=myLoad$beta, pheno = myLoad$pd$Sample_Group)
myNorm <- champ.norm(beta=myLoad$beta, plotBMIQ = T,  arraytype="450K",cores= 8)
champ.SVD()#####有批次效应的才需要做，没有就不用做了

###另外800多个病人，速度会比较慢。。。。
# ##标准化，我们使用文章中使用的BMIQ法进行标准化。core值根据自己分析环境设置
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores= 8)
# ##同样及时保存，如果是自己的电脑就不要跑这一步了，会崩的,所以我就不演示了
save(myNorm,file="meth_norm.rds")
 # #标准化后得到新的数据
##标准化后数据质控
champ.QC(beta=myLoad$beta, pheno = myLoad$pd$Sample_Group)
####甲基化差异位点分析
myDMP <- champ.DMP(beta = myLoad$beta ,pheno=myLoad$pd$Sample_Group)
#查看我们的分析结果
head(myDMP[[1]])
DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=myLoad$pd$Sample_Group)
####差异甲基化区域，DMR
myDMR <- champ.DMR(beta=myLoad$beta,pheno=myLoad$pd$Sample_Group,method="Bumphunter")
DMR.GUI(DMR=myDMR,beta=myNorm,pheno=myLoad$pd$Sample_Group)

#查看差异甲基化区域
head(myDMR$BumphunterDMR)
myGSEA <- champ.GSEA()
myGSEA.ebGSEA <- champ.ebGSEA(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="450K")
save(myGSEA, myGSEA.ebGSEA, file = "champ_gsea.rds")
myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group)
## Not run: 
# myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))
# myNorm <- champ.norm()
# myDMP <- champ.DMP()
# myDMR <- champ.DMR()
# myGSEA <- champ.GSEA()
# 文章中的条件筛选包括：
# 
# （1）adj.P<0.05,delta |β|>0.2
# 
# （2）位于启动子区域（5’-UTR, TSS200, TSS1500 and 1stExon），那我们开始把。

##########差异甲基化位点提取#############
feature_pro=c("1stExon","5'UTR","TSS1500","TSS200")
data_dif <- myDMP[[1]] %>% filter(adj.P.Val<0.05&abs(deltaBeta)>0.2)
# data_dif <- data_dif[data_dif$feature%in%feature_pro,] # 位于启动子区域
data_dif$group <- ifelse(data_dif$deltaBeta > 0, "up", "down")
write.csv(data_dif, "outdata/DMP_champ.csv")
#加载包
library(pheatmap)
#设置注释列名的文件ann_col
ann_col <- data_samp
ann_col1 <- ann_col %>% 
  column_to_rownames(var = "Sample_Name") %>% 
  arrange(Sample_Group)
##根据差异甲基化位点提取beta矩阵子集
data_heat=data_beta[rownames(data_dif),rownames(ann_col1)]
pheatmap(data_heat,#color =colorRampPalette(c("navy", "white","firebrick3"))(10),
         cluster_rows = T,cluster_cols = F,legend = T,show_rownames = T,
         show_colnames = F,annotation_col = ann_col1)
graph2ppt(file = "output/差异甲基化位点热图.pptx")
graph2pdf(file = "output/差异甲基化位点热图.pdf")

# 甲基化的注释
# test <- data_heat
# rownames(test) <- hm450.manifest.hg19$gene[match(rownames(test), hm450.manifest.hg19$probeID)]
# 
# 
# hm450.manifest.hg19$probeID == ""
# id <- rownames(test)
# hm450.manifest.hg19$gene_HGNC[match(id, hm450.manifest.hg19$probeID)]
# probe.features$gene[rownames(probe.features)==id]
# 
# probe.features[id,"gene"]
# gtfmethy <- fread(file = "data/illuminaMethyl450_hg38_GDC")
# gtfmethy$gene[match(id,gtfmethy$`#id`)]
# 
library(DOSE)
data(geneList)
x <- gseDO(geneList)
gseaplot(x, geneSetID=1)
gseaplot()

#####################生信星球的代码 1.差异分析######################
library(ChAMP)
library(tibble)
# 差异分析
group_list <- data_samp$Sample_Group
myDMP2 <- champ.DMP(beta = myLoad$beta ,pheno=group_list)
df_DMP <- myDMP2$mut_to_wild
df_DMP=df_DMP[df_DMP$gene!="",] # 去掉基因名为空的69539 - 48446
logFC_t <- 0.45     #推荐用deltabeta值替代logFC，就是甲基化信号值的差值。
P.Value_t <- 0.05
df_DMP$change <- ifelse(df_DMP$adj.P.Val < P.Value_t & abs(df_DMP$deltaBeta) > 0.2,
                        ifelse(df_DMP$deltaBeta > 0.2 ,'UP','DOWN'),'NOT') 
table(df_DMP$change) 

##############差异甲基化位点火山图###############
library(dplyr)
library(ggplot2)
dat  = df_DMP
for_label <- dat%>% head(3)
p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(adj.P.Val))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p 
#####火山图太烂了。差异几乎没有 logFC基本在0.2附近
volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = rownames(for_label) ),
    data = for_label,
    color="black"
  )
volcano_plot



#############差异甲基化的热图##############
cg <-  rownames(df_DMP[df_DMP$change != "NOT",])
plot_matrix <- myNorm[cg,]
annotation_col <- data.frame(Sample=pd$group_list) 
rownames(annotation_col) <- colnames(plot_matrix)
ann_colors = list(Sample = c(Normal="#4DAF4A", Tumor="#E41A1C"))

library(pheatmap)
pheatmap(plot_matrix,show_colnames = T,
         annotation_col = annotation_col,
         border_color=NA,
         color = colorRampPalette(colors = c("white","navy"))(50),
         annotation_colors = ann_colors,show_rownames = F)
##############将差异甲基化的位点拿来做生存分析#############

rm(list = ls())
library(data.table)
library(stringr)
library(survival)
library(survminer)
load("./Rdata/step2_filtered_pd_myNorm.Rdata")
load("./Rdata/step3.df_DMP.Rdata")

cg <- rownames(df_DMP[df_DMP$change != "NOT",])
myNorm_tumor <- myNorm[cg,]
##########logranktest批量生存分析##########
suv_dat <- data.table::fread("./raw_data/TCGA-HNSC.survival.tsv.gz")
suv_dat$sample = str_sub(suv_dat$sample,1,15)
suv_dat <- suv_dat[suv_dat$sample %in% colnames(myNorm_tumor),]
suv_dat <- suv_dat[str_sub(suv_dat$sample,14,15)=="01",] 
suv_dat = merge(pd,suv_dat,by.x = "sampleID",by.y = "sample")
myNorm_tumor <- myNorm_tumor[,suv_dat$sample]
identical(colnames(myNorm_tumor),suv_dat$sample)
#> [1] TRUE
library(survival)
logrankP <- apply(myNorm_tumor, 1, function(x){
  #x <- myNorm_tumor[1,]
  suv_dat$group <- ifelse(x>mean(x),"High","Low")
  res <- coxph(Surv(OS.time, OS)~group, data=suv_dat)
  beta <- coef(res)
  se <- sqrt(diag(vcov(res)))
  p.val <- 1 - pchisq((beta/se)^2, 1)
})
table(logrankP<0.05) #110个CpG位点
#> 
#> FALSE  TRUE 
#>  1676   110

##########得到110个对生存影响显著的差异甲基化位点，取前20个画图。####
surv_gene <- names(sort(logrankP))[1:20] 
choose_matrix <- myNorm[surv_gene,]
annotation_col <- data.frame(Sample=pd$group_list) 
rownames(annotation_col) <- colnames(choose_matrix)
ann_colors = list(Sample = c(Normal="#4DAF4A", Tumor="#E41A1C"))

library(pheatmap)
pheatmap(choose_matrix,show_colnames = T,
         annotation_col = annotation_col,
         border_color=NA,
         color = colorRampPalette(colors = c("white","navy"))(50),
         annotation_colors = ann_colors)

############也可以画画生存分析的图##########

gs=head(surv_gene,4)
exprSet = myNorm_tumor
meta = suv_dat
splots <- lapply(gs, function(g){
  meta$gene=ifelse(exprSet[g,]>median(exprSet[g,]),'high','low')
  sfit1=survfit(Surv(OS.time, OS)~gene, data=meta)
  ggsurvplot(sfit1,pval =TRUE, data = meta)
}) 
arrange_ggsurvplots(splots, print = TRUE,  
                    ncol = 2, nrow = 2)


############富集分析######
# 利用ChAMP包对过滤后的数据做了差异甲基化位点分析。如果是肿瘤数据的话，可以加一步生存分析。
rm(list = ls())
load(file = 'Rdata/step3.df_DMP.Rdata')
library(ggplot2)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
##ID转换
{length(unique(df_DMP$gene))
#> [1] 16338
s2e <- bitr(unique(df_DMP$gene), fromType = "SYMBOL",
            toType = c( "ENTREZID"),
            OrgDb = org.Hs.eg.db)

df_DMP=merge(df_DMP,s2e,by.y='SYMBOL',by.x='gene')
table(!duplicated(df_DMP$ENTREZID))
#> 
#> FALSE  TRUE 
#> 83778 14188

gene_up= unique(df_DMP[df_DMP$change == 'UP','ENTREZID'] )
gene_down=unique(df_DMP[df_DMP$change == 'DOWN','ENTREZID'] )
gene_diff=c(gene_up,gene_down)
gene_all=unique(df_DMP$ENTREZID)
}
#####富集分析

kkgo_file = "./Rdata/kkgo_file.Rdata"
if(!file.exists(kkgo_file)){
  kk <- enrichKEGG(gene         = gene_diff,
                   universe     = gene_all,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)
  go <- enrichGO(gene_diff, OrgDb = "org.Hs.eg.db", ont="all") 
  save(kk,go,file = kkgo_file)
}
load(kkgo_file)

barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=45))
dotplot(kk)

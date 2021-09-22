rm(list = ls())
p_load(tidyverse,ggplot2, limma, export)

library(pheatmap)
# 准备数据
allid <- read.csv("outdata/allid_493.txt") %>% .[,2]
rt <- read_csv(file = "outdata/mrna_fpkm_finish.txt",col_names = T) %>% 
  column_to_rownames("X1")
gtf_gene <- read.csv(file = "output/gtfmRNA22.txt", header = T,sep = "\t")
## 选择proteincoding的蛋白质
gtf_gene_protein <- gtf_gene %>% 
  dplyr::rename( Ensembl_ID = gene_id,
                 Symbol  = gene_name,
                 Biotype = gene_type  ) 
rownames(rt) <- substring(text = rownames(rt), first = 1, last = 15)
stopifnot(
  {
    duplicated(rownames(rt))== F
  }
)
# 查看行和列数
dim(rt) #57502 493(mut111,wild382)
matchgene <- intersect(rownames(rt), gtf_gene_protein$Ensembl_ID)
mrna <- rt[matchgene, ] %>% as.matrix()### 57502
# 可以转换ID
rownames(mrna) <- gtf_gene_protein$Symbol[match(x = rownames(mrna),table = gtf_gene_protein$Ensembl_ID)]
# 重复的基因名字进行合并且取均值
mrna <- avereps(mrna) 
dim(mrna)#55796
mrna <- t(mrna)
mrna[1:3,1:3]
#(2)固有免疫功能、抗原提呈能力的比较
# immugene <- c('IRF3','MYD88','TICAM1','TLR3','TLR5','TLR7',
#               'TLR8','DDX58','IFIH1','MAVS','CLEC7A','CLEC4E',
#               'CD209','CLEC10A','NLRP3','AIM2','PYCARD', #innate免疫共17个分子initiation of innate immunity
#               'HLA-A','HLA-B','HLA-C','TAP1','TAP2','B2M','HLA-DPA1',
#               'HLA-DPB1','HLA-DQA1','HLA-DQB1','HLA-DQB2',
#               'HLA-DRB4','HLA-DRB6','HLA-E','HLA-F','HLA-J', #抗原提呈系列MHC-I/II antigen-presenting process的HLA共16个分子
#               'CTLA4','CD48','PDCD1','LAG3','HAVCR1','BTN2A2',
#               'LAIR1','BTN3A1','PDCD1LG2','BTN1A1','VTCN1','BTNL2', #免疫检查点12个immune co-inhibitors
#               'ICOS','TNFRSF9','CD70','CD80','TNFRSF13B','CD27',
#               'SLAMF1','TNFRSF14','TNFRSF8','ICOSLG','CD58','TNFSF13','TNFSF15',
#               'TNFRSF4','HAVCR1','TNFSF18','TNFSF4','CD86' ) #co-stimulators
# 数据读取
{
immugeneMHC <- data.table::fread(
"Antigen presentation
HLA-A
HLA-B
HLA-C
HLA-DPA1
HLA-DPB1
HLA-DQA1
HLA-DQA2
HLA-DQB1
HLA-DQB2
HLA-DRA
HLA-DRB1
HLA-DRB3
HLA-DRB4
HLA-DRB5
MICA
MICB
") %>% 
    pull()
  
  } # 读取抗原呈递的MHC
{
  immugeneINHI <- data.table::fread(
  "inhibitory
  ADORA2A
ARG1
BTLA
CD274
CD276
CTLA4
EDNRB
HAVCR2
IDO1
IL10
IL13
IL4
KIR2DL1
KIR2DL2
KIR2DL3
LAG3
PDCD1
SLAMF7
TGFB1
TIGIT
VEGFA
VEGFB
C10orf54
VTCN1
"
) %>% pull()
  } # 读取免疫检查点 抑制免疫作用的
{immugeneSTIM <- data.table::fread(
  "Stimulatory
  BTN3A1
BTN3A2
CCL5
CD27
CD28
CD40
CD40LG
CD70
CD80
CX3CL1
CXCL10
CXCL9
ENTPD1
GZMA
HMGB1
ICAM1
ICOS
ICOSLG
IFNA1
IFNA2
IFNG
IL12A
IL1A
IL1B
IL2
IL2RA
ITGB2
PRF1
SELP
TLR4
TNF
TNFRSF14
TNFRSF18
TNFRSF4
TNFRSF9
TNFSF4
TNFSF9
"
) %>% pull()} # 读取免疫检查点 促进免疫作用的
immugene2 <- rbind(immugeneINHI,immugeneSTIM,use.names = F) %>% pull() %>% dput()
immugene <- c('IRF3','MYD88','TICAM1','TLR3','TLR5','TLR7',
              'TLR8','DDX58','IFIH1','MAVS','CLEC7A','CLEC4E',
              'CD209','CLEC10A','NLRP3','AIM2','PYCARD', #innate免疫共17个分子initiation of innate immunity
              'HLA-A','HLA-B','HLA-C','TAP1','TAP2','B2M','HLA-DPA1',
              'HLA-DPB1','HLA-DQA1','HLA-DQB1','HLA-DQB2',
              'HLA-DRB4','HLA-DRB6','HLA-E','HLA-F','HLA-J', #抗原提呈系列MHC-I/II antigen-presenting process的HLA共16个分子
              "ADORA2A", "ARG1", "BTLA", "CD274", "CD276", "CTLA4", "EDNRB", 
"HAVCR2", "IDO1", "IL10", "IL13", "IL4", "KIR2DL1", "KIR2DL2", 
"KIR2DL3", "LAG3", "PDCD1", "SLAMF7", "TGFB1", "TIGIT", "VEGFA", 
"VEGFB", "C10orf54", "VTCN1", # 免疫检查点抑制基因 24个
"BTN3A1", "BTN3A2", "CCL5", "CD27", 
"CD28", "CD40", "CD40LG", "CD70", "CD80", "CX3CL1", "CXCL10", 
"CXCL9", "ENTPD1", "GZMA", "HMGB1", "ICAM1", "ICOS", "ICOSLG", 
"IFNA1", "IFNA2", "IFNG", "IL12A", "IL1A", "IL1B", "IL2", "IL2RA", 
"ITGB2", "PRF1", "SELP", "TLR4", "TNF", "TNFRSF14", "TNFRSF18", 
"TNFRSF4", "TNFRSF9", "TNFSF4", "TNFSF9" ) #免疫检查点激活基因37个

Innateimmune <- mrna[, match(immugene,colnames(mrna),nomatch = 0)] %>% as.data.frame()
is.na(Innateimmune) %>% table()

Innateimmune$group <- "Wild"
Innateimmune$group[rownames(Innateimmune) %in% allid[1:111]] <- "Mut"
table(Innateimmune$group)
Innateimmune[, 59] <- NULL
# 过滤掉第一组，保留第二、第三组内的数据
# Innateimmune<-Innateimmune[match(GSEpred$X,rownames(Innateimmune)),]
# Innateimmune$TMEcluster<-GSEpred$GSEpred
################ 画直方图，比较两组的差异###############
bodat <- Innateimmune %>% 
  select(group, everything())
bodat2 <- bodat %>% 
  pivot_longer(cols = c(2:ncol(bodat)))
x <- compare_means(value ~ group, group.by = "name",data = bodat2)
xnosig <- x[which(x$p > 0.05),]

ggboxplot(bodat2,x="group",y="value",fill = "group",size = 0.1,
          facet.by = "name",
         palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("InfiltrationScore")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  # scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6))
graph2ppt(file = "output/plots/immuneAI",append = T,asp = 1,width = 12)
######### 只看免疫检查点的基因########
Innateimmune <- mrna[, match(x = c(immugeneINHI,immugeneSTIM),colnames(mrna),nomatch = 0)] %>% as.data.frame()
is.na(Innateimmune) %>% table()
allid <- read.csv("outdata/allid_493.txt") %>% .[,2]
Innateimmune$group <- "Wild"
Innateimmune$group[rownames(Innateimmune) %in% allid[1:111]] <- "Mut"
table(Innateimmune$group)

bodat <- Innateimmune %>% 
  select(group, everything())
bodat2 <- bodat %>% 
  pivot_longer(cols = c(2:ncol(Innateimmune)))
x <- compare_means(value ~ group, group.by = "name",data = bodat2)
xnosig <- x[which(x$p > 0.001),]
bodatsig <- bodat2 %>% 
  filter(!name %in% xnosig$name)
ggboxplot(bodatsig,x="group",y="value",fill = "group",size = 0.1,
          facet.by = "name",
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("InfiltrationScore")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  # scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6))


bodatpd <- bodat2 %>% 
  filter(name %in% c("PDCD1","CTLA4","CD274","LAG3","TGFB1","VEGFA","VEGFB"))
ggboxplot(bodatpd,x="name",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("log2(FPKM+1)")+
  theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 0, vjust = 0))+
  # scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 0, vjust = 0.6))
graph2ppt(file = "output/plots/6免疫基因")
############## 画热图 求的是分组均数到总中位数的距离###################
InnateimmuneA <- Innateimmune %>% filter(group == "Mut")
InnateimmuneA$group <- NULL
InnateimmuneA[] <- lapply(InnateimmuneA, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
MeanA <- as.data.frame(apply(InnateimmuneA,2,mean))

InnateimmuneB <- Innateimmune%>%filter(group == "Wild")
InnateimmuneB$group <- NULL
InnateimmuneB[] <- lapply(InnateimmuneB, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
MeanB <- as.data.frame(apply(InnateimmuneB,2,mean))

Mean <- cbind( MeanA, MeanB)
colnames(Mean)<-c("Mut","Wild")

Innateimmune$group <- NULL
Innateimmune[] <- lapply(Innateimmune, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
Median <- as.data.frame(apply(Innateimmune,2,median))
Mean$Mut <- Mean$Mut - Median$`apply(Innateimmune, 2, median)`
Mean$Wild <- Mean$Wild - Median$`apply(Innateimmune, 2, median)`

# Mean2 <- cbind( MeanA, MeanB)
# colnames(Mean2)<-c("Mut","Wild")
# pheatmap(Mean2,scale = "column",cluster_rows = F, cluster_cols = F,
#          show_colnames =T,show_rownames = T, 
#          color = colorRampPalette(color.key)(50),
#          cutree_cols = 0,
#          cellwidth = 8,
#          cellheight = 8,)
#          
         #z化
#Mean_Z=t(scale(t(Mean)))
#Mean_Z<-as
#Mean_Z<-as.data.frame(Mean_Z)
#Mean_Z[Mean_Z>  2]=2 #限定上限，使表达量大于0.05的等于0.05
#Mean_Z[Mean_Z< -2]= -2 #限定下限，使表达量小于-2的等于-2
#Mean_log<-10*Mean
# dd1 <- Mean
# dd2=apply(dd1,2,rescale )        ##归一化
# dd2=t(dd2)    
# pheatmap(dd2)
#热图绘制
color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pheatmap(Mean,cluster_rows = F, cluster_cols = F,
         show_colnames =T,show_rownames = T, 
         color = colorRampPalette(color.key)(50),
         cutree_cols = 0,
         #annotation_colors = color.annotation,
         #annotation_col = annotation_col,
         cellwidth = 50,
         cellheight = 18,
         gaps_row = c(17,32,67),
         filename = "output/plots/mianyijinrun.pdf"
)
dev.off()
graph2pdf(file = "output/plots/mianyijinrun.pdf",width = 20)

tmean <- t(Mean)
pheatmap(t(Mean[!rownames(Mean) %in% xnosig$name,]),
         cluster_rows = F, cluster_cols = F,
         show_colnames =T,show_rownames = T, 
         color = colorRampPalette(color.key)(50),
         cutree_cols = 0,
         #annotation_colors = color.annotation,
         #annotation_col = annotation_col,
         cellwidth = 18,
         cellheight = 50,
         gaps_col =  c(13,25,41), 
         filename = "output/plots/mianyijinrun2.pdf")
# 筛选出有意义的
xnosig$name
Mean[!rownames(Mean) %in% xnosig$name,]
pheatmap(Mean[!rownames(Mean) %in% xnosig$name,],
         cluster_rows = F, cluster_cols = F,
         show_colnames =T,show_rownames = T, 
         color = colorRampPalette(color.key)(50),
         cutree_cols = 0,
         #annotation_colors = color.annotation,
         #annotation_col = annotation_col,
         cellwidth = 50,
         cellheight = 18,
         gaps_row = c(13,25,41),
         filename = "output/plots/mianyijinrun.pdf"
)
# pheatmap(Mean,cluster_rows = F, cluster_cols = F,treeheight_col = 0,
#          show_colnames =F,show_rownames = T, 
#          color = colorRampPalette(color.key)(50),
#          cutree_cols = 0,
#          #annotation_colors = color.annotation,
#          #annotation_col = annotation_col,
#          cellwidth = 50,
#          cellheight = 18,
#          gaps_row = c(17,33,45),
#          file = "/Users/biguoshu/desktop/pheatmap.pdf",plot = "pdf"
# )




# 免疫相关基因-by SQH -----------------------------------------------------------
imgsqh <- read.csv(file = "cibersort/immugene.csv")
Innateimmune <- mrna[, match(x = unique(imgsqh$Symbol),colnames(mrna),nomatch = 0)] %>% as.data.frame()
is.na(Innateimmune) %>% table()
Innateimmune <- Innateimmune[,colMeans(Innateimmune) > 0.5]
allid <- read.csv("outdata/allid_493.txt") %>% .[,2]
Innateimmune$group <- "Wild"
Innateimmune$group[rownames(Innateimmune) %in% allid[1:111]] <- "Mut"
table(Innateimmune$group)

bodat <- Innateimmune %>% 
  select(group, everything())
bodat2 <- bodat %>% 
  pivot_longer(cols = c(2:ncol(Innateimmune)))
x <- compare_means(value ~ group, group.by = "name",data = bodat2)
xnosig <- x[which(x$p > 0.05),]
xhsig <- x[which(x$p < 0.001),]
xhhsig <- x[which(x$p.adj < 0.001),]
xhsig$func <- imgsqh$Category[match(x = xhsig$name, imgsqh$Symbol,nomatch = 0) ]
xhhsig$func <- imgsqh$Category[match(x = xhhsig$name, imgsqh$Symbol,nomatch = 0) ]
table(xhsig$func)
table(xhhsig$func)
table(imgsqh$Category)
xhhsig <- xhhsig[order(xhhsig$p.adj, decreasing = F),]

quhua <- xhhsig %>% 
  filter(func == "Chemokine_Receptors"|
           func == "Chemokines") %>% 
  select(name)
bodatsig <- bodat2 %>% 
  filter(name %in% quhua$name)

all_diff <- read.csv(file = "outdata/mrna_all_diff.csv",row.names = 1)
all_diff$group = "not-significant"
#将adj.P.Val<0.05,logFC>0.5的基因设置为显著上调基因
#将adj.P.Val<0.05,logFC<0.5的基因设置为显著下调基因
all_diff$group[which((all_diff$adj.P.Val<0.05) & (all_diff$logFC> 0.5))] = "up-regulated"
all_diff$group[which((all_diff$adj.P.Val<0.05) & (all_diff$logFC< -0.5))] = "down-regulated"
#查看上调和下调基因的数目
table(all_diff$group)
all_diff$group<-as.factor(all_diff$group)
all_diff$logP<- -log10(all_diff$adj.P.Val)
top10 <- intersect(all_diff$ID[all_diff$group != "not-significant"], xhhsig$name[1:11])
### 选差异top10
bodattop <- bodat2[bodat2$name %in% top10,]

 ggboxplot(bodattop,x="name",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("log2(FPKM+1)")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6))
graph2ppt(file = "output/plots/6免疫基因.pptx", append= T)
xhhsig$func[which(xhhsig$name %in% top10)]
bodatpd <- bodat2 %>% 
  filter(name %in% c("PDCD1","CTLA4","CD274","TGFB1","VEGFA","VEGFB"))
ggboxplot(bodatpd,x="name",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("log2(FPKM+1)")+
  theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 0, vjust = 0))+
  # scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 0, vjust = 0.6))
############### 相关性图############
aidat <- read.table(file = "cibersort/immucellAI/ImmuCellAI_icb_result.txt",header = T)
aidat <- aidat %>% 
  select(Th17,nTreg, InfiltrationScore)
TCGAtmb<-read.csv("outdata/TCGA_LUAD_tmv.csv", row.names = 1)

tmb <- TCGAtmb %>% 
  filter(Tumor_Sample_Barcode %in% allid) %>% 
  select(Tumor_Sample_Barcode, total) %>% 
  rename(SampleID = Tumor_Sample_Barcode,
         tmb =  total)
rownames(tmb) <- tmb[,1]  
tmb <- tmb[allid, -1] 
repla <- mean(tmb,na.rm = T)
tmb <- replace_na(data = tmb,replace = repla)

corgene <- c("PDCD1","CTLA4","CD274","TGFB1","VEGFA","VEGFB", 
             "HLA-DQB2","CD40","CX3CL1","ICAM1","ITGB2")
cordat <- mrna[, match(corgene,colnames(mrna),nomatch = 0)]
cordat <- cbind(aidat,tmb,cordat)

cormat <- round(cor(cordat), 2) 

ggcorrplot(cormat, type = "lower",
           outline.col = "white",
           ggtheme = ggplot2::theme_bw(),
           )
graph2ppt(file = "output/plots/6免疫基因.pptx",append = T)
### 样式复杂的相关性图
p_load(PerformanceAnalytics)
my_data <- cordat
chart.Correlation(my_data, histogram=TRUE, pch=19)



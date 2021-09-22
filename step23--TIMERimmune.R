rm(list = ls())
library(tidyverse)
# fpkm转tpm --------------------------------------------------------------------
rt <- data.table::fread(file = "outdata/mrna_fpkm_finish.txt",
                        check.names = F,
                        header = T) %>% 
  column_to_rownames("V1")
fpkm <- rt
# log（fpkm+1）的表达数据去log化。
fpkm1 <- apply(fpkm, 2, FUN = function(x){2^x - 1}) # 去log化
# 把data.frame转换成matrix,因为这样就允许行名重复了
rt <- as.matrix(fpkm1)        

### 转换ID
gtf_gene <- read.csv(file = "output/gtfmRNA22.txt", header = T,sep = "\t")
## 选择proteincoding的蛋白质
gtf_gene_all <- gtf_gene %>% 
  dplyr::rename( Ensembl_ID = gene_id,
                 Symbol  = gene_name,
                 Biotype = gene_type  ) 
rownames(rt) <- substring(text = rownames(rt), first = 1, last = 15)
# 确认没重复的基因ENSEMBL ID
stopifnot(duplicated(rownames(rt)) == F)
rownames(rt) <- gtf_gene_all$Symbol[match(rownames(rt),gtf_gene_all$Ensembl_ID,nomatch = 0)]
duplicated(rownames(rt)) %>% table()
# 换成symbol后，有很多重复的symbol名字 
# 注意此时RT是matrix，很多操作无法做
mean_ofgene <- rowMeans(rt)
index <- order(mean_ofgene,decreasing = T)
rt <- rt[index,]
rt1 <- rt[unique(rownames(rt)),];anyDuplicated(rownames(rt1))
# rt2 <- rt[!duplicated(rownames(rt)),]
rt1[1:4,1:4]

tpm_luad = t(t(rt1)/colSums(rt1))*10^6
tpm_luad[1:4,1:4]
colSums(tpm_luad) # tpm的colSum一定都是10的6次方
load("outdata/sample_group.Rdata")
tpm_luad <- tpm_luad[,sample_group$sample]
colnames(tpm_luad)

write.csv(tpm_luad,file = "outdata/TPMofKEAP1LUAD.txt")
save(tpm_luad,file = "outdata/TPM_of_TCGA_luad.Rdata")


# tide网站下载的TCGA病人的免疫浸润数据 --------------------------------------------------
rm(list = ls())
tide_immune_data <- data.table::fread("../data/IMMUNE_data/infiltration_estimation_for_tcga.csv.gz")
dim(tide_immune_data)
colnames(tide_immune_data)
load("outdata/sample_group.Rdata")
sample_group$pid <- sample_group$sample %>% substring(1,15)
keap_immu_df2 <- tide_immune_data[tide_immune_data$cell_type %in% sample_group$pid,]
keap_immu_df <- tide_immune_data[match(sample_group$pid,tide_immune_data$cell_type,nomatch = 0),]
dim(keap_immu_df)
# [1] 492 120 
setdiff(sample_group$pid,keap_immu_df$cell_type) # 突变组少一人 TCGA-44-3917-01
test <- colnames(keap_immu_df)[-1]
test <- test %>% str_split(pattern = "_",n = 2,simplify = T)
test[,2] %>% unique() %>% dput
# TIDE数据库提供了几乎所有TCGA病人的几个主流预测免疫浸润指数的数据
# 包括"TIMER", "CIBERSORT", "CIBERSORT-ABS", "QUANTISEQ", "MCPCOUNTER", "XCELL", "EPIC"
table(test[,2])

# CIBERSORT CIBERSORT-ABS   EPIC   MCPCOUNTER QUANTISEQ TIMER XCELL 
# 22            22             8      11      11        6     39 



# EPIC-------------------------------------------------------------------------

 data <- keap_immu_df %>% column_to_rownames("cell_type")
data <- data %>% select(ends_with("EPIC"))
data1 <- data
sample_group$groupw <- factor(sample_group$group,labels = c("Mut",'Mut','Mut','Wild'))
k <- apply(data1,2,function(x) {sum(x == 0) < nrow(data1)/2})
table(k)

data2 <- as.data.frame(t(data1[,k]))

an_col <- data.frame(group = sample_group$groupw[match(colnames(data2),sample_group$pid)], 
                     row.names = colnames(data2))
# 热图不好看
pheatmap(data2[-8,],
         scale = "column",
         show_colnames = F,
         annotation_col = an_col,
         cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
         )
# 宽转长数据
dat2 <- data1 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

dat2$group <- an_col$group[match(dat2$Sample, rownames(an_col))]


# 画图数据

mypalette <- colorRampPalette(RColorBrewer::brewer.pal(8,"Set1"))


ggplot(dat2,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell_type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(8))


ggplot(dat2,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))
graph2ppt(file = "output/plots/EPIC.pptx",append = T,)

# 加上顺序
a2 = dat2 %>% 
  group_by(Cell_type) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cell_type)

dat2$Cell_type = factor(dat2$Cell_type,levels = a2)

ggplot(dat2,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))


library(ggpubr)
ggplot(dat2,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,vjust = 0.5))+
  scale_fill_manual(values = mypalette(8)[c(6,1)])+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),hide.ns = T,method = "kruskal.test")

graph2ppt("output/plots/EPIC.pptx",append = T,aspectr = 1)

compare_means(Proportion ~ group,group.by = "Cell_type", data = dat2)


# 只选有意义的作图
dat2sig <- dat2 %>% 
  filter(Cell_type %in% c("B cell_EPIC", 
                          "Cancer associated fibroblast_EPIC","T cell CD4+_EPIC",
                          "T cell CD8+_EPIC",
                          "Macrophage_EPIC"
                         ))

ggplot(dat2sig,aes(Cell_type,Proportion,fill = group)) +  
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  scale_y_continuous(limits = c(0,0.15))+
  theme(axis.text.x = element_text(angle=0,vjust = 0.5))+
  scale_fill_manual(values = mypalette(8)[c(6,1)])+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),hide.ns = T,method = "kruskal.test")

graph2ppt(file = "output/plots/EPIC.pptx",append = T,)




# QUANTISEQ --------------------------------------------------------------------

data <- keap_immu_df %>% column_to_rownames("cell_type")
data <- data %>% select(ends_with("QUANTISEQ")) 
data1 <- data
sample_group$groupw <- factor(sample_group$group,labels = c("Mut",'Mut','Mut','Wild'))
k <- apply(data1,2,function(x) {sum(x == 0) < nrow(data1)/2})
table(k)

data2 <- as.data.frame(t(data1[,k]))

an_col <- data.frame(group = sample_group$groupw[match(colnames(data2),sample_group$pid)], 
                     row.names = colnames(data2))
# 热图不好看
pheatmap(data2[-8,],
         scale = "row",
         show_colnames = F,
         annotation_col = an_col,
         cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
)
# 宽转长数据
dat2 <- data1 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

dat2$group <- an_col$group[match(dat2$Sample, rownames(an_col))]


# 画图数据

mypalette <- colorRampPalette(RColorBrewer::brewer.pal(8,"Set1"))


ggplot(dat2,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell_type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(11))


ggplot(dat2,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))
graph2ppt(file = "output/plots/EPIC.pptx",append = T,)

# 加上顺序
a2 = dat2 %>% 
  group_by(Cell_type) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cell_type)

dat2$Cell_type = factor(dat2$Cell_type,levels = a2)

ggplot(dat2,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))


library(ggpubr)
ggplot(dat2,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,vjust = 0.5))+
  scale_fill_manual(values = mypalette(8)[c(6,1)])+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),hide.ns = T,method = "kruskal.test")

graph2ppt("output/plots/EPIC.pptx",append = T,aspectr = 1)

x <- compare_means(Proportion ~ group,group.by = "Cell_type", data = dat2)
sigcell <- x$Cell_type[x$p.signif != "ns"] %>% as.character()
sigcell <- sigcell[!sigcell %in% "uncharacterized cell_QUANTISEQ"]
sigcell
# 只选有意义的作图
dat2sig <- dat2 %>% 
  filter(Cell_type %in% sigcell)
dat2sig$Cell_type <- dat2sig$Cell_type %>% str_remove(pattern = "_QUANTISEQ")

ggplot(dat2sig,aes(Cell_type,Proportion,fill = group)) +  
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  scale_y_continuous(limits = c(0,0.15))+
  theme(axis.text.x = element_text(angle=0,vjust = 0.5))+
  scale_fill_manual(values = mypalette(8)[c(5,2)])+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),hide.ns = T,method = "kruskal.test")

graph2ppt(file = "output/plots/EPIC.pptx",append = T,)



# 保存最重要的免疫数据 --------------------------------------------------------------


save(tide_immune_data,file = "../data/IMMUNE_data/TIDE.Rdata")


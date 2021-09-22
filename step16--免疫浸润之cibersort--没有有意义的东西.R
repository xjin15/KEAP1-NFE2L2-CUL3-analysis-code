
# 数据准备 --------------------------------------------------------------------

rm(list = ls())
p_load(tidyverse, export, pheatmap,limma)
fpkm <- read_csv("outdata/mrna_fpkm_id_finish.txt",col_names = T,) %>% column_to_rownames(names(.)[1])
fpkm1 <- apply(fpkm, 2, FUN = function(x){2^x - 1}) # 去log化
str(fpkm1)
write.table(fpkm1,file = "cibersort/DATA.txt",sep = "\t",col.names = T,row.names = T )
# 写出fpkm1文件，做表达矩阵用

# write.table(fpkm,file = "cibersort/DATA.txt",sep = "\t",col.names = F,row.names = F,fileEncoding="UTF16")
# fpkm1 <- read.table(file = "cibersort/DATA.txt",header = T, sep = "\t", row.names = 1,check.names = F)

# sig <- read.table(file = "cibersort/LM22.txt", header = T, sep = "\t", row.names = 1, 
#                   check.names = F,fileEncoding="UTF16") # 这段代码是为了修改原cibersort出错而改的




# Cibersort
source("cibersort/Cibersort.R")

# Define LM22 file
LM22.file <- "cibersort/LM22.txt"
exp.file <- "cibersort/DATA.txt"

TME.results <-  CIBERSORT(LM22.file, exp.file, perm = 1000, QN = TRUE) 
#其中nperm给的是置换的次数，
# QN如果是芯片设置为T，如果是测序就设置为F，测序数据最好是TPM
# 前面是22种免疫细胞的在每个病人中的占比，后面还有三列代表p值，相关性，以及RMSE（均方根误差）
# TME.results <- t(TME.results)
# 每一行一个样本，每一列一种细胞，总共有22种细胞，
# 这里的数值代表的是免疫细胞所占的比例，
# 比如CD8+ T Cell在第一个肿瘤样本中是0.2282，
# 那就代表着在该样本中CD8+ T Cell占总的免疫细胞的22.82%，
# 所以，CIBERSORT输出的是一个比例，22种免疫细胞的比例加起来就等于1。
# # output CIBERSORT results
re <- TME.results[,-(23:25)]

write.table(TME.results, "TME.results.output.txt", 
            sep = "\t", row.names = T, col.names = T, quote = F)

# tme结果作热图 -----------------------------------------------------------------
TME.results <- read.table("TME.results.output.txt",
                          sep = "\t",row.names = 1,
                          header = T,check.names = F) %>% t()
  
re <- TME.results[,-(23:25)]

library(pheatmap)
pheatmap(re,show_rownames = F,cluster_rows = F) 
# 热图丑得不堪入目，要去掉在大部分样本中都为0的免疫细胞
k <- apply(re,2,function(x) {sum(x == 0) < nrow(TME.results)/2})
table(k)

re2 <- as.data.frame(t(re[,k]))

allid <- read.csv("outdata/allid_493.txt") %>% .[,2]
phe3 <- read.csv("outdata/phe_finish.csv",row.names = 1)
rownames(phe3) <-  NULL
phe3 <- phe3[match(x = allid, table = phe3$sample), ]
rownames(phe3) <- phe3$sample
phe3_anno <- phe3 %>% dplyr::select(group, age, gender, smoke_group,ajcc_T, ajcc_N, radiotherapy)

an <- phe3_anno

# an = data.frame(group = Group,
#                 row.names = colnames(exp))
pheatmap(re2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))


# 直方图每个患者的免疫细胞比例 ---------------------------------------------------------------------
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

dat <- re %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
dat$group <- an$group[match(dat$Sample, rownames(an))]



ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell_type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))

  

# 箱线图展示免疫细胞之间的比较 -------------------------------------------------------------------
ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))


# 加上顺序
a = dat %>% 
  group_by(Cell_type) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cell_type)

dat$Cell_type = factor(dat$Cell_type,levels = a)

ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))

######比较两组##########

# dat$Group = ifelse(as.numeric(str_sub(dat$Sample,14,15))<10,"tumor","normal")
library(ggpubr)
ggplot(dat,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),method = "kruskal.test")

graph2ppt(file = "output/plots/cibersort22cells",append = T,aspect = 1)


# 用现成的cibersort -----------------------------------------------------------
data <- read_csv("cibersort/cibersortscore.csv")
data <- data %>% 
  filter(CancerType == "LUAD") %>% 
  select(-c(25:27))
data <- data %>% 
  mutate(SampleID = substring(data$SampleID,first = 1,last = 16)) 
data$SampleID <- str_replace_all(string = data$SampleID,pattern = "\\.",replacement = "-")
data1 <- data[data$SampleID %in% allid, ] 
names <- data1[,1] %>% as.data.frame() %>% .[,1]
data1 <- data1[, -c(1:2)] %>% as.matrix()
rownames(data1) <- names
data1 <- avereps(data1)
data1 <- data1[allid, ]
k <- apply(data1,2,function(x) {sum(x == 0) < nrow(data1)/2})
table(k)

data2 <- as.data.frame(t(data1[,k]))

pheatmap(data2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

dat2 <- data1 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
dat2$group <- an$group[match(dat2$Sample, rownames(an))]

ggplot(dat2,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell_type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))


ggplot(dat2,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))
graph2ppt(file = "output/plots/cibersorts22cells",append = T,)

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
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),hide.ns = T,method = "kruskal.test")

graph2ppt("output/plots/cibersorts22cells",append = T,aspectr = 2)

compare_means(Proportion ~ group,group.by = "Cell_type", data = dat)


# 只选有意义的作图
dat2sig <- dat2 %>% 
  filter(Cell_type %in% c("Dendritic.cells.resting", 
                          "Mast.cells.resting","Monocytes",
                          "T.cells.CD4.memory.resting",
                          "T.cells.follicular.helper",
                          "T.cells.regulatory..Tregs."))
  
ggplot(dat2sig,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=0,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),hide.ns = T,method = "kruskal.test")

graph2ppt(file = "output/plots/cibersorts22cells",append = T,)



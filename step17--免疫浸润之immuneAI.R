
# 数据准备 --------------------------------------------------------------------

rm(list = ls())
p_load(tidyverse, export, pheatmap,limma)

# 用现成的immuneai -----------------------------------------------------------
data <- read.table(file = "cibersort/immucellAI_luad.txt",header = T)
allid <- read.csv("outdata/allid_493.txt") %>% .[,2]
data <- data %>% rownames_to_column("SampleID")
data <- data %>% 
  mutate(SampleID = substring(data$SampleID,first = 1,last = 16)) 
data$SampleID <- str_replace_all(string = data$SampleID,pattern = "\\.",replacement = "-")
names <- intersect(allid, data$SampleID)
data <- data %>% column_to_rownames("SampleID")
setdiff(allid,names) #"TCGA-44-3917-01A" 第9属于突变组
data1 <- data[names, -ncol(data1)]

# #########免疫浸润分数对比 箱线图、小提琴图###########
scdata <- data1 %>% 
  mutate( group = ifelse(rownames(data1) %in% allid[1:111], "Mut","Wild")) %>% 
  rownames_to_column("SampleID") %>% 
  select("SampleID","group",InfiltrationScore)
str(scdata)


library(ggpubr)
ggboxplot(scdata,x="group",y="InfiltrationScore",fill = "group",size = 0.1,
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

## 小提琴图
p <- ggviolin(scdata,x="group",y="InfiltrationScore",fill = "group",size = 0.1,
              palette = c("#00CCFF","#FF3333"),add = "boxplot",position = position_dodge(1), 
              add.params = list(group = "group"))+
  xlab("")+
  ylab("")+
  # scale_y_continuous(limits =c(-200, 1300),breaks = seq(-200,1300,500))+ 
  stat_compare_means(aes(group=group),label = "p.format", hide.ns = T,
                     bracket.size = 20,size = 9)+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "Immune Infiltration Score")+ #Mutation Load Plot
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0)) 
p+  theme(legend.text = element_text(size = 14, face = "bold"), 
          legend.title = element_text(size = 14, face = "bold")) +
    labs(x = NULL, y = NULL)
graph2ppt(file = "output/plots/immuneAI",append = T,asp = 1,width = 12)


######################## 24中细胞的免疫成分图#############
k <- apply(data1,2,function(x) {sum(x == 0) < nrow(data1)/2})
table(k)

data2 <- as.data.frame(t(data1[,k]))
phe3 <- read.csv("outdata/phe_finish.csv",row.names = 1)
rownames(phe3) <-  NULL
phe3 <- phe3[match(x = allid, table = phe3$sample), ]
rownames(phe3) <- phe3$sample
phe3_anno <- phe3 %>% dplyr::select(group, age, gender, smoke_group,ajcc_T, ajcc_N, radiotherapy)

an <- phe3_anno
pheatmap(data2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

dat2 <- data1 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
dat2$group <- an$group[match(dat2$Sample, rownames(an))]

library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

ggplot(dat2,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell_type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(24))


ggplot(dat2,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(24))
graph2ppt(file = "output/plots/immuneAI",append = T)

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
  scale_fill_manual(values = mypalette(24))
graph2ppt(file = "output/plots/immuneAI",append = T)

library(ggpubr)
ggplot(dat2,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),hide.ns = T,method = "kruskal.test")

graph2ppt(file = "output/plots/immuneAI",append = T)

compare_means(Proportion ~ group,group.by = "Cell_type", data = dat2)


# 只选有意义的作图
dat2sig <- dat2 %>% 
  filter(Cell_type %in% c("NK","Macrophage","DC","CD4_T","Th17",
                          "MAIT","Tr1","Neutrophil","CD8_naive",
                          "NKT","Tem"))
ggplot(dat2sig,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=0,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),hide.ns = T,method = "kruskal.test")

graph2ppt(file = "output/plots/immuneAI",append = T)



#

# 用immuneAI的分析结果 ----------------------------------------------------------
data <- read.table(file = "cibersort/immucellAI/ImmuCellAI_icb_result.txt",header = T)
allid <- read.csv("outdata/allid_493.txt") %>% .[,2]
data <- data %>% rownames_to_column("SampleID")
data <- data %>% 
  mutate(SampleID = substring(data$SampleID,first = 1,last = 16)) 
names <- intersect(allid, data$SampleID)
data <- data %>% column_to_rownames("SampleID")
data1 <- data[names, -c(1,ncol(data))]
scdata <- data %>% 
  mutate( group = ifelse(rownames(data1) %in% allid[1:111], "Mut","Wild")) %>% 
  rownames_to_column("SampleID") %>% 
  select("SampleID","group",Response,InfiltrationScore)
str(scdata)
table(scdata$group)
# 展示两个组别的免疫治疗的预测的成功可能性###########
tab = as.data.frame(table(scdata$group,scdata$Response))
tab$Var1=factor(tab$Var1,levels=c("Mut","Wild"))
fix(tab)
colnames(tab) <- c("Group","ICB Response","Number")
library(ggpubr)
ggbarplot(tab,x="Group",y= "Number",fill="ICB Response", size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6))
graph2ppt(file = "output/plots/immuneAI.pptx", append = T)
graph2ppt(file = "output/plots/immuneAI",append = T,asp = 1,width = 12)



library(ggpubr)
ggboxplot(scdata,x="group",y="InfiltrationScore",fill = "group",size = 0.1,
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

# 小提琴图
p <- ggviolin(scdata,x="group",y="InfiltrationScore",fill = "group",size = 0.1,
              palette = c("#00CCFF","#FF3333"),add = "boxplot",position = position_dodge(1), 
              add.params = list(group = "group"))+
  xlab("")+
  ylab("")+
  # scale_y_continuous(limits =c(-200, 1300),breaks = seq(-200,1300,500))+ 
  stat_compare_means(aes(group=group),label = "p.format", hide.ns = T,
                     bracket.size = 20,size = 9)+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "Immune Infiltration Score")+ #Mutation Load Plot
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0)) 
p+  theme(legend.text = element_text(size = 14, face = "bold"), 
          legend.title = element_text(size = 14, face = "bold")) +
  labs(x = NULL, y = NULL)
graph2ppt(file = "output/plots/immuneAI",append = T,asp = 1,width = 12)


#### 24种细胞的分布和对比#####
k <- apply(data1,2,function(x) {sum(x == 0) < nrow(data1)/2})
table(k)

data2 <- as.data.frame(t(data1[,k]))
phe3 <- read.csv("outdata/phe_finish.csv",row.names = 1)
rownames(phe3) <-  NULL
phe3 <- phe3[match(x = allid, table = phe3$sample), ]
rownames(phe3) <- phe3$sample
phe3_anno <- phe3 %>% dplyr::select(group, age, gender, smoke_group,ajcc_T, ajcc_N, radiotherapy)

an <- phe3_anno
pheatmap(data2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

dat2 <- data1 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
dat2$group <- an$group[match(dat2$Sample, rownames(an))]

library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

ggplot(dat2,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell_type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(24))


ggplot(dat2,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(24))
graph2ppt(file = "output/plots/immuneAI",append = T)

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
  scale_fill_manual(values = mypalette(24))
graph2ppt(file = "output/plots/immuneAI",append = T)

library(ggpubr)
ggplot(dat2,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),hide.ns = T,method = "kruskal.test")

graph2ppt(file = "output/plots/immuneAI",append = T)

compare_means(Proportion ~ group,group.by = "Cell_type", data = dat2)


# 只选有意义的作图
dat2sig <- dat2 %>% 
  filter(Cell_type %in% c("Macrophage","DC","Tfh","Th17","NK",
                          "Neutrophil","iTreg","CD4_T","Tr1","Th1",
                          "B_cell","Gamma_delta","nTreg","CD4_naive",
                          "NKT","Effector_memory"))
psig <- ggplot(dat2sig,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_classic() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=0,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),hide.ns = T,method = "kruskal.test")

graph2ppt(file = "output/plots/immuneAI",append = T)

######## 比较CD8+和nTreg以及iTreg 的比值######

dataratio <- data %>% 
  select(CD8_T,nTreg,iTreg) %>% 
  mutate(ratioN = CD8_T / nTreg,
         ratioI = CD8_T / iTreg,
         ratio = CD8_T / (nTreg + iTreg)  )
dataratio$group <- ifelse(rownames(dataratio) %in% allid[1:111], "Mut","Wild")


dataratio1 <- pivot_longer(dataratio,cols = c("ratioN","ratioI","ratio"))


ggviolin(x="name",y="value",fill = "group",data = dataratio1,size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(fill = "group"))+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits =c(0, 3),breaks = seq(0,3,0.5))+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = F,
                     bracket.size = 20,label.x = 1.35,label.y = 1.3,size = 15)+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "Mutation Load Plot")+
  theme(plot.title = element_text(size = 55, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))
ggboxplot(dataratio1,x="name",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  scale_y_continuous(limits = c(0, 3),breaks = seq(0,3,0.5))+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6)) 

ggviolin(x="group",y="ratio",fill = "group",data = dataratio,size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(fill = "group"))+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits =c(0, 3),breaks = seq(0,3,0.5))+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = F,
                     bracket.size = 20,label.x = 1.35,label.y = 1.3,size = 15)+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "Mutation Load Plot")+
  theme(plot.title = element_text(size = 55, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))

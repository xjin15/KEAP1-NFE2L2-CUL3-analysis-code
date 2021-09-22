############### 导入mirna_mature的数据并做整理构建表达矩阵 ##############
rm(list = ls())
p_load(limma,edgeR, ggplot2,ggthemes, ggpubr, tidyverse,ggrepel)
############# mirna 数据加载##########
mirna <- read.csv(file = "data/miRNA_HiSeq_gene", 
                  header = T, sep = "\t", stringsAsFactors = F, 
                  row.names = 1,check.names = F)
allid <- read.csv(file = "outdata/allid_493.txt", row.names = 1) %>% .[, 1]
# TCGA的名称.改成-
# colnames(mirna_stemloop) <- gsub('\\.', '-', colnames(mirna_stemloop))
mirna_mature <- mirna
# 确认没有重复的列名（样本名称）
stopifnot({length(colnames(mirna_mature))==length(unique(colnames(mirna_mature)))})
# 只保留fpkm样品中结尾01（初发癌）的样品 564减少到510
mirna_mature <- mirna_mature[,substring(text = colnames(mirna_mature), 
                          first = 14) == "01"]
colnames(mirna_mature) <- paste0(colnames(mirna_mature), "A")
dim(mirna_mature) # 去之前有mirna_mature 2228个,448个样本 
# 去掉在所有样品都不表达的基因（mirna_mature=0）
# 先转换成array 矩阵而不 是dataframe
rt <- as.matrix(mirna_mature) 
rt <- replace_na(rt,0)
rt <- rt[rowSums(rt) > 0, ]  
dim(rt) # 去以后剩mirna_mature 2221个
### 挑选出需要的sample
intersample <- intersect(allid,colnames(rt))
rt1 <-rt[ , intersample] # 426个样本
length(intersect(intersample,allid[1:111])) #前96个是mut
length(intersect(intersample,allid[112:493]))# 后330是wild

############## id转换############
name <- rownames(rt1)
p_load(miRBaseVersions.db)
items <- select(miRBaseVersions.db,
                keys = name,
                keytype = "MIMAT",
                columns = c("ACCESSION","NAME","VERSION"))
id_name <- items[items$VERSION == 20.0, c("ACCESSION","NAME")]
id_name <- id_name[match(x=rownames(rt1), table = id_name$ACCESSION), ]
rt2 <- cbind(rt1,id_name) %>% dplyr::select(ACCESSION, NAME, everything())
stopifnot({rownames(rt2)==rt2[,1]})
rownames(rt1) <- rt2[,2]
write.csv(rt1, file = "outdata/miRNA_rpm.txt", row.names = T)
############   用limma做差异分析      ###########
dt <- read.csv(file = "outdata/miRNA_rpm.txt", header = T,
               check.names = F, 
               row.names = 1,
               stringsAsFactors = F)
rt <- as.matrix(dt)
######## 创建分组信息
{
  group_list <- c(rep(1,96), rep(0, 330))
  group <- factor(group_list, 
                  levels = c(0,1), 
                  labels = c("Wild","Mut"))
  group
  design <- model.matrix(~0+group)
  rownames(design) <- colnames(dt)
  colnames(design) <- levels(group)
  design
  contrast.matrix <- makeContrasts(Mut - Wild,levels = design)
  contrast.matrix
  }

#######进行limma 差异分析######

  {
    fit <- lmFit(rt,design)
    fit1 <- contrasts.fit(fit, contrast.matrix)
    fit1 <- eBayes(fit1)
    qqt(fit1$t, df = fit1$df.prior+fit1$df.residual, pch = 16, cex = 0.2)
    abline(0,1) #QQ图看正态分布?
    all_diff <- topTable(fit1, adjust.method = 'fdr',coef=1,p.value = 1,lfc <- log(1,2),number = Inf,sort.by = 'logFC')
    head(all_diff)
    # 所有差异加一列ID改成gene名字
    all_diff$ID <- rownames(all_diff)
    # 加一列表示logP
    all_diff$logP<- -log10(all_diff$adj.P.Val)
    #保存为csv
    write.csv(all_diff,file="outdata/miRNA_all_diff.csv")
  
}

##############火山图绘制####################
{ 
  ## 得出所有的差异基因
  all_diff<-read.csv("outdata/miRNA_all_diff.csv",row.names = 1)
  all_diff$Group <- NULL
  all_diff$group = "not-significant"
  #将adj.P.Val<0.05,logFC>0.5的基因设置为显著上调基因
  #将adj.P.Val<0.05,logFC<0.5的基因设置为显著下调基因
  all_diff$group[which((all_diff$adj.P.Val<0.05) & (all_diff$logFC> 0.5))] = "up-regulated"
  all_diff$group[which((all_diff$adj.P.Val<0.05) & (all_diff$logFC< -0.5))] = "down-regulated"
  #查看上调和下调基因的数目
  table(all_diff$group)
  all_diff$group<-as.factor(all_diff$group)
  all_diff$logP<- -log10(all_diff$adj.P.Val)
  #火山图，加上颜色
  ggscatter(all_diff,x="logFC",y="logP",
            color = "group",
            palette = c("green","gray","red"),
            size = 1.5)+theme_base()
  #再加上辅助线
  p <- ggscatter(all_diff,x="logFC",y="logP",
                 color = "group",
                 palette = c("green","gray","red"),
                 size = 1.5)+theme_base() + 
    geom_hline(yintercept = 1.30,linetype = "dashed") + 
    geom_vline(xintercept = c(-0.5,0.5),linetype = "dashed")
  p
  # dev.size("px")
  # ggsave(p,filename = "output/KEAP1_diffgene_number_volcanoplot.pdf")      ##,width = ,height = )
}
#加上排名前十的基因
{
  all_diff$label = ""
  #对差异基因的p值由小到大排序，学会一个order算法！
  all_diff<-all_diff[order(all_diff$adj.P.Val),]
  #
  all_diff$X <- rownames(all_diff)
  all_diff$X <- NULL
  #高表达的基因中，选择adj.P.val最小的10个
  up.mirs<-head(all_diff$ID[which(all_diff$group == "up-regulated")],10)
  #低表达的基因中，选择adj.P.val最小的10个
  down.mirs<-head(all_diff$ID[which(all_diff$group == "down-regulated")],10)
  #以上两步合并加入label中
  all_diff.top10.mirs <- c(as.character(up.mirs),as.character(down.mirs))
  all_diff$label[match(all_diff.top10.mirs,all_diff$ID)] <- all_diff.top10.mirs
}

# 画火山图,repel=T表示字体不会重叠#   一运行就直接Rsession破灭
{
  # ggscatter(all_diff,x="logFC",y="logP",
  #           color = "group",
  #           palette = c("#00468BCC","gray","#ED0000CC"),
  #           label = all_diff$label, #repel=T
  #           label.select = all_diff.top10.genes,
  #           font.label = 10,
  #           size = 2,
  #           alpha = 0.7)+
  #   geom_hline(yintercept = 1.30,linetype = "dashed") +
  #   geom_vline(xintercept = c(-0.5,0.5),linetype = "dashed")+
  #   theme_minimal()+
  #   theme(axis.text.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
  #         axis.text.y = element_text(size = 15),
  #         axis.title.x = element_text(size = 15),
  #         axis.title.y = element_text(size = 15),legend.position = "right")
  # 
  #   
  p <- ggplot(data = all_diff, 
              aes(x = logFC, y = logP)) +
    geom_point(alpha=0.7, size=2, aes(color=group) ) +
    scale_color_manual(values = c("#00468BCC","gray","#ED0000CC"))+
    geom_vline(xintercept=c(-0.5,0.5),lty="dashed",col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty="dashed",col="black",lwd=0.8) +
    geom_text_repel(aes(label = label), box.padding = unit(1, "lines"), 
                    point.padding = unit(1, "lines"), show.legend = F, 
                    segment.color = 'black', size = 3,max.overlaps = 30)+
    theme_minimal()+
    theme(axis.text.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.position = "right"
    )
  p
  library(export)
  graph2ppt(file="output/plots/miRNA_volcanoplot.pptx")
  # ggsave(file="output/plots/miRNA_volcanoplot.pdf", p, width = 6, height = 6)
  
}


###### 绘制热图#####
#data presentation 数据准备好了
{
  DEG <- read.csv("outdata/miRNA_all_diff.csv")
  DEG <- DEG[,c(2,3,7)] #要求DEG有4列，geneid、FC、regulation、pval
  DEG <- DEG[is.finite(DEG$logFC),]#去掉无限小数的行
  DEG <- DEG[abs(DEG$logFC)>0.5 & DEG$adj.P.Val<0.05,] # 筛选logFC大于1.5的
  names(DEG) <- c("genes","fold_change","p_value")
  DEG$regulation <- "up"
  DEG$regulation[DEG$fold_change<0] <- "down"
}

#pheatmap 做热图
{
  library(scales)
  library(pheatmap)
  dd <- read.csv(file = "outdata/miRNA_rpm.txt", header = T,
                 check.names = F, ## 不检查列名！
                 row.names = 1,
                 stringsAsFactors = F) # dd表示总的表达矩阵
  DEG1<- DEG[order(abs(DEG$fold_change),decreasing = T),]#按照FC排序
  DEG1 <- DEG1[1:19,]#只取前40个基因，要求DEG有4列，geneid、FC、regulation、pval
  dd1=dd[rownames(dd) %in% DEG1$genes, ]# 在表达矩阵中选取这些基因
  dd2=apply(dd1,1,rescale)         ##归一化
  dd2=t(dd2)                     ##转置回来
  pheatmap(dd2,cluster_rows = T,cluster_cols = F,show_colnames = F)  # 做热图
  graph2ppt(file="output/plots/miRNA_top40genes_heatmap.pptx")
  # pdf(file="output/Top_40genes.pdf",width = 3,height = 5.5)
  # pheatmap(dd2,cutree_rows = 2,cutree_cols = 2,fontsize_row = 8 )
  # dev.off()
}




#################用R包预测miRNA的靶基因，然后对靶基因进行富集分析#################
# BiocManager::install("multiMiR",ask = F,update = F)
library(multiMiR)
# db.ver = multimir_dbInfoVersions()
# db.ver[,1:3]

# The default is to search validated interactions in human
example1 <- get_multimir(mirna = 'hsa-miR-18a-3p', summary = TRUE)
names(example1)
# Check which types of associations were returned
table(example1@data$type)

# Detailed information of the validated miRNA-target interaction
head(example1@data)
dim(example1@data)

# Which interactions are supported by Luciferase assay?
example1@data[grep("Luciferase", example1@data[, "experiment"]), ]
example1@summary[example1@summary[,"target_symbol"] == "KRAS",]
summary(example1)
#查询多个mirna的靶基因
multimir_results <- get_multimir(org = "hsa",
                                 mirna   = all_diff.top10.mirs,
                                  table   = 'validated',
                                 summary = T,
                                )
table(multimir_results@data$type)
dim(multimir_results@data)
head(multimir_results@data)
multimir_results@data[grep("Luciferase", multimir_results@data[, "experiment"]), ]

example1 <- get_multimir(org = "hsa",
                mirna   = up.mirs,
                target = down.genes,
                table   = 'all',
                summary = T ,
                predicted.cutoff = 10,
                predicted.cutoff.type = "p",
                use.tibble = T)
table(example1@data$type)
result <- select(example1, keytype = "type", keys = "validated", columns = columns(example1))
unique_pairs <- 
  result[!duplicated(result[, c("mature_mirna_id", "target_entrez")]), ]

result


example2 <- get_multimir(org = "hsa",
                         mirna   = down.mirs,
                         target = up.genes,
                         table   = 'all',
                         summary = T ,
                         predicted.cutoff = 10,
                         predicted.cutoff.type = "p",
                         use.tibble = T)
table(example2@data$type)
result2 <- select(example2, keytype = "type", keys = "validated", columns = columns(example2))
unique_pairs2 <- 
  result2[!duplicated(result2[, c("mature_mirna_id", "target_entrez")]), ]

result2

# ############富集分析

#enrichment analysis 富集分析：可以用R语言的clusterprofiler或者用david网站

{
  #   if (!requireNamespace("BiocManager", quietly = TRUE))
  #     install.packages("BiocManager")
  #   
  #   BiocManager::install("clusterProfiler")
  #   if (!requireNamespace("BiocManager", quietly = TRUE))
  #     install.packages("BiocManager")
  #   
  #   BiocManager::install("org.Hs.eg.db")
  #   BiocManager::install("DO.db")
  ##安装clusterProfiler,org.Hs.eg.db,
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(cowplot)
  keytypes(org.Hs.eg.db)
  
  
  gene <- bitr(DEG$genes, fromType ="SYMBOL",##把差异表达基因名称转换成geneID和SYMBOL
               toType = "ENTREZID", 
               OrgDb = org.Hs.eg.db)
  geneList <- bitr(rownames(dd), fromType = "SYMBOL",##dd里面是所有的基因
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  ego <- enrichGO(gene           = gene$ENTREZID,
                  universe       = names(geneList$ENTREZID),
                  OrgDb          = org.Hs.eg.db,
                  ont            = "BP", 
                  pAdjustMethod  = "BH",
                  pvalueCutoff   = 0.01,
                  qvalueCutoff   = 0.05,
                  readable       = T )              ##GO分析
  
  kk <- enrichKEGG(gene          = gene$ENTREZID,
                   organism      = 'hsa',
                   pvalueCutoff  = 0.05,
                   pAdjustMethod = "BH", 
                   qvalueCutoff  = 0.05)            ##KEGG分析  
  
  p1 <- dotplot(object = ego, showCategory = 10, orderBy = "x") + ggtitle("dotplot for GOBP")
  p2 <- dotplot(object = kk,  showCategory = 10, orderBy = "x") + ggtitle("dotplot for KEGG")
  p1
  graph2ppt(file="output/plots/mrna_go_plot.pptx")
  p2
  graph2ppt(file="output/plots/mrna_kegg_plot.pptx")
  
  # ggsave(pp,filename = "output/plots/DEG_enrichment.pdf",width = 12,height = 3.8)
  write.table(DEG,file = "output/DEG_keap1.xls",sep = "\t",quote = F,row.names = F)##输出差异表达基因
}
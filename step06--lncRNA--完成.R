#############lncRNA提取####################
rm(list = ls())
library(rtracklayer)
library(tidyverse)
p_load(limma,edgeR,ggplot2,ggthemes,ggpubr,ggrepel,export)

AnnoData = import('data/gencode.v22.long_noncoding_RNAs.gtf')
index = which(AnnoData$type == 'gene')
Target = data.frame(Ensembl_ID = AnnoData$gene_id[index], Symbol = AnnoData$gene_name[index], Biotype = AnnoData$gene_type[index])
Target$Ensembl_ID = gsub('\\..*', '', Target$Ensembl_ID)

rt <- read.csv(file = "outdata/mrna_fpkm_finish.txt", header = T,
               check.names = F, ## 不检查列名！
               row.names = 1,
               stringsAsFactors = F)
# 把data.frame转换成matrix,因为这样就允许行名重复了
rt <- as.matrix(rt)            # 57502个基因
rt[1:3,1:3]
# 去掉不表达的基因
rt <- rt[rowSums(rt) > 0, ]    # 57481个基因  
# 基因ID去掉点号后面的版本标识，准备转换ID。
rownames(rt) <- substring(text = rownames(rt), first = 1, last = 15)
# 确认没重复的基因ENSEMBL ID
stopifnot(
  {
    duplicated(rownames(rt))== F
  }
)
# 查看行和列数
dim(rt) #57481 493(mut111,wild382)
common = intersect(Target$Ensembl_ID, rownames(rt))#####15900中的15337个lncrna
lnc =rt[common, ]
rownames(lnc) <- Target$Symbol[match(x = rownames(lnc),table = Target$Ensembl_ID)]
# 重复的基因名字进行合并且取均值
dim(lnc) # 15337
lnc <- avereps(lnc) #15328
lnc[1:5,1:2]
write.csv(x = lnc, file = "outdata/lncRNA_fpkm_id_finish.txt")

####################### 正式开始limma分析####################
rm(list = ls())
dt <- read.csv(file = "outdata/lncRNA_fpkm_id_finish.txt", header = T,
               check.names = F, 
               row.names = 1,
               stringsAsFactors = F)
rt <- as.matrix(dt)


########### 设计limma比较分组矩阵################
allid <- read.csv("outdata/allid_493.txt") %>% .[,2]
group_list <- c(rep(1,111), rep(0, 382))
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


#######进行limma 差异分析######
{
  fit <- lmFit(rt,design)
  fit1 <- contrasts.fit(fit, contrast.matrix)
  fit1 <- eBayes(fit1)
  qqt(fit1$t, df = fit1$df.prior+fit1$df.residual, pch = 16, cex = 0.2)
  abline(0,1) #QQ图看正态分布?
  all_diff <- topTable(fit1, adjust.method = 'fdr',coef=1,p.value = 1,lfc <- log(1,2),number = Inf,sort.by = 'logFC')
  all_diff$ID <- rownames(all_diff)   # 所有差异加一列ID改成gene名字
  all_diff$logP<- -log10(all_diff$adj.P.Val)   # 加一列表示logP
  head(all_diff)
  write.csv(all_diff,file="outdata/lncRNA_all_diff.csv")   #保存为csv 
}





##############火山图绘制####################
{ 
  all_diff<-read.csv("outdata/lncrna_all_diff.csv",row.names = 1)  ## 读入所有的差异基因
  all_diff$group = "not-significant"  
  #将adj.P.Val<0.05,logFC>0.5的基因设置为显著上调基因
  #将adj.P.Val<0.05,logFC<0.5的基因设置为显著下调基因
  all_diff$group[which((all_diff$adj.P.Val<0.05) & (all_diff$logFC> 0.5))] = "up-regulated"
  all_diff$group[which((all_diff$adj.P.Val<0.05) & (all_diff$logFC< -0.5))] = "down-regulated"
  table(all_diff$group)      #查看上调和下调基因的数目
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
  up.lncs<-head(all_diff$ID[which(all_diff$group == "up-regulated")],10)
  #低表达的基因中，选择adj.P.val最小的10个
  down.lncs<-head(all_diff$ID[which(all_diff$group == "down-regulated")],10)
  #以上两步合并加入label中
  all_diff.top10.lncs <- c(as.character(up.lncs),as.character(down.lncs))
  all_diff$label[match(all_diff.top10.lncs,all_diff$ID)] <- all_diff.top10.lncs
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
  graph2ppt(file="output/plots/lncrna_volcanoplot.pptx",width = 8, height = 6)
  # ggsave(file="output/plots/lncrna_volcanoplot.pdf", p, width = 6, height = 6)
  
}


###### 绘制热图#####
#data presentation 数据准备好了
{
  DEG <- read.csv("outdata/lncrna_all_diff.csv")
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
  dd <- read.csv(file = "outdata/lncrna_fpkm_id_finish.txt", header = T,
                 check.names = F, ## 不检查列名！
                 row.names = 1,
                 stringsAsFactors = F) # dd表示总的表达矩阵
  DEG1<- DEG[order(abs(DEG$fold_change),decreasing = T),]#按照FC排序
  DEG1 <- DEG1[1:40,]#只取前40个基因，要求DEG有4列，geneid、FC、regulation、pval
  dd1=dd[rownames(dd) %in% DEG1$genes, ]# 在表达矩阵中选取这些基因
  dd2=apply(dd1,1,rescale)         ##归一化
  dd2=t(dd2)                     ##转置回来
  pheatmap(dd2,cluster_rows = T,cluster_cols = F,show_colnames = F)  # 做热图
  graph2ppt(file="output/plots/lncrna_top40genes_heatmap.pptx")
  # pdf(file="output/Top_40genes.pdf",width = 3,height = 5.5)
  # pheatmap(dd2,cutree_rows = 2,cutree_cols = 2,fontsize_row = 8 )
  # dev.off()
}

# ############富集分析做不了， 因为lncrna不是编码蛋白的

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
  graph2ppt(file="output/plots/lncrna_go_plot.pptx")
  p2
  graph2ppt(file="output/plots/lncrna_kegg_plot.pptx")
  
  # ggsave(pp,filename = "output/plots/DEG_enrichment.pdf",width = 12,height = 3.8)
  
}
write.table(DEG,file = "output/DEGlnc_keap1.xls",sep = "\t",quote = F,row.names = F)##输出差异表达基因

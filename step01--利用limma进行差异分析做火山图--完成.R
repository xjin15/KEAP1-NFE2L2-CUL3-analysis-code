#####准备工作####
rm(list = ls())
# 设定阈值
folfchange <- 1
padj <- 0.05
p_load(limma,edgeR,ggplot2,ggthemes,ggpubr,ggrepel,export,data.table,tidyverse)
# library(limma)
# library(edgeR)
# library(ggplot2)
# library(ggthemes)
# library(ggpubr)
# library(ggrepel)
# library(export)
rt <- fread(file = "outdata/mrna_fpkm_finish.txt",check.names = F,header = T)
rt <- rt %>% column_to_rownames("V1")

# rt <- read.csv(file = "outdata/mrna_fpkm_finish.txt", header = T,
#                 check.names = F, ## 不检查列名！
#                row.names = 1,
#                 stringsAsFactors = F)
fpkm <- rt
# 把data.frame转换成matrix,因为这样就允许行名重复了
rt <- as.matrix(rt)            # 57502个基因
rt[1:3,1:3]
# 去掉不表达的基因
rt <- rt[rowSums(rt) > 0, ]    # 57481个基因  

# 基因ID转换###
{
gtf_gene <- read.csv(file = "output/gtfmRNA22.txt", header = T,sep = "\t")
  ## 选择proteincoding的蛋白质
gtf_gene_protein <- gtf_gene %>% 
  dplyr::rename( Ensembl_ID = gene_id,
          Symbol  = gene_name,
          Biotype = gene_type  ) %>% 
  dplyr::filter(Biotype %in% c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", 
                               "IG_M_gene", "IG_V_gene", "IG_Z_gene", 
                               "nonsense_mediated_decay", "nontranslating_CDS", 
                               "non_stop_decay", "polymorphic_pseudogene", 
                               "protein_coding", "TR_C_gene", "TR_D_gene", "TR_gene", 
                               "TR_J_gene", "TR_V_gene"))
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
matchgene <- intersect(rownames(rt), gtf_gene_protein$Ensembl_ID)
mrna <- rt[matchgene, ]### 20039
# 可以转换ID
rownames(mrna) <- gtf_gene_protein$Symbol[match(x = rownames(mrna),table = gtf_gene_protein$Ensembl_ID)]
# 重复的基因名字进行合并且取均值
mrna <- avereps(mrna) #19971 
}
dim(mrna)
# write.table(x = rt, file = "outdata/mrna_fpkm_id_finish.txt", 
#             sep = "\t", col.names = T, row.names = T, 
#             quote = F)
write.csv(x = mrna, file = "outdata/mrna_fpkm_id_finish.txt")
write.csv(x=colnames(mrna), file = "outdata/allid_493.txt")
####################### 正式开始limma分析####################
# group=c(rep("mut",111),rep("wild",382))               ####按照癌症和正常样品数目修改
# design <- model.matrix(~group)    
dt <- read.csv(file = "outdata/mrna_fpkm_id_finish.txt", header = T,
               check.names = F, 
               row.names = 1,
               stringsAsFactors = F)
dt <- read_csv(file = "outdata/mrna_fpkm_id_finish.txt",col_names = T) %>% column_to_rownames('X1' ) %>% as.matrix()
rt <- as.matrix(dt)
 # data <- dt
# data<-as.matrix(data)
# y <-DGEList(counts=data,group=group)             ####转化成R擅长处理的格式，构建基因表达列表
# y <- calcNormFactors(y)                           ####＃标准化数据／归一化/计算样本内标准化因子
# y <- estimateCommonDisp(y)                      ####计算普通的离散度
# y <- estimateTagwiseDisp(y)                       ####计算基因间范围内的离散度
# et <- exactTest(y,pair = c("mut","wild"))              ####进行精确检验
# topTags(et)                                      ####输出排名靠前的差异表达基因信息
# ordered_tags <- topTags(et, n=100000)               
# top100<-topTags(et,100)  
# write.table(top100,file="top100.txt",sep="\t",quote=F,col.names=F)
# v <- voom(exprSet,design,normalize="quantile")
## 下面的代码如果你不感兴趣不需要运行，免得误导你
## 就是看看normalization前面的数据分布差异
# png("RAWvsNORM.png")
# 
# exprSet_new=v$E
# par(cex = 0.7)
# n.sample=ncol(exprSet)
# if(n.sample>40) par(cex = 0.5)
# cols <- rainbow(n.sample*1.2)
# par(mfrow=c(2,2))
# boxplot(exprSet, col = cols,main="expression value",las=2)
# boxplot(exprSet_new, col = cols,main="expression value",las=2)
# hist(exprSet)
# hist(exprSet_new)
#dev.off()

########### 设计limma比较分组矩阵################
{
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
write.csv(all_diff,file="outdata/mrna_all_diff.csv")
}





##############火山图绘制####################
{ 
## 得出所有的差异基因
all_diff<-read.csv("outdata/mrna_all_diff.csv",row.names = 1)
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
up.genes<-head(all_diff$ID[which(all_diff$group == "up-regulated")],10)
#低表达的基因中，选择adj.P.val最小的10个
down.genes<-head(all_diff$ID[which(all_diff$group == "down-regulated")],10)
#以上两步合并加入label中
all_diff.top10.genes <- c(as.character(up.genes),as.character(down.genes))
all_diff$label[match(all_diff.top10.genes,all_diff$ID)] <- all_diff.top10.genes
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
graph2ppt(file="output/plots/mrna_volcanoplot.pptx", width = 7, aspectr = 1.5)
# ggsave(file="output/plots/mrna_volcanoplot.pdf", p, width = 6, height = 6)

}


###### 绘制热图#####
#data presentation 数据准备好了
{
  DEG <- read.csv("outdata/mrna_all_diff.csv")
  DEG <- DEG[,c(2,3,7)] #要求DEG有4列，geneid、FC、regulation、pval
  DEG <- DEG[is.finite(DEG$logFC),]#去掉无限小数的行
  DEG <- DEG[abs(DEG$logFC)>0.5 & DEG$adj.P.Val<0.05,] # 筛选logFC大于0.5的
  names(DEG) <- c("genes","fold_change","p_value")
  DEG$regulation <- "up"
  DEG$regulation[DEG$fold_change<0] <- "down"
  table(DEG$regulation)
}
#bar plot条形图
{
  library(ggplot2)
  tab=as.data.frame(table(DEG$regulation))
  tab$Var1=factor(tab$Var1,levels = c("up","down"))
  p <- ggplot(tab,aes(x=Var1,y=Freq,label=Freq,fill=Var1))+geom_bar(stat="identity")  ##ggplot的使用
  p=p+theme_classic(8)+xlab("diffrential expression")+ylab("number of genes")+ggtitle("KEAP1cluste_Mut_diffgenes")
  p=p+theme(plot.title = element_text(hjust = 0.5))
  p=p+geom_text(position = position_dodge(0.9),vjust=0,size=3)+ylim(0,max(tab$Freq)*1.1)
  p
  graph2ppt(file="output/plots/KEAP1cluster_Mut_diffgenes")
  # dev.size("px")
   ggsave(p,filename = "output/KEAP1cluste_Mut_diffgenes_barplot.pdf")      ##,width = ,height = )
  # 
}

#pheatmap 做热图
{
  library(scales)
  library(pheatmap)
  # dd <- read.csv(file = "outdata/mrna_fpkm_id_finish.txt", header = T,
  #                        check.names = F, ## 不检查列名！
  #                        row.names = 1,
  #                        stringsAsFactors = F) # dd表示总的表达矩阵
  dd <- read_csv(file = "outdata/mrna_fpkm_id_finish.txt") %>% column_to_rownames("X1")
  DEG1<- DEG[order(abs(DEG$fold_change),decreasing = T),]#按照FC排序
  DEG1 <- DEG1[1:40,]#只取前40个基因，要求DEG有4列，geneid、FC、regulation、pval
  dd1=dd[rownames(dd) %in% DEG1$genes, ]# 在表达矩阵中选取这些基因
  dd2=apply(dd1,1,rescale)         ##归一化
  dd2=t(dd2)                     ##转置回来
  pheatmap(dd2,cluster_rows = T,cluster_cols = F,show_colnames = F)  # 做热图
  
#############给热图加上注释信息，也就是survival信息###############
  phe3 <- read.csv("outdata/phe_finish.csv",row.names = 1)
  rownames(phe3) <-  NULL
  phe3 <- phe3[match(x = allid, table = phe3$sample), ]
  rownames(phe3) <- phe3$sample
  phe3_anno <- phe3 %>% dplyr::select(group, age, gender, smoke_group,ajcc_T, ajcc_N, radiotherapy)
  pheatmap(dd2,cluster_rows = T,cluster_cols = F,show_colnames = F,annotation_col = phe3_anno)
  
  # annotation_col = data.frame( CellType = factor(rep(c("CT1", "CT2"), 5)), Time = 1:5 )
  # rownames(annotation_col) = paste("Test", 1:10, sep = "")
  # annotation_row = data.frame( GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6))))
  # rownames(annotation_row) = paste("Gene", 1:20, sep = "")
  graph2ppt(file="output/plots/mrna_top40genes_heatmap.pptx", width = 7, aspectr = 1)
  # pdf(file="output/Top_40genes.pdf",width = 3,height = 5.5)
  # pheatmap(dd2,cutree_rows = 2,cutree_cols = 2,fontsize_row = 8 )
  # dev.off()
}

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
  # library(clusterProfiler)
  # library(org.Hs.eg.db)
  # library(cowplot)
  p_load(clusterProfiler,org.Hs.eg.db,cowplot)
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
  
  p1 <- dotplot(object = ego, showCategory = 5, orderBy = "x") + ggtitle("dotplot for GOBP")
  p2 <- dotplot(object = kk,  showCategory = 5, orderBy = "x") + ggtitle("dotplot for KEGG")
  plot_grid(p1, p2, labels = c("A", "B"),ncol = 1)
  plot_grid(p1,p2,ncol = 1)
  p_load(GOSemSim,enrichplot)
  d <- godata('org.Hs.eg.db', ont="BP")
  ego2 <- pairwise_termsim(ego, method="Wang", semData = d)
  emapplot(ego2,showCategory=10) 
  graph2ppt(file="output/plots/mrna_go_plot.pptx", aspectr=1 , append = T)
  upsetplot(ego,5)
  graph2ppt(file="output/plots/mrna_go_plot.pptx", aspectr=1.5,append = T)
  graph2pdf(file="output/plots/mrna_go_plot.pdf", width = 15,)
  
p1
graph2ppt(file="output/plots/mrna_go_plot.pptx", append = T)
graph2pdf(file="output/plots/mrna_go_plot.pdf",width = 20,)


p2
graph2ppt(file="output/plots/mrna_kegg_plot.pptx",width = 10,append = T,aspectr=1)

  # ggsave(pp,filename = "output/plots/DEG_enrichment.pdf",width = 12,height = 3.8)
}
write.table(DEG,file = "output/DEGmrna_keap1.xls",sep = "\t",quote = F,row.names = F)##输出差异表达基因

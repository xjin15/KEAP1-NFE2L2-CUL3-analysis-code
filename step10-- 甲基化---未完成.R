###################甲基化分析流程##################

Methydata = read.table('data/',sep = '\t', header = T, quote = '', check.names = F,
                       fill = T,  comment.char = "!", stringsAsFactors = FALSE)
rownames(Methydata)<-Methydata[,1]
Methydata<-Methydata[,-1]
Methydata<-Methydata[,colnames(Methydata) %in% rownames(Methy_group)]
Methy_group<-Methy_group [rownames(Methy_group) %in% colnames(Methydata),]
Methydata<-Methydata[,match(rownames(Methy_group),colnames(Methydata))]
#这样就对齐了
Methydata<-as.matrix(Methydata)
myLoad$beta<-Methydata
myLoad$pd<-Methy_group
print("successful") 
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=4)
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$score_group,adjPVal = 0.05,adjust.method = "BH",arraytype = "450K")
diff_Methy<-myDMP

Methydata<-lapply(Methydata,BIrownames)

#对齐Methygroup和data
Methydata<-Methydata[,colnames(Methydata) %in% rownames(Methy_group)]
Methy_group<-Methy_group[rownames(Methy_group) %in% colnames(Methydata),]
Methydata<-Methydata[,match(rownames(Methy_group),colnames(Methydata))]
Methydata<-lapply(Methydata,as.matrix)
Methy_group<-lapply(Methy_group,as.matrix)
Methydata<-lapply(Methydata,function(df){
  rownames(df)<-Methyids$gene}
  #进行Methy差异分析
  #正式
  fit <- lmFit(Methydata,Methy_group)
  contrast.matrix <- makeContrasts(High - Low,levels=Methy_group)
  fit2 <- contrasts.fit(fit,contrast.matrix)
  fit2 <- eBayes(fit2)
  Methy_diff <- topTable(fit2, adjust.method = 'fdr',coef=1,p.value = 1,lfc <- log(1,2),number = 1881,sort.by = 'logFC')
  Methy_diff$logP<- -log10(Methy_diff$adj.P.Val)
  MethyDEG<-Methy_diff
  DEGhahaha<-DEGbar
  #从这里开始！绘制纵向条形图
  DEGtest<-DEG
  #保留adj.P<0.05者
  MethyDEG_test<-MethyDEG
  MethyDEG<-lapply(MethyDEG,function(df){
    df<-df%>%filter(df$adj.P.Val<0.05)
    #再保留|logFC|>0.2的,abs函数表示绝对值
    MethyDEG<-lapply(MethyDEG,function(df){
      df<-df%>%filter(abs(df$logFC)>0.2)
    })
    #加入high/low
    MethyDEG<-within(MethyDEG,{
      group<-NA
      group[logFC > 0 ] = 'high'
      group[logFC < 0 ] = 'low'
      MethyDEG<-lapply(MethyDEG,function(df){
        df$group<-as.factor(df$group)
        return(df)
library(ChAMP)
Methylation_high_to_low<-diff_Methy$high_to_low
Methylation_low_to_high<-diff_Methy$low_to_high
#删除delta<0.2的以及没有P>0.05的
Methylation_high_to_low<-lapply(Methylation_high_to_low,function(df){
df<-df%>%filter(df$adj.P.Val<0.05 & abs(df$deltaBeta)>0.2)
        })
        
Methylation_low_to_high<-lapply(Methylation_low_to_high,function(df){
df<-df%>%filter(df$adj.P.Val<0.05 & abs(df$deltaBeta)>0.2)
        })
#以上部分（主要是差异分析）需要用实验室电脑完成，下面的部分对电脑配置没有要求
#save.image("~/Documents/Research/增殖评分/DATA/Methy/Methy_filter_result.RData")
        
#调整正负号，以high组为基准，即deltaBeta=high-low,logFC=high/low
Methylation_high_to_low<-lapply(Methylation_high_to_low,function(df){
          df$logFC<-df$logFC*-1
          df$deltaBeta<-df$deltaBeta*-1
          return(df)
        })
        
#改变high/low列的位置
    for (i in 1:length(names(Methylation_low_to_high))){
      print(i)
      if(length(Methylation_low_to_high[[i]]$logFC) == 0) print("弹无虚发的大浦洞神枪手") else{
         Methylation_low_to_high[[i]]<-Methylation_low_to_high[[i]][,c(1,2,3,4,5,6,8,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)]
          }}
        
#双数据融合
    Methylation<-c(Methylation_high_to_low,Methylation_low_to_high)
    Methylation$OV<-Methylation$ACC
        
#######九.3 补绘制前面的条图：（四）中的关于甲基化的条图没画好，在这里重新画######
#加入high/low
Methylation_test<-Methylation
    for (i in 1:30){
    print(i)
    Methylation[[i]]<-within(Methylation[[i]],{
    group<-NA
    group[logFC > 0 ] = 'high'
    group[logFC < 0 ] = 'low'
    })
    }
Methylation<-lapply(Methylation,function(df){
    df$group<-as.factor(df$group)
    return(df)
        })
        
#合并数据
Methylationbar<-data.frame()
    for (i in 1:30){
      print(i)
  Methylationbar[2*i-1,1]<-names(Methylation)[i]
  Methylationbar[2*i-1,2]<-summary(Methylation[[i]]$group)[1]
  Methylationbar[2*i-1,3]<-'High'
  Methylationbar[,2][is.na(Methylationbar[,2])] <- 0
  Methylationbar[2*i,1]<-cancername[i]
  Methylationbar[2*i,2]<-summary(Methylation[[i]]$group)[2]
  Methylationbar[2*i,3]<-'Low'
  Methylationbar[,2][is.na(Methylationbar[,2])] <- 0
  Methylationbar[2*i-1,4]<-Methylationbar[2*i-1,2]+Methylationbar[2*i,2]
  Methylationbar[2*i,4]<-Methylationbar[2*i-1,2]+Methy
  Methylationbar[2*i,4]<-Methylationbar[2*i-1,2]+Methylationbar[2*i,2]
    }
colnames(Methylationbar)<-c('cancer','value','Group','sum')
Methylationbar$identify<-paste0(Methylationbar$cancer,Methylationbar$Group)
DEGbar$identify<-paste0(DEGbar$cancer,DEGbar$Group)
Methylationbar<-Methylationbar[match(DEGbar$identify,Methylationbar$identify),]

Methylationbar$rank<-DEGbar$sum

Methylationbar%>%mutate(cancer = fct_reorder(cancer, desc(rank)))%>% #此处表示按反向排列！管道连接符
  ggplot(aes(x=cancer, y=value,fill = Group)) + 
  geom_bar(stat="identity",alpha=0.7, width=0.8,position=position_dodge(0.8)) +
  theme_bw()+
  scale_fill_manual(values=c('High' = '#FF3333','Low' = '#0099FF'))+ #表示手动调节颜色！
  coord_flip() +
  xlab("") + ylab("")+
  theme(axis.text.x = element_text(size = 13, angle = 0, hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, angle = 0, hjust = 0, vjust = 0))







      




##############甲基化gsea图__________做不出来##############


rm(list = ls())
p_load(clusterProfiler,enrichplot,DOSE)
load( "outdata/champ_gsea.rds", verbose = T )
barplot(myGSEA) 
gseaplot(x = myGSEA, geneSetID = 1 )
# 不能用这些常规函数，因为数据结构是不一样的
### 下面是示例数据的gsea分析
data(geneList)
x <- gseDO(geneList) # 全部基因都用来做gsea，导入的数据必须是FDR值排好序的，gene列表
dotplot(x)
gseaplot(x, geneSetID=1)
gseaplot2(x, geneSetID=1)
gseaplot2(x, geneSetID=1:2)
dotplot(myGSEA)

barplot(myGSEA$DMP)

# barplot函数内部赋值处理
# layers = list(<environment>), scales = <environment>, mapping = structure(list(
#   size = ~Count, colour = ~p.adjust, x = ~GeneRatio, y = ~Description), class = "uneval"), 
# theme = structure(list(line = structure(list(colour = "black", 
#                                              size = 0.5, linetype = 1, lineend = "butt", arrow = FALSE, 
#                                         inherit.blank = TRUE)


test <- myGSEA   
# test <- reorder(test$DMP$Gene_List,test$DMP$fRep)
test1 <- test$DMP %>% 
  rename(Description = Gene_List,
         GeneRatio = fRep,
         p.adjust = adjPval,
         Count = nOVLAP) 
library(ggplot2)
# reorder(test$DMP$Gene_List,test$DMP$fRep)
# ggplot(data = test[1:10,], 
#        mapping = aes(size = Count, colour = p.adjust, x = GeneRatio, y = Description)
#        )+
#   geom_dotplot()
# Gene_list# MSigDB数据库中定义的基因集合
# nList# 每个基因集合包括的基因个数
# nRep# 基因集合的基因与所有输入的gene list 中overlap的基因个数
# fRep# overlap的基因的比例
# nOVLAP# 位于该基因集合下的基因与输入的gene list 中overlap的个数
# OR# 费舍尔精确检验的odds ratio
# Pvalue# 单尾fisher exact test检验的p值，具体代码如下
# listPV.v_2 <- t(apply(fisher.lm_2,1,function(x) unlist(fisher.test(matrix(x,2,2),alternative=”greater”)[c(1,3)])))
# adjPval# 多重假设检验校正之后的P值，默认采用”BH”方法
# Genes # gene symbol， 个数和nOVLAP相同
# 需要注意一点，对于fRep < 0.6 的基因集合，会过滤掉。官方是这样解释的 remove lists with less than 60% representation on array。
# test1 <- test %>% reorder(test$Description,test$fRep)
test1 <- test1 %>% 
  arrange(desc(GeneRatio))
# test1 <- test %>% 
#   arrange(desc(fRep))

ggplot(test1[1:10,],aes(x = GeneRatio,y = Description ))+
  geom_point(aes(color = p.adjust,
                 size = Count))+
  scale_color_gradient(low = "red", high = "blue")+
  xlab("Fold Enrichment")+
  theme_bw()+
  #edit legends
  guides(    #reverse color order (higher value on top)
    color = guide_colorbar(reverse = TRUE))
#reverse size order (higher diameter on top) 
#size = guide_legend(reverse = TRUE))
test2 <- test$DMR %>% 
  rename(Description = Gene_List,
         GeneRatio = fRep,
         p.adjust = adjPval,
         Count = nOVLAP) 
test2 <- test2 %>% 
  arrange(desc(GeneRatio))
# test1 <- test %>% 
#   arrange(desc(fRep))

ggplot(test2[1:10,],aes(x = GeneRatio,y = Description ))+
  geom_point(aes(color = p.adjust,
                 size = Count))+
  scale_color_gradient(low = "red", high = "blue")+
  xlab("Fold Enrichment")+
  theme_bw()+
  #edit legends
  guides(
    #reverse color order (higher value on top)
    color = guide_colorbar(reverse = TRUE))







############

library(ggplot2)
library(ggthemes)
# tab<-read.table("leaf.ori.IDs.keggEnrich.tab.final.xls",header=T,sep="\t")
tab <- TEST
colnames(tab)
tab <- tab  %>% 
  rename(Term.Name = Gene_List,
         richFactor = OR)
tab<-tab[,c("Term.Name",#"MainClass",
            # "GeneHitsInSelectedSet",
            # "GeneHitsInBackground"
            "richFactor")]
head(tab)
# tab$richFactor<-tab$GeneHitsInSelectedSet/tab$GeneHitsInBackground

tab<-tab[order(tab$richFactor,decreasing = T),]
# tab<-tab[order(tab$MainClass,decreasing = F),]
tab <- tab[!is.infinite(tab$richFactor),]
tab$Term.Name<-factor(tab$Term.Name,levels=unique(as.character(tab$Term.Name)))
p<-ggplot(tab,aes(x=Term.Name,y=richFactor))
p+geom_bar(stat="identity",width=0.1)+
# +geom_point(aes(color=MainClass),size=10)+
  # geom_text(aes(label=GeneHitsInSelectedSet),alpha=I(0.8))+
  theme_bw() +
  theme(
    
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), # element_line(size = 0.8,color="darkgray"), # element_blank(),
    axis.line.x = element_line(colour = "black", size = 0.8),
    axis.line.y = element_line(colour = "black", size = 0.8),
    axis.ticks.x = element_line(size = 0.8),
    axis.ticks.y = element_line(size = 0.8),
    axis.text.x = element_text(
      angle = 90, hjust = 0, vjust = 0
    ),
    #legend.position="NA",
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    legend.background = element_rect(fill = "transparent"),
    strip.background = element_rect(
      colour = "white", fill = "white",
      size = 0.2
    ),
    strip.text.x = element_text(size = 14),
    strip.text.y = element_text(size = 14),
    
    text = element_text(
      size = 14, 
      #family = "arial",
      face = "bold"
    ),
    plot.title = element_text(
      size = 16, 
      #family = "arial",
      face = "bold"    ))+scale_color_pander()+xlab("KEGG Pathway")+ylab("Rich Factor")

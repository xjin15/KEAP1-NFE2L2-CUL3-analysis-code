library(pacman)
library(limma)
library(export)
library(dplyr)
allid <- read.csv(file = "outdata/allid_493.txt") %>% .[,2]
clndata <- read.csv("outdata/cln_xena.csv") %>% arrange(desc(group))
Methydata = read.table(file = "data/TCGA-LUAD.methylation450.tsv", 
                       sep = '\t', header = T, quote = '', fill = T, 
                       comment.char = "!", stringsAsFactors = FALSE)
rownames(Methydata)<-Methydata[,1]
Methydata<-Methydata[,-1]
colnames(Methydata) <- gsub(pattern = "\\.",replacement = "-",x = colnames(Methydata),fixed = F)
rt <- Methydata ##rt 为原始文件，做备份，不然重新读入要半小时。
# save.image(file = "methy.RData")


Methydata <- Methydata[ , colnames(Methydata) %in% clndata$sample ]
allid1 <- allid [allid %in% colnames(Methydata)]
clndata <- clndata[clndata$sample %in% colnames(Methydata), ]
Methydata <- Methydata[,match(clndata$sample,colnames(Methydata))]
#这样就对齐了
Methydata<-as.matrix(Methydata)
Methy_group <- clndata[, c(1,5)]
rownames(Methy_group) <- Methy_group[, 1] 
myLoad <- list()
Methydata1<-na.omit(Methydata) #### 去掉有NA的行
myLoad$beta<-Methydata1
myLoad$pd<-Methy_group
print("successful") 



library(ChAMP)
####chamP
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores = 8)  ###标准化
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$group,adjPVal = 0.05,adjust.method = "BH",arraytype = "450K")
diff_Methy<-myDMP
write.csv(diff_Methy, file = "outdata/methy_diff_champ.csv")
# Methydata<-lapply(Methydata, BIrownames)
Methydata1 <- na.omit(Methydata)
#对齐Methygroup和data
Methydata1<-Methydata1[ , colnames(Methydata1) %in% rownames(Methy_group)]
Methy_group<-Methy_group[rownames(Methy_group) %in% colnames(Methydata1),]
Methydata1<-Methydata1[, match(rownames(Methy_group),colnames(Methydata1))]
# Methydata1<-lapply(Methydata1,as.matrix)
# Methy_group<-lapply(Methy_group,as.matrix) 做泛癌研究需要的,下面一行是做循环的函数
# Methydata1<-lapply(Methydata1,function(df){
#   rownames(df)<-Methyids$gene}
  #进行Methy差异分析
  #正式
  library(limma)
  group_list <- Methy_group$group
  group <- factor(group_list, levels = c("wild","mut"), labels = c("Wild","Mut") )
  group
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  rownames(design) <- colnames(Methydata1)
  design
  fit <- lmFit(Methydata1,design)
  head(coef(fit))
  contrast.matrix <- makeContrasts(Mut - Wild, Wild - Mut, levels=design)
  contrast.matrix
  fit1 <- contrasts.fit(fit,contrast.matrix)
  fit1 <- eBayes(fit1)

  
  head(coef(fit1)) ### 查看相关系数
  Methy_diff <- topTable(fit1, adjust.method = 'fdr',coef = 1, p.value = 1, lfc = log(1,2), number = Inf,sort.by = 'logFC')
  Methy_diff$logP<- -log10(Methy_diff$adj.P.Val)
  MethyDEG<-Methy_diff
  # DEGhahaha<-DEGbar
  #从这里开始！绘制纵向条形图
  # DEGtest<-DEG ### 这里放差异基因，而不是差异甲基化
  #保留adj.P<0.05者
  MethyDEG_test<-MethyDEG
  # MethyDEG<-lapply(MethyDEG,function(df){
    df <- MethyDEG
    df<- df %>% filter(df$adj.P.Val<0.05)
    df<-df%>%filter(abs(df$logFC)> 0.2)
    df <- within(df,{
      group <- NA
      group[logFC > 0.2] <-  "high"
      group[logFC < -0.2] <-"low"
    })   #得到了差异甲基化的表，进行过滤了。下一步做热图吧
    table(df$group)
    # #再保留|logFC|>0.2的,abs函数表示绝对值
    # MethyDEG<-lapply(MethyDEG,function(df){
    #   df<-df%>%filter(abs(df$logFC)>0.2)
    # })
    #加入high/low
    # MethyDEG<-within(MethyDEG,{
    #   group<-NA
    #   group[logFC > 0 ] = 'high'
    #   group[logFC < 0 ] = 'low'
      # MethyDEG<-lapply(MethyDEG,function(df){
      #   df$group<-as.factor(df$group)
      #   return(df)
write.csv(df, file = "outdata/methy_diff_limma.csv")

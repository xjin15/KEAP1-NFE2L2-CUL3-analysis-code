########## cnv分析，最后得到差不多是cnv的箱线图或者小提琴图，和snv分析流程相似
####### ZYS's scripts for CNV analysis######
rm(list = ls())
p_load(TCGAbiolinks,dplyr,DT,tibble,SummarizedExperiment,reshape2,maftools,limma,export)
cln_survival <- read.csv("outdata/cln_xena.csv")

# TCGAbiolinks::getGDCprojects()$project_id    #  查看projectid有哪些
# TCGAbiolinks:::getProjectSummary("TCGA-LUAD")#  查看TCGA-LUAD的category有哪些
#### TCGAbiolink 下载cnv gistic score 数据 
{
query_cnv <- GDCquery(project = "TCGA-LUAD",
                      data.category = "Copy Number Variation",
                      data.type = "Gene Level Copy Number Scores",
                      access = "open"
                      )
GDCdownload(query = query_cnv)
cnv_bio <- GDCprepare(query = query_cnv)
}
#### 导入xena上面下载的gistic数据。

TCGAcnv_xena<-read.table("../maf分组/data/TCGA-LUAD.gistic.tsv",
                    sep = '\t',
                    header = T,
                    quote = '',
                    fill = T, 
                    check.names = F,
                    stringsAsFactors = FALSE)
# view以后对比看了一下，xena和tcga官网的的数据一样的。
# 整理数据 （还是以 TCGA biolink 的为准）
# 行名改成sample id 相当于进行基因ID转换
colnames(cnv_bio) <- substr(x = colnames(cnv_bio), start = 1,stop = 16)

names(TCGAcnv_xena)[1] <- "Gene.Symbol"
TCGAcnv_xena$Gene.Symbol<- base::gsub("\\..*", "", TCGAcnv_xena$Gene.Symbol)#删去点后的内容
cnv_bio$`Gene Symbol`<- base::gsub("\\..*", "", cnv_bio$`Gene Symbol`)#删去点后的内容


###############ID转换######################
{
  # genename
  # library(org.Hs.eg.db)
  # library(clusterProfiler)
  # keytypes(org.Hs.eg.db)  # ID转换会失败的
  # genename1 <- bitr(genename, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")#将基因名转换为ENTREZID格式
  # genename2 <- genename1[match(x = TCGAcnv$Gene.Symbol, table = genename1$ENSEMBL), ]
# 有部分ID转换失败，并且都是有数值的。
gtf_gene <- read.csv(file = "outdata//gtfmRNA22.txt", header = T,sep = "\t")
genename <- as.character(TCGAcnv_xena$Gene.Symbol) 
genename2 <- gtf_gene[match(x=genename, table = gtf_gene$gene_id), ]
genename2$gene_name %>% unique() %>% length() #有重复的symbol名称
TCGAcnv_xena$Gene.Symbol <- genename2$gene_name
}


#去除缺失值
# TCGAcnv=TCGAcnv %>% filter(is.na(Gene.Symbol)==FALSE)
# genename=genename %>% filter(is.na(SYMBOL)==FALSE)
cnv<-as.matrix(TCGAcnv_xena)
rownames(cnv)<-cnv[,1]
cnv<-cnv[,-1]
dim(cnv)
cnv <- avereps(cnv)
dim(cnv) # 19729 去重后19645
### 删除重复的样本，并且只保留后缀为01A的, 还把名字改成‘-’的
cnv <- cnv[, nchar(colnames(cnv)) == 16 ]
cnv <- cnv[, substring(text = colnames(cnv), 
                               first = 14, last = 16) == "01A"]
write.csv(cnv,file = "outdata/cnv_id_finish.txt") #整理完了写出去


#####################表达矩阵整理好了，接下来做分析 ################
#此时，将数据改成长格式
TCGAcnv_xena <- read.csv(file = "outdata/cnv_id_finish.txt", row.names = 1, check.names = F) %>% as.matrix()
TCGAcnv_xena<-t(TCGAcnv_xena)
TCGAcnvmelt<-as.data.frame(TCGAcnv_xena)
TCGAcnvmelt<-reshape2::melt(TCGAcnv_xena)
#去除0
class(TCGAcnvmelt$value)
TCGAcnvmelt$value<-as.character(TCGAcnvmelt$value)
TCGAcnvmelt$value<-as.numeric(TCGAcnvmelt$value)
TCGAcnvmelt <-  TCGAcnvmelt%>%filter(value != 0)
TCGAcnvmelt$CN[which(TCGAcnvmelt$value ==1)] = "Amp"
TCGAcnvmelt$CN[which(TCGAcnvmelt$value ==-1)] = "Del"
#重新排序
colnames(TCGAcnvmelt)<-c("Sample_name","Gene","Value", "CN")
TCGAcnvmelt<-TCGAcnvmelt[,c("Gene","Sample_name", "CN")]

write.csv(TCGAcnvmelt,file = "outdata/TCGAcnvmelt.csv")

TCGAcnvmelt<-read.csv("outdata/TCGAcnvmelt.csv",row.names = 1)
#整合入data
class(TCGAcnvmelt$Sample_name)
TCGAcnvmelt$Sample_name<-as.character(TCGAcnvmelt$Sample_name)#注意这里要改为文字格式，要不然下面会分成两类！


####################    mut组的cnv分析   #############
TCGAcnvmelt<-read.csv("outdata/TCGAcnvmelt.csv",row.names = 1)
#整合入data
class(TCGAcnvmelt$Sample_name)
TCGAcnvmelt$Sample_name<-as.character(TCGAcnvmelt$Sample_name)#注意这里要改为文字格式，要不然下面会分成两类！
TCGAdataMAF<-read.maf("outdata/TCGAdata_mut_MAF_maftools.maf")
#通过写出data再读回去的方式分亚组
dataA<-TCGAdataMAF@data
dataA<-as.data.frame(dataA)
# dataA$Tumor_Sample_Barcode<-gsub('-','.',dataA$Tumor_Sample_Barcode)
# dataA$Tumor_Sample_Barcode<-gsub('............$','',dataA$Tumor_Sample_Barcode)
dataA = dataA[dataA$Tumor_Sample_Barcode %in% cln_survival$sample,]

TCGAcnvmeltA = TCGAcnvmelt[TCGAcnvmelt$Sample_name %in% dataA$Tumor_Sample_Barcode,]

TCGAdataMAFA.plus.cn <- read.maf(dataA,cnTable = TCGAcnvmeltA,isTCGA = T)

#展示重点变量的总结信息
#Shows sample summry.
getSampleSummary(TCGAdataMAF)
SampleSummaryA<-getSampleSummary(TCGAdataMAFA.plus.cn)
#Shows gene summary.
getGeneSummary(TCGAdataMAF)
GeneSummary<-getGeneSummary(TCGAdataMAFA.plus.cn)
#Shows all fields in MAF
getFields(TCGAdataMAF)
Fields<-getFields(TCGAdataMAFA.plus.cn)
#shows clinical data associated with samples
getClinicalData(TCGAdataMAF)
ClinicalData<-getClinicalData(TCGAdataMAFA.plus.cn)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = TCGAdataMAFA.plus.cn, basename = 'outdata/TCGAdataMAFA(mut).plus.cn')
# plotmafSummary
#这东西不能有cnv在里面


TCGAdataMAFA.plus.cn<-read.maf(maf = 'outdata/TCGAdataMAFA(mut).plus.cn_maftools.maf')
plotmafSummary(maf = TCGAdataMAFA.plus.cn, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
graph2ppt(file = "output/plots/突变组maf全景图.pptx")
graph2pdf(file = "output/plots/突变组maf全景图.pdf")
#瀑布图
oncoplot(maf = TCGAdataMAFA.plus.cn, top = 20,genesToIgnore = "KEAP1")#,
graph2pdf(file = "output/plots/突变组top20oncoplot")
# genes = c("TP53","TTN","MUC16","CSMD3","RYR2","LRP1B","USH2A","ZFHX4","KRAS","FLG","SPTA1","XIRP2","NAV3","ZNF536","ANK2","CSMD1","FAT3","COL11A1","MUC17","PCDH15"))#显示特定基因 
#最终确定的基因(根据total的数目选择，或者根据与研究相关的基因进行选择)
# oncoplot(maf = TCGAdataMAFA.plus.cn, 
#          keepGeneOrder = T,
#          genes = c("TTN","TP53","CSMD3","RYR2","USH2A","MUC16","ZFHX4","SPTA1",
#                    "LRP1B","NAV3","KRAS","FLG","XIRP2","STK11","COL11A1",
#                    "APOB","TNR","RYR3","PCLO","FAT3","ADGRG4"))#显示特定基因 
# oncoplot(maf = TCGAdataMAFA.plus.cn, 
#          genes = c( "AKR1C1", "AKR1C2", "ZNF724P", "ZNF730", "TMPRSS7", 
#   "ZNF98", "SCN1A", "AKR1C3", "ZNF676", "SMARCA4", "SNX19", 
#   "RPSAP58", "SPICE1", "ZNF681", "CD200R1L", "ABHD10", "ATG3", 
#   "TSPAN16", "FAM181B", "PRSS55", "STK11", "DDIAS", "ZNF728"))
graph2ppt(file = "output/plots/突变组加上cnv的oncoplot.pptx")
graph2pdf(file = "output/plots/突变组加上cnv的oncoplot.pdf")


####################  wild组的cnv分析 ####################


TCGAdataMAF <- read.maf("outdata/TCGAdata_wild_MAF_maftools.maf")
#通过写出data再读回去的方式分亚组
dataB<-TCGAdataMAF@data
dataB<-as.data.frame(dataB)
# dataB$Tumor_Sample_Barcode<-gsub('a-','.',dataB$Tumor_Sample_Barcode)
# dataB$Tumor_Sample_Barcode<-gsub('............$','',dataB$Tumor_Sample_Barcode)
dataB <- dataB[dataB$Tumor_Sample_Barcode %in% cln_survival$sample,]

TCGAcnvmeltB = TCGAcnvmelt[TCGAcnvmelt$Sample_name %in% dataB$Tumor_Sample_Barcode,]

TCGAdataMAFB.plus.cn<-read.maf(dataB,cnTable = TCGAcnvmeltB)

#展示重点变量的总结信息
#Shows sample summry.
getSampleSummary(TCGAdataMAF)
SampleSummaryB<-getSampleSummary(TCGAdataMAFB.plus.cn)
#Shows gene summary.
getGeneSummary(TCGAdataMAF)
GeneSummary<-getGeneSummary(TCGAdataMAFB.plus.cn)
#Shows all fields in MAF
getFields(TCGAdataMAF)
Fields<-getFields(TCGAdataMAFB.plus.cn)
#shows clinical data associated with samples
getClinicalData(TCGAdataMAF)
ClinicalData<-getClinicalData(TCGAdataMAFB.plus.cn)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = TCGAdataMAFB.plus.cn, basename = 'outdata/TCGAdataMAFB.plus.cn')
#plotmafSummary
#这东西不能有cnv在里面
TCGAdataMAFB.plus.cn<-read.maf(maf = 'outdata/TCGAdataMAFB.plus.cn_maftools.maf')
plotmafSummary(maf = TCGAdataMAFB.plus.cn, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
graph2ppt(file = "output/plots/野生组maf全景图.pptx")
graph2pdf(file = "output/plots/野生组maf全景图.pdf")
#瀑布图
oncoplot(maf = TCGAdataMAFB.plus.cn, top = 20)#显示特定基因
graph2pdf(file = "output/plots/野生组top20oncoplot")

# genes = c("TP53","TTN","MUC16","CSMD3","RYR2","LRP1B","USH2A","ZFHX4","KRAS","FLG","SPTA1","XIRP2","NAV3","ZNF536","ANK2","CSMD1","FAT3","COL11A1","MUC17","PCDH15")
#最终确定的基因
oncoplot(maf = TCGAdataMAFB.plus.cn, 
         keepGeneOrder = F,
         genes = 
         c( "AKR1C1", "AKR1C2", "ZNF724P", "ZNF730", "TMPRSS7",
                    "ZNF98", "SCN1A", "AKR1C3", "ZNF676", "SMARCA4", "SNX19",
                    "RPSAP58", "SPICE1", "ZNF681", "CD200R1L", "ABHD10", "ATG3",
                    "TSPAN16", "FAM181B", "PRSS55", "STK11", "DDIAS", "ZNF728"))#显示特定基因
graph2ppt(file = "output/plots/野生组加上cnv的oncoplot.pptx")
graph2pdf(file = "output/plots/野生组加上cnv的oncoplot.pdf")


################ 两个MAF相互比较 ############
compareMAF_AB <- mafCompare(m1 = TCGAdataMAFA.plus.cn, m2 = TCGAdataMAFB.plus.cn, 
                            m1Name = 'MUT', m2Name = 'WILD', minMut = 10, useCNV = T)#最小mut数默认为5的基因
forestPlot(mafCompareRes = compareMAF_AB, color = c('royalblue', 'maroon'), fdr = 0.3, geneFontSize = 0.8)
## 设置了fdr以后，设置P值就没用了
graph2ppt(file = "output/plots/野生vs突变8基因森林图.pptx") # 记得把KEAP1和 NFE2L2删除 

print(compareMAF_AB)
summary_compareMAF<-as.data.frame(compareMAF_AB$results)
write.csv(summary_compareMAF,"outdata/summary_compareMAF.csv")

# deSNP <- compareMAF_AB$results %>%
#   dplyr::arrange(., order(or,decreasing = T)) %>% 
#   dplyr::filter(is.finite(or)) %>% 
#   dplyr::select(Hugo_Symbol) %>% head( , n=25) %>% as.data.frame() %>% .[,1] 
# deSNP
##得到差异明显的前8个基因，但是这些基因是突变频率比较低的
deSNP <- c("STK11","EGFR","GRIN2B","SPEF2", "SNTG2","BRWD3","OR6N1","ADGRB1")
###这里我确定了10个基因，是前10突变的，看起来有差距的10个基因。
genes <- c("TTN","RYR2","CSMD3","USH2A","SPTA1","ZFHX4","NAV3","EGFR","SPEF2","SNTG2")
deSNP <- genes
## 画出瀑布图
oncoplot(maf = TCGAdataMAFA.plus.cn, 
         keepGeneOrder = T,
         genes = deSNP)
graph2ppt(file = "output/plots/突变组8低频基因的oncoplot.pptx")
graph2pdf(file = "output/plots/突变组8低频基因的oncoplot")

oncoplot(maf = TCGAdataMAFB.plus.cn, 
         keepGeneOrder = T,
         genes = deSNP)#显示特定基因
graph2ppt(file = "output/plots/野生组8低频基因的oncoplot")
graph2pdf(file = "output/plots/野生组8低频基因的oncoplot.pdf")

coOncoplot(genes = deSNP, m1 = TCGAdataMAFA.plus.cn, m2 = TCGAdataMAFB.plus.cn, m1Name = 'MUT', m2Name = 'WILD', removeNonMutated = TRUE)
graph2pdf(file = "output/plots/野生vs突变8低频基因的cooncoplot.pdf")

coBarplot(genes = deSNP, m1 = TCGAdataMAFA.plus.cn, m2 = TCGAdataMAFB.plus.cn, 
          m1Name = 'MUT', m2Name = 'WILD', yLims = c(30,30) )
graph2pdf(file = "output/plots/野生vs突变8低频基因的cobarplot.pdf")

#将比较结果可视化，绘制森林图 keap1肯定明显，因为这是我的分组条件
forestPlot(mafCompareRes = compareMAF_AB, pVal = 0.1, 
           color = c('royalblue', 'maroon'),
           fdr = 0.1,
           geneFontSize = 0.8)


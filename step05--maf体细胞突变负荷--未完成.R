#########体细胞突变负荷###########
############################# 1 数据准备##########################
rm(list = ls())

p_load(maftools, GenVisR, ggplot2, reshape2, data.table, ggpubr, export)
cln_survival <- read.csv("outdata/cln_xena.csv")
allid <- read.csv("outdata/allid_493.txt") %>% .[,2]

load("outdata/sample_group.Rdata")
TCGAdataMAF<-read.maf("data/LUAD.maf")

somaticInteractions(maf = TCGAdataMAF,top = 15,)
INTERACT <- somaticInteractions(maf = TCGAdataMAF, genes = c("KEAP1","EGFR","SPEF2","GRIN2B","SNTG2",
                                                 "RYR2","STK11","TP53","MUC16",
                                                 "TTN"), pvalue = c(0.05, 0.1))
# graph2ppt(file = "output/plots/8SNP相关性图",append= T)
# graph2pdf(file = "output/plots/8SNP相关性图")
# INTERACT

#通过写出data再读回去的方式分亚组
data <- TCGAdataMAF@data
data <- as.data.frame(data)

data$Tumor_Sample_Barcode <- gsub('-','.',data$Tumor_Sample_Barcode)
data$Tumor_Sample_Barcode <- gsub('............$','',data$Tumor_Sample_Barcode)
data$Tumor_Sample_Barcode <- gsub( '\\.','-', x = data$Tumor_Sample_Barcode)

mutname <- cln_survival$sample[which(cln_survival$group == "mut")]
wildname <- cln_survival$sample[which(cln_survival$group == "wild")]

datamut <- data[data$Tumor_Sample_Barcode %in% mutname, ]
datawild <- data[data$Tumor_Sample_Barcode %in% wildname, ]
# 看相关性
somaticInteractions(maf = TCGAdataMAF, top = 25, pvalue = c(0.05, 0.1))


##################### 2 突变组做瀑布图#######################
TCGAdataMAF<-read.maf(datamut)
somaticInteractions(maf = TCGAdataMAF, top = 25, pvalue = c(0.05, 0.1))


#展示重点变量的总结信息
#Shows sample summry.
getSampleSummary(TCGAdataMAF)
SampleSummaryA <- getSampleSummary(TCGAdataMAF)
#Shows gene summary.
getGeneSummary(TCGAdataMAF)
#Shows all fields in MAF
getFields(TCGAdataMAF)
#shows clinical data associated with samples
getClinicalData(TCGAdataMAF)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = TCGAdataMAF, basename = 'outdata/TCGAdata_mut_MAF')
# 瀑布图里面ignore KEAP1
oncoplot(maf = TCGAdataMAF, top = 20, genesToIgnore = "KEAP1")


MAF <- as.data.frame(TCGAdataMAF@data)
waterfall(MAF[!MAF$Hugo_Symbol=="KEAP1",], mainRecurCutoff=.23,#[0,1]范围内的数值，基于样本中特定基因的突变频率选择plot的基因。
          mainXlabel=F,#布尔型，是否显示X轴的sample name；默认为FALSE
          mainDropMut=F,#布尔型，是否从mutation type legend中去除不使用的”mutation type”，默认为FALSE
          plotMutBurden = T,#布尔型，是否在顶端展示mutBurden的sub-plot，默认为TRUE
          #plotSamples = ,#指定plot的样本
          #maxGenes=100,#指定最多可以plot的基因数
          rmvSilent = F,#是否去除沉默突变，默认为FALSE
          fileType = "MAF",#读入X的数据格式，分为：”MGI”, “MAF”, “Custom”
          out = "plot") #输出的类型,包含 “data”, “grob”, 和”plot”三种,默认为 “plot”
## 也可以选择自己想展示的基因来画图
waterfall(MAF, #mainRecurCutoff=.17,#[0,1]范围内的数值，基于样本中特定基因的突变频率选择plot的基因。
            mainXlabel=F,#布尔型，是否显示X轴的sample name；默认为FALSE
            mainDropMut=F,#布尔型，是否从mutation type legend中去除不使用的”mutation type”，默认为FALSE
            plotMutBurden = T,#布尔型，是否在顶端展示mutBurden的sub-plot，默认为TRUE
            #plotSamples = ,#指定plot的样本
            #maxGenes=100,#指定最多可以plot的基因数
            plotGenes = c("TTN","TP53","CSMD3", "USH2A","MUC16","ZFHX4","SPTA1",
                          "LRP1B","NAV3","KRAS","FLG","XIRP2","STK11","COL11A1",
                          "APOB","TNR","RYR3","PCLO","FAT3","ADGRG4"),
            # geneOrder = c("TTN","TP53","CSMD3", "USH2A","MUC16","ZFHX4","SPTA1",
            #               "LRP1B","NAV3","KRAS","FLG","XIRP2","STK11","COL11A1",
            #               "APOB","TNR","RYR3","PCLO","FAT3","ADGRG4"),
            rmvSilent = F,#是否去除沉默突变，默认为FALSE
            fileType = "MAF",#读入X的数据格式，分为：”MGI”, “MAF”, “Custom”
            out = "plot") #输出的类型,包含 “data”, “grob”, 和”plot”三种,默认为 “plot” 
# graph2ppt(file="output/plots/突变组瀑布图.pptx")
graph2pdf(file="output/plots/突变组瀑布图2", width = 12, aspectr = 0.7)
# 这行代码有问题
{waterfalldata <- waterfall (MAF, mainRecurCutoff=.17,#[0,1]范围内的数值，基于样本中特定基因的突变频率选择plot的基因。
                             mainXlabel=F,#布尔型，是否显示X轴的sample name；默认为FALSE
                             mainDropMut=F,#布尔型，是否从mutation type legend中去除不使用的”mutation type”，默认为FALSE
                             plotMutBurden = T,#布尔型，是否在顶端展示mutBurden的sub-plot，默认为TRUE
                             #plotSamples = ,#指定plot的样本
                             #maxGenes=100,#指定最多可以plot的基因数
                             rmvSilent = F)#是否去除沉默突变，默认为FALSE
                                }

        
################### 3 野生组做瀑布图###############
TCGAdataMAF<-read.maf(datawild)
somaticInteractions(maf = TCGAdataMAF, genes = c("KEAP1","EGFR"), pvalue = c(0.05, 0.1))
OncogenicPathways(maf = TCGAdataMAF)

#展示重点变量的总结信息
#Shows sample summry.
getSampleSummary(TCGAdataMAF)
SampleSummaryB <- getSampleSummary(TCGAdataMAF)
#Shows gene summary.
getGeneSummary(TCGAdataMAF)
#Shows all fields in MAF
getFields(TCGAdataMAF)
#shows clinical data associated with samples
getClinicalData(TCGAdataMAF)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = TCGAdataMAF, basename = 'outdata/TCGAdata_wild_MAF')
oncoplot(maf = TCGAdataMAF, top = 20)

MAF<-as.data.frame(TCGAdataMAF@data)
waterfall(MAF, mainRecurCutoff=.17,#[0,1]范围内的数值，基于样本中特定基因的突变频率选择plot的基因。
          mainXlabel=F,#布尔型，是否显示X轴的sample name；默认为FALSE
          mainDropMut=F,#布尔型，是否从mutation type legend中去除不使用的”mutation type”，默认为FALSE
          plotMutBurden = T,#布尔型，是否在顶端展示mutBurden的sub-plot，默认为TRUE
          #plotSamples = ,#指定plot的样本
          #maxGenes=100,#指定最多可以plot的基因数
          rmvSilent = F,#是否去除沉默突变，默认为FALSE
          fileType = "MAF",#读入X的数据格式，分为：”MGI”, “MAF”, “Custom”
          out = "plot") #输出的类型,包含 “data”, “grob”, 和”plot”三种,默认为 “plot”
# graph2ppt(file="output/plots/野生组瀑布图.pptx")
 graph2pdf(file="output/plots/野生组瀑布图2", width = 12, aspectr = 0.7)


################ mutational burden  肿瘤突变负荷 ###########

TCGAtmb<-rbind(SampleSummaryA,SampleSummaryB) #只有111+381=492个sample
write.csv(TCGAtmb,file = "outdata/TCGA_LUAD_tmv.csv")
TCGAtmb<-read.csv("outdata/TCGA_LUAD_tmv.csv", row.names = 1)
cln_survival <- read.csv("outdata/cln_xena.csv")
# TCGAsanky<- read.csv("data/cln_xena.csv")

#### 给TCGAtmb加上分组情况
# TCGAsanky<-TCGAsanky[match(TCGAtmb$Tumor_Sample_Barcode,TCGAsanky$X),]

TCGAtmb$group <- cln_survival$group[match(x = TCGAtmb$Tumor_Sample_Barcode, 
                                          table = cln_survival$sample, nomatch = 0)]

#### 只保留突变负荷数和分组信息两列
rownames(TCGAtmb)<-TCGAtmb$Tumor_Sample_Barcode
TCGAtmb<-TCGAtmb[,-1]
TCGAtmb<-TCGAtmb[c(10,11)]
library(tidyverse)
TCGAtmbbox %>% filter(group == "mut") %>% summarise(mean = mean(value))
mutval <- TCGAtmbbox %>% filter(group == "mut") %>%select(value)  
wildval <- TCGAtmbbox %>% filter(group == "wild") %>%select(value)  
apply(mutval, 2, median)
apply(wildval, 2, median)
#变成长格式
TCGAtmbbox <- reshape2::melt(data = TCGAtmb,id.vars=c("group"))
TCGAtmbbox$group <- as.factor(TCGAtmbbox$group)
#箱线图绘制
ggplot(TCGAtmbbox)+  
  geom_boxplot(aes(x= variable,y= value,fill=group),
               width = 0.6,        #宽度
               position = position_dodge(0.8),
               outlier.size = 1, outlier.color = "black"#箱外点大小颜色
               )+
  scale_fill_manual(values = c("#FF3333","#00CCFF"),
                    breaks = c("mut","wild"),#分组的东西
                    labels = c("Mut","Wild"))+
  xlab("Mut  VS.  Wild")+
  ylab("Total Mutation Burden")+
  theme(axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 0.5))+
  scale_y_continuous(limits =c(0, 150),
                     breaks = seq(0,150,50))


####小提琴图绘制
#  scale_y_continuous(limits =c(0,1600),breaks = seq(0,1600,500))+ 

p <- ggviolin(x="variable",y="value",fill = "group",data = TCGAtmbbox,size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot",position = position_dodge(1), 
         add.params = list(group = "group"))+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits =c(-200, 1300),breaks = seq(-200,1300,500))+ 
  stat_compare_means(aes(group=group),label = "p.format", hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y = 1200,size = 9)+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "Tumor Mutation Burden")+ #Mutation Load Plot
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0)) 
p+  theme(legend.text = element_text(size = 14, 
    face = "bold"), legend.title = element_text(size = 14, 
    face = "bold")) +labs(x = NULL, y = NULL)


graph2ppt(file="output/plots/TMBvalue_violinplot.pptx", width = 8, aspectr = 1.5)
############ 更好看的箱线图
# TCGAtmbbox$value <- as.numeric(TCGAtmbbox$value)
# compare_means(value ~ group, data = TCGAtmbbox)
# ggboxplot(TCGAtmbbox,x="variable",y="value",fill = "group",size = 0.1,
#           palette = c("#00CCFF","#FF3333"))+
#   xlab("")+ylab("Log(Copy number)")+
#   theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 0.5, vjust = 0.5))+
#   scale_y_continuous(limits = c(-50, 200),breaks = seq(-50,200,50))+
#   stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
#                      bracket.size = 10)+
#   theme_gray()

ggviolin(x="variable",y="value",fill = "group",data = TCGAtmbbox,size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot",position = position_dodge(1), 
         add.params = list(group = "group"))+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits =c(-200, 1300),breaks = seq(-200,1300,500))+ 
  stat_compare_means(aes(group=group),label = "p.format", hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y = 1200,size = 9)+
  theme_cleveland()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "Tumor Mutation Burden")+ #Mutation Load Plot
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0)) +
  theme(legend.text = element_text(size = 14, 
                                    face = "bold"), 
        legend.title = element_text(size = 14,face = "bold")) +
  labs(x = NULL, y = NULL)
graph2ppt(file="output/plots/TMBvalue_violinplot.pptx", append=T, width = 8, aspectr = 1.5)
graph2pdf(file="output/plots/TMBvalue_violinplot", )


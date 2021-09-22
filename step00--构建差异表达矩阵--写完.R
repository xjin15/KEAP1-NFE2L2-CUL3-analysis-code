###### 整理数据
rm(list = ls())
library(limma)
library(dplyr)
library(tidyverse)
library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)
# 因为是根据三基因突变来分组,所以先用到了snv数据，之后还需要再用到snv。
snv <- read.csv(file = "data/TCGA-LUAD.mutect2_snv.tsv", 
             header = T, sep = "\t", stringsAsFactors = F)
# 导入mrna——fpkm数据
fpkm <- read.csv(file = "data/TCGA-LUAD.htseq_fpkm.tsv", 
             header = T, sep = "\t", stringsAsFactors = F, 
             row.names = 1, check.names = F)#第一列做行名
fpkm1 <- read_tsv(file = "data/TCGA-LUAD.htseq_fpkm.tsv", col_names = T) 
rownames(fpkm1) <- fpkm1$Ensembl_ID


# # TCGA的名称.改成-
# colnames(fpkm) <- gsub('\\.', '-', colnames(fpkm))
# 确认没有重复的列名（样本名称）
length(colnames(fpkm))==length(unique(colnames(fpkm)))

# 只保留fpkm样品中结尾01（初发癌）的样品585个
fpkm <- fpkm[,substring(text = colnames(fpkm), 
          first = 14, last = 16) == "01A"]
dim(fpkm) # 去之前有60483基因 510个样本
# 去掉在所有样品都不表达的基因（fpkm=0）
fpkm <- fpkm[rowSums(x = fpkm) > 0 , ] 
dim(fpkm) # 去以后剩57502基因 510

# 行和列进行转置 # 转置的意义不大，这里注释了 2021/1/24 by jinx
# fpkm_t <- t(fpkm)

# 导入xena的cln数据 sample样本，OS生存，OStime生存时间，X_PATIENT病人id
cln_xena <- read.csv(file = "data/LUAD_survival.csv") #738个样品
length(unique(cln_xena$X_PATIENT)) # 509 个病人

# snv加一列patientid
snv <- transform(snv, patientid = substr(x = Sample_ID, start = 1, stop = 12))
snv <- snv[,  c(1, 2, 13, 3:12)] # 变换一下列的顺序，利于查看

length(unique(snv$patientid)) # snv里面有567个病人

length(substring(text = snv$Sample_ID, first = 14, 
                 last = 16) == "01A")
length(substring(text = snv$Sample_ID, first = 14, 
                 last = 14) == "0")
length(snv$Sample_ID)
## 来自于567个病人的208180个snv样品都是01A结尾（原位癌）


# 筛选到有KEAP1,CUL3,NFE2L2突变的SAMPLEid

# patientfilter1 <- snv$patientid[snv$mutation == "1"]
patientfilter1 <- snv$Sample_ID[snv$gene == "KEAP1" | 
                                snv$gene == "CUL3"  | 
                                snv$gene == "NFE2L2"]

# 有KEAP1、CUL3、NFE2L2突变筛选出132个样本
patientfilter1 <- unique(patientfilter1)
#实际只有126个sample id，因为有些病人有多个基因突变。

# patientfilter2 <- snv$patientid[snv$mutation == "0"]
# patientfilter2 <- unique(patientfilter2)
## 这样子是筛不出来的，因为有KEAP1等基因突变的病人也有其他基因突变。
patientfilter2 <- unique(snv$Sample_ID[!(snv$Sample_ID %in% patientfilter1)])
# 从snv中减去有突变的病人id得到没有三个突变的病人的id

patientfilter3 <- intersect(patientfilter1, colnames(fpkm))
# 筛选在fpkm里面有突变的sampleid
patientfilter4 <- intersect(patientfilter2, colnames(fpkm))
# 在fpkm中没有keap1,cul3,nfe2l2突变的sampleid
patientfilter5 <- intersect(patientfilter3, cln_xena$sample)
# 在cln、fpkm、snv三个数据集都有的(KEAP1,CUL3,NFE2L2)突变sampleid
patientfilter6 <- intersect(patientfilter4, cln_xena$sample)
# 在cln、fpkm、snv三个数据集都有的(KEAP1,CUL3,NFE2L2)野生sampleid
mut111 <- fpkm[, patientfilter5]
wild382 <- fpkm[, patientfilter6]
fpkm_finish <- cbind(mut111, wild382)
# 输出mrna_fpkm_finish.txt文件
# write.table(x = fpkm_finish, file = "outdata/mrna_fpkm_finish.txt", 
#             sep = "\t", col.names = T, row.names = T, 
#             quote = F)
write.csv(x = fpkm_finish, file = "outdata/mrna_fpkm_finish.txt")



keap1id <- snv$Sample_ID[snv$gene == "KEAP1" ] %>% unique()
nfe2l2id <- snv$Sample_ID[snv$gene == "NFE2L2" ] %>% unique()
Cul3id <- snv$Sample_ID[snv$gene == "CUL3" ] %>% unique()
intersect(keap1id,nfe2l2id)
intersect(keap1id,Cul3id)
intersect(Cul3id,nfe2l2id)

data <- c(k=keap1id,n=nfe2l2id,c=Cul3id); duplicated(data) %>% table()

data1 <- unique(data)
data2 <- intersect(data1,colnames(fpkm))
allid <- intersect(data2,cln_xena$sample);allid %>% length()
anyDuplicated(allid)

kid <- intersect(allid,keap1id);length(kid)
nid <- intersect(allid,nfe2l2id);length(nid)
cid <- intersect(allid,Cul3id);length(cid)

intersect(kid,nid)
nid <- nid[nid != "TCGA-55-8302-01A"]
intersect(kid,nid)
## 88，14，9
allid <- c(kid,nid,cid,colnames(wild382))
sample_group <- data.frame(sample = allid,
                           group = factor(x = c(rep(1,88),rep(2,14),rep(3,9),rep(4,382)),
                                          labels = c('KEAP1','NFE2L2','CUL3','Wild'),
                                          levels = c(1:4)
                                          )
                           )

sample_group <- sample_group %>% 
  mutate(
    groupk = factor(group,labels = c('KEAP1','KEAP1-WT','KEAP1-WT','KEAP1-WT')),
    groupn = factor(group,labels = c('NFE2L2-WT','NFE2L2','NFE2L2-WT','NFE2L2-WT')),
    groupc = factor(group,labels = c('CUL3-WT','CUL3-WT','CUL3','CUL3-WT'),)
  )
## relevel 函数，把突变放在前面，factor的顺序不要乱
str(sample_group)
# levels(sample_group$groupn) <- c( "NFE2L2","NFE2L2-WT") 不行
sample_group$groupn <- relevel(x = sample_group$groupn,ref = "NFE2L2")
sample_group$groupc <- relevel(x = sample_group$groupc,ref = "CUL3")
str(sample_group)


save(sample_group,file = "outdata/sample_group.Rdata")

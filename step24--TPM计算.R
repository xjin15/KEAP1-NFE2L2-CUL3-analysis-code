# 首先模拟一个数据（之后的演示都基于这个数据mat）
set.seed(123)
mat=matrix(sample(200:300,36,replace=T),6,6)
rownames(mat)=paste0('gene',1:6)
colnames(mat)=paste0('sample',1:6)
mat[5,]=sample(1:10,6,replace = T)
mat[6,]=sample(80:100,6,replace = T)
mat
# sample1 sample2 sample3 sample4 sample5 sample6
# gene1     230     249     290     292     208     275
# gene2     278     242     268     298     282     214
# gene3     250     300     290     271     235     231
# gene4     213     213     256     225     277     206
# gene5      10       7       5       7       5       6
# gene6      81      84      87      91      92      97

# 手动计算gene1在sample1的cpm
colSums(mat)
# sample1 sample2 sample3 sample4 sample5 sample6 
# 1062    1095    1196    1184    1099    1029 
230/colSums(mat)[1]*10^6
# sample1 
# 216572.5

# edger的cpm函数自动计算
edgeR::cpm(mat)[1,1]
# [1] 216572.5
edgeR::cpm(mat) %>% colSums()
# sample1 sample2 sample3 sample4 sample5 sample6 
# 1e+06   1e+06   1e+06   1e+06   1e+06   1e+06


###### RPKM##### 
# 假设拿到基因长度
lengths <- c(1000,2000,500,1500,3000,2500)

# 各个样本文库大小
colSums(mat) 
# sample1 sample2 sample3 sample4 sample5 sample6 
# 1062    1095    1196    1184    1099    1029 

# RPKM需要先对文库标准化（以第一列为例演示）
x1 = mat[,1]*10^6/colSums(mat)[1]
# 再对长度标准化
x1*10^3/lengths
# gene1      gene2      gene3      gene4      gene5      gene6 
# 216572.505 130885.122 470809.793 133709.981   3138.732  30508.475

# 写成函数就是：
countToFpkm <- function(x){
  N <- sum(x)
  x*10^3*10^6/(N*lengths)
}
fpkm <- apply(mat,2,countToFpkm) 

fpkm[1:4,1:4]
# sample1  sample2  sample3  sample4
# gene1 216572.5 227397.3 242474.9 246621.6
# gene2 130885.1 110502.3 112040.1 125844.6
# gene3 470809.8 547945.2 484949.8 457770.3
# gene4 133710.0 129680.4 142697.9 126689.2
colSums(fpkm)  # FPKM值的每列之和就是每个样本的文库大小
# sample1   sample2   sample3   sample4   sample5   sample6 
# 985624.6 1048340.9 1012653.3  989639.6  948256.0  993326.9

###### TPM########
# 以第一列为例
# 先对长度标准化
x1 = mat[,1]*10^3/lengths
# 再对文库标准化
x1*10^6/sum(x1)
# gene1      gene2      gene3      gene4      gene5      gene6 
# 219731.227 132794.090 477676.581 135660.149   3184.511  30953.442

# 写成函数就是：
countToTpm <- function(x){
  normLen <- x*10^3/lengths
  normLen*10^6/sum(normLen)
}

tpm=apply(mat, 2, countToTpm)

tpm[1:4,1:4]
# sample1  sample2  sample3  sample4
# gene1 219731.2 216911.6 239445.1 249203.5
# gene2 132794.1 105406.8 110640.2 127162.0
# gene3 477676.6 522678.4 478890.3 462562.6
# gene4 135660.1 123700.6 140914.8 128015.5
colSums(tpm) # 每一列的长度之和都相等
# sample1 sample2 sample3 sample4 sample5 sample6 
# 1e+06   1e+06   1e+06   1e+06   1e+06   1e+06 


#### 指标四 TMM 适用于分组比较###
# 需要制定一个分组
group <- factor(c('c','c', 'c', 't', 't', 't'))
library(edgeR)

y <- DGEList(counts=mat, group=group)
# normalize for library size by cacluating scaling factor using TMM (default method)
y <- calcNormFactors(y)

# 可以看到每个样本根据各自文库大小，计算了norm.factors 
y
{
# An object of class "DGEList"
# $counts
# Sample1 Sample2 Sample3 Sample4 Sample5 Sample6
# gene1     276     213     250     217     262     233
# gene2     258     261     267     293     272     288
# gene3     286     228     257     294     240     300
# gene4     261     210     257     223     271     294
# gene5       4       5       9       8       1       9
# gene6      92      99      97      91      98      80
# 
# $samples
# group lib.size norm.factors
# Sample1     c     1177    0.9880091
# Sample2     c     1016    1.0158008
# Sample3     c     1137    0.9877519
# Sample4     t     1126    1.0290219
# Sample5     t     1144    0.9829446
# Sample6     t     1204    0.9973071
}

#Obtain the TMM values of your data 
# TMM和cpm不一样，但还是比较接近的。
cpm(y); cpm(mat)
# sample1    sample2   sample3    sample4    sample5    sample6
# gene1 209788.154 231597.913 241569.42 247519.357 195093.769 262820.562
# gene2 253570.030 225087.128 223243.47 252605.371 264502.130 204522.183
# gene3 228030.603 279033.630 241569.42 229718.307 220418.441 220769.272
# gene4 194282.073 198113.877 213247.49 190725.532 259812.376 196876.494
# gene5   9121.224   6510.785   4164.99   5933.683   4689.754   5734.267
# gene6  73881.915  78129.416  72470.83  77137.882  86291.475  92703.980

##### 指标5 DESeq median-of ratios
# 不管是DESeq还是DESeq2，校正方法都和TMM类似，
# 同样假设：大部分基因并非差异表达
# 构建一个样本比较矩阵
cond = data.frame(sample = colnames(mat), 
                  condition = group)
cond
# sample condition
# 1 sample1         c
# 2 sample2         c
# 3 sample3         c
# 4 sample4         t
# 5 sample5         t
# 6 sample6         t

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = mat, 
                              colData = cond, 
                              design = ~ condition)

dds
{
# class: DESeqDataSet 
# dim: 6 6 
# metadata(1): version
# assays(1): counts
# rownames(6): gene1 gene2 ... gene5 gene6
# rowData names(0):
#   colnames(6): sample1 sample2 ... sample5 sample6
# colData names(2): sample condition
}
dds <- estimateSizeFactors(dds)
y = counts(dds, normalized = TRUE)

y
# sample1    sample2   sample3    sample4    sample5    sample6
# gene1 244.59887 258.846568 272.35173 275.824071 215.221126 301.977645
# gene2 295.64560 251.569757 251.69057 281.491689 291.790180 234.993513
# gene3 265.86834 311.863335 272.35173 255.987408 243.158484 253.661221
# gene4 226.51983 221.422968 240.42084 212.535671 286.616596 226.208708
# gene5  10.63473   7.276811   4.69572   6.612221   5.173585   6.588603
# gene6  86.14134  87.321734  81.70552  85.958871  95.193960 106.515751

sizeFactors(dds)
# sample1   sample2   sample3   sample4   sample5   sample6 
# 0.9403150 0.9619598 1.0647995 1.0586458 0.9664479 0.9106634 

##### 指标6：GeTMM####兼具样本间和样本内的比较.
# Smid et al., 2018在BMC Bioinformatics发表了这个指标，特色是兼容样本内+样本间的比较

# calculate reads per Kbp of gene length (corrected for gene length)
# gene length is in bp in exppression dataset and converted to Kbp
rpk <- ( (mat*10^3 )/lengths)

y <- DGEList(counts=rpk, group=group)
y <- calcNormFactors(y)
# count per million read (normalized count)
norm_counts <- cpm(y)

norm_counts
# 小结论:
# 样本间比较用CPM；样本内比较TPM优于RPKM




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


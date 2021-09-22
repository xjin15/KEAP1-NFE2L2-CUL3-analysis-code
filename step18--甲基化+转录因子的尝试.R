##### 甲基化和转录因子的尝试####
# https://bioconductor.org/packages/devel/bioc/vignettes/ELMER/inst/doc/analysis_regulatory_tf.html
rm(list = ls())
p_load(MultiAssayExperiment,ELMER,ELMER.data,sesameData,tidyverse)

# step1--数据导入 --------------------------------------------------------------------
# 需要导入的文件 前两个必须
# 以探针为行，以样本为列的β值甲基化矩阵
# 以基因为为行，以样本为列并且标准化处理后的转录组矩阵（FPKM等）
# 3.甲基化样本名和转录组样本名对应的转化矩阵
# 4.样本本身的临床信息（metadata）
# 先筛选远端探针（distal probe）的id。
# 示范的时候只选取了第一条染色体的
allid <- read.csv("outdata/allid_493.txt") %>% .[,2]
distal.probes <- get.feature.probe(genome = "hg19", 
                                   met.platform = "450K", 
                                   rm.chr = paste0("chr",c("X","Y")))
#这里同样可以获取promoter区域内的探针

promoter.probes <- get.feature.probe(genome = "hg19", 
                                     met.platform = "450K", 
                                     TSS.range = list(upstream = 2000, downstream = 2000),
                                     promoter = T,
                                     rm.chr = paste0("chr",c("X","Y")))                          
head(distal.probes)  #返回的是Grange格式的
# 导入示例数据

mrna <- read_csv("outdata/mrna_fpkm_finish.txt")
mrna$X1 <- substring(mrna$X1,1,15)
mrna <- mrna %>% column_to_rownames("X1")
mrna <- as.matrix(mrna)
aa <- apply(mrna,1,sum) 
mrna <- mrna[rowSums(mrna) > 5, ]
mrna[1:5,1:5]

met <- data.table::fread("data/TCGA-LUAD.methylation450.tsv")
met <- met %>% column_to_rownames(names(met)[1])
met <- met %>% na.omit(met) %>% as.matrix()
met[1:5,1:5]
phe3 <- read.csv("outdata/phe_finish.csv",row.names = 1)
rownames(phe3) <-  NULL
phe3 <- phe3[match(x = allid, table = phe3$sample), ]
rownames(phe3) <- phe3$sample
phe3_anno <- phe3 %>% dplyr::select(group, age, gender, smoke_group,ajcc_T, ajcc_N, radiotherapy)

metadata <- phe3_anno
#开始创建MAE对象
mae <- createMAE(exp = mrna, 
                 met = met,
                 save = TRUE,
                 linearize.exp = TRUE,
                 save.filename = "mae.rda",#同时保存数据为mae.rda
                 filter.probes = distal.probes,
                 met.platform = "450K",
                 genome = "hg19",
                 TCGA = TRUE)
# 创建出MAE对象后我们可以对这个对象进行适当了解，
# 看一下我们如何获取MEA对象中相应的数据。
#获取甲基化表达谱数据
getMet(mae) %>% assay() %>% as.data.frame() %>% .[1:5,1:5]
#获取MEA对象转录组数据
getExp(mae) %>% assay() %>% as.data.frame() %>% .[1:5,1:5]
#获取MEA对象的metadata信息
colData(mae) %>% as.data.frame() %>% .[1:5,1:5]


# step2--筛选出差异甲基化探针 -------------------------------------------------------


# 我们通过比较病变组vs.正常组样本中所有远端探针（distal probe）的甲基化水平，
# 并筛选出adj.P value<0.01，并且两组之间的Δβ＞0.3的甲基化位点，
# 并且是在肿瘤组织低甲基化水平的位点。
# （1）supervised，即我们定义好的疾病组和正常组中所有样本进行分组，
# 从而进行差异分析筛选出在疾病组低甲基化远端探针。
# （2）unsupervised，即对甲基化探针在所有样本中甲基化水平进行排序，
# 选取最高甲基化水平20%以及最低甲基化水平20%进行差异分析，
# 从而得到在疾病组低甲基化远端探针。
# 这一类做法的原因是：疾病（尤其是肿瘤）具有明确的异质性，
# 通过提取20%的子集能够进一步找到疾病中特定的分子亚型。
sig.diff <- get.diff.meth(data = mae, 
                          group.col = "definition",
                          group1 =  "Primary solid Tumor",
                          group2 = "Solid Tissue Normal",
                          minSubgroupFrac = 0.2, 
                          sig.dif = 0.3,
                          diff.dir = "hypo", # 在group 1为低甲基化水平的差异探针
                          cores = 1, 
                          dir.out ="result", 
                          pvalue = 0.01)


#展示一下结果
head(sig.diff)
# 工作目录会生成一个result文件夹，里面是火山图


# step3--确定probe-gene对 ----------------------------------------------------
# 这一步我们将会把远端探针（distal probe）的甲基化水平和靶基因的表达水平进行相关性分析，
# 从而构建出probe-gene对。这个包具体的做法是分别找到差异甲基化远端探针（distal probe）
# 的上游和下游最近的10个基因，分析其余对应探针甲基化水平是否存在负相关，
# 从而筛选出来的潜在probe-gene对。
# 这一步ELMER包同样提供了两类方式：
# （1）supervised，即通过上游分析疾病组vs.正常组差异远端探针，
# 进一步纳入所有样本通过非参数检验的方式，找出同时在疾病组高表达的基因；
# （2）unsupervised，根据基因的distal probe甲基化水平进行排序，
# 提取出甲基化水平分别在前20%和后20%的样本，分别为M和U组，
# 通过U检验方式比较两组间基因表达水平,如果在M组的表达水平显著低于U组，
# 并通过迭代的方式计算出所有远端探针（distal probe）和靶基因的相关性P value,
# 从而筛选出具有意义的probe-gene关系对。

# 加载我们之前的数据
mae <- get(load("mae.rda"))
#加载我们筛选得到的差异甲基化远端探针（distal probe）
sig.diff <- read.csv("result/getMethdiff.hypo.probes.significant.csv")
#找到差异甲基化远端探针（distal probe）上下游10个基因
nearGenes <- GetNearGenes(data = mae, 
                          probes = sig.diff$probe, 
                          numFlankingGenes = 20) # 分别是上游10个和下游10基因
#因为前面使用unsupervised，我们这而也用unsupervised
Hypo.pair <- get.pair(data = mae,
                      group.col = "definition",
                      group1 =  "Primary solid Tumor",
                      group2 = "Solid Tissue Normal",
                      nearGenes = nearGenes,
                      mode = "unsupervised",
                      permu.dir = "result/permu",
                      permu.size = 100, # 一般要设置迭代次数为10000，处于演示目的减少次数
                      raw.pvalue = 0.05,   
                      Pe = 0.01, 
                      filter.probes = TRUE, 
                      filter.percentage = 0.05,
                      filter.portion = 0.3,
                      dir.out = "result",#结果同样输出到result目录
                      cores = 1,
                      label = "hypo")
#展示一下probe-gene对
head(Hypo.pair)


# step4--四、筛选后探针的motif富集分析 -------------------------------------------------------------------
# 这一步我们通过上一步分析得到的probe-gene对中的探针用以富集分析。
# 提取探针上下游250bp区域的碱基序列，从而找到富集的motif。
# 接下来再通过转录因子（TF）结合motif数据库，从而预测出motif结合的转录因子。
enriched.motif <- get.enriched.motif(data = mae,
                                     probes = Hypo.pair$Probe, 
                                     dir.out = "result", 
                                     label = "hypo",
                                     min.incidence = 10,
                                     lower.OR = 1.1)
# 这一步里，我们需要输入我们的MAE对象，
# 筛选出prob-gene对中的probes，并且还需要设置motif最低95%置信区间
# OR阈值以及最少几个probe进行富集motif。
#展示一下名称
names(enriched.motif)
#展示一下每个motif是通过哪些探针进行富集的
enriched.motif[[1]]
# # 除此之外，ELMER包还提供了两个表：
# “getMotif.hypo.enriched.motifs.rda” 
# “getMotif.hypo.motif.enrichment.csv”，大家可以在result目录看一下。



# step5--确定调控作用的转录因子 ------------------------------------------------------
# 这一步，ELMER包对motif以及上游转录因子的关系对进行筛选，
# 从而得到具有调控作用的转录因子。如果某一特定亚组中基因的enhancer发生改变，
# 其上有调控的转录因子同样也会发生改变。基于此
# ，ELMER包将上一步我们通过富集分析找到富集的motif以及对应motif的转录因子进行分析
# ，筛选出对应转录水平发生改变的转录因子。
# （1）unsupervised，即将存在相同motif的远端探针（distal probe）
# 根据甲基化水平分为高甲基化组（一般取前20%）和低甲基化组（后20%），
# 比较其两组间对应TF的表达值是否存在差异。
# （2）Supervised，即对疾病组vs.对照组相同motif对应TF表达值进行分析。
# 接下来筛选出前5% P value最小的TFs，并视为潜在的上游调控转录因子。

##找到调控motif对应的转录因子
TF <- get.TFs(data = mae, 
              group.col = "definition",
              group1 =  "Primary solid Tumor",
              group2 = "Solid Tissue Normal",
              mode = "unsupervised",
              enriched.motif = enriched.motif,
              dir.out = "result", 
              cores = 1, 
              label = "hypo")
#展示一下结果
head(TF)
dir(path = "result", pattern = "getTF")  
# TF ranking plot based on statistics will be automatically generated.
dir(path = "result/TFrankPlot/", pattern = "pdf")


# step6--可视化 --------------------------------------------------------------
mae <- get(load("mae.rda"))

scatter.plot(data = mae,
             byProbe = list(probe = c("cg19403323"), numFlankingGenes = 20), 
             category = "definition", 
             lm = TRUE, # Draw linear regression curve
             save = FALSE) 
# 每个散点图显示了样本探针cg19403323在所有LUSC样本中的甲基化水平与20个相邻基因中的一个表达的关系
scatter.plot(data = mae,
             byPair = list(probe = c("cg19403323"), gene = c("ENSG00000143469")), 
             category = "definition", save = TRUE, lm_line = TRUE) 
# Scatter plot shows the methylation level of an example probe cg19403323 in all LUSC samples plotted against the expression of the putative target gene SYT14.
######## 展示转录因子和甲基化表达水平的散点图
load("result/getMotif.hypo.enriched.motifs.rda")
names(enriched.motif)[1]
scatter.plot(data = mae,
             byTF = list(TF = c("TP53","SOX2"),
                         probe = enriched.motif[[names(enriched.motif)[1]]]), 
             category = "definition",
             save = TRUE, 
             lm_line = TRUE)
# Each scatter plot shows the average methylation level of sites with the first enriched motif in all LUSC samples plotted against the expression of the transcription factor TP53, SOX2 respectively.

# step7--可视化 示意图 ----------------------------------------------------------
# Load results from previous sections
mae <- get(load("mae.rda"))
pair <- read.csv("result/getPair.hypo.pairs.significant.csv")
# 基因组图
# 展示邻近基因和探针的关系
# Generate schematic plot for one probe with 20 nearby genes and label the gene significantly linked with the probe in red.
schematic.plot(pair = pair, 
               data = mae,
               group.col = "definition",
               byProbe = pair$Probe[1],
               save = FALSE)
# 生成一个探针与20个邻近基因的示意图
# ，用红色标记与探针显著相关的基因。

# 生成和基因显著关联的附近的探针的图 
schematic.plot(pair = pair, 
               data = mae,   
               group.col = "definition", 
               byGene = pair$GeneID[1],
               save = FALSE)
# 图中显示了红色的基因和所有蓝色的探针，它们与该基因的表达显著相关。


# step8-- motif富集图 --------------------------------------------------------
# motif 转录因子结合基序
# 通俗来说就是转录因子结合的一个序列
# 类似于森林图，还有OR值
motif.enrichment.plot(motif.enrichment = "result/getMotif.hypo.motif.enrichment.csv", 
                      significant = list(OR = 1.5,lowerOR = 1.3), 
                      label = "hypo", 
                      save = FALSE)  

motif.enrichment.plot(motif.enrichment = "result/getMotif.hypo.motif.enrichment.csv", 
                      significant = list(OR = 1.5,lowerOR = 1.3), 
                      label = "hypo", 
                      summary = TRUE,
                      save = FALSE)

####### 调节TF图########
# 所有TF按照-logp值排序，
# p值是按照和对应motif甲基化水平反相关算来的的
# 最不相关的TF排在最前面
load("result/getTF.hypo.TFs.with.motif.pvalue.rda")
motif <- colnames(TF.meth.cor)[1]
TF.rank.plot(motif.pvalue = TF.meth.cor, 
             motif = motif,
             save = FALSE) 
# 对于给定的motif，所有转录因子进行排名
# 最相关的前三个标记出来


# step9--可视化：热图 -----------------------------------------------------------
mae <- get(load("mae.rda"))
# Generate a heatmap of paired probes. 表达矩阵和DNA甲基化的热图
pair <- read.csv("result/getPair.hypo.pairs.significant.csv")

heatmapPairs(data = mae, 
             group.col = "definition",
             group1 = "Primary solid Tumor", 
             annotation.col = c("years_smoked","gender"),
             group2 = "Solid Tissue Normal",
             pairs = pair,
             filename =  NULL)

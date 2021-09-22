# 构建好表达矩阵就做生存分析曲线。
rm(list = ls())
p_load(ggplot2, ggpubr, dplyr, survival, survminer, stringr,export)

###### 准备数据，包括读入fpkm_finish 和 survival 和 phenotype 以及加载好包
cln_xena<-read.csv("data/LUAD_survival.csv") # 738个sample
allid <- read.csv(file = "outdata/allid_493.txt", row.names = 1) %>% .[,1]
# 
# fpkm_id_finish <- read.csv("outdata/mrna_fpkm_id_finish.txt",header = T,
#                            row.names = 1,
#                            check.names = F,
#                            stringsAsFactors = F)
# 从cln中选择fpkm中存在的样品ID做生存分析
cln_xena1 <- cln_xena[cln_xena$sample %in% allid, ]
# 查看有没有na值,有就去掉na值
table(is.na(cln_xena1))
cln_xena1 <- na.omit(cln_xena1)


########   对clinical数据 进行生存分析数据的准备   #######


# 生存时间除以30得到月份数
cln_xena1$OS.time <- cln_xena$OS.time[cln_xena$sample %in% allid] / 30
# cln_xena1$OS.time <- cln_xena$OS.time[as.numeric(rownames(cln_xena1))] / 30
# cln_xena1$OS.time <- cln_xena$OS.time[cln_xena$sample %in% cln_xena1$sample] / 30 方法2

# cln加一列mut代表分组信息 
cln_xena1$group <- 'wild'
cln_xena1$group[ cln_xena1$sample %in% allid[1:111]  ] <- 'mut'
# cln_xena1$group[match(x = colnames(fpkm_finish)[1:111], table = cln_xena1$sample,nomatch = 0)] <- 'mut'

# 检查group分组信息情况
table(cln_xena1$group)
# group必须是因子，代表有无突变
str(cln_xena1)
cln_xena1$group <- as.factor(cln_xena1$group)
write.csv(cln_xena1, file = "outdata/cln_xena.csv", row.names = F)
################## 读入整理好的生存数据，做KM图###############
rm(list = ls())
cln_xena1 <- read.csv("outdata/cln_xena.csv")
# 以下代码直接运行即可
data1 <- cln_xena1 %>% filter(OS.time<=62)

fit<-survfit(Surv(OS.time,OS)~group, data = data1) 

surv_summary(fit) #查看生存率及其标准误
surv_median(fit = fit,combine = F) # 查看中位生存时间
pp <- ggsurvplot(fit, data=data1, linetype = 1, palette = c("jco"),size=1,
                 surv.scale = c("percent"),pval = TRUE,
                 legend.title = '', legend.labs=c("Mut","Wild"),
                 break.time.by =12,
                 xlim = c(0,60),risk.table = TRUE,
                 risk.table.title = "Patients at risk",
                 ylab = "Overall Survival, %",xlab = "Months",
                 font.x = c(20,"plain","black"),font.y = c(20,"plain","black"),
                 font.tickslab = c(20,"plain","black"),
                 risk.table.fontsize = 6,font.legend =c(20,"plain","black"),
                 font.main = c(20,"plain","black"),pval.size = 8)

######risk table 的修改
pp$table <- ggpar(
  pp$table,
 font.title    = c(13, "bold.italic", "green"),  ###risk table 标题的修改
font.subtitle = c(15, "bold", "pink"),  ###risk table小标题的修改
font.caption  = c(11, "plain", "darkgreen"), ####插入字的修改
  font.x        = c(18, "plain", "black"), ### risk table x的修改
  font.y        = c(18, "plain", "black"),### risk table y的修改
  font.xtickslab = c(18, "plain", "black"),### risk table x 坐标轴的修改
  font.ytickslab = c(18),  ### risk table y 坐标轴的修改
  legend=c(0.8,0.88),
  censor.size=3
)
print(pp)
graph2ppt(file="output/plots/kmplot_os_by_xena_survival.pptx")


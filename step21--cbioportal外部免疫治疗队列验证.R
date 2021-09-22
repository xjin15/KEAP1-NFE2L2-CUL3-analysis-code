# 免疫治疗队列 ------------------------------------------------------------------
# 读取数据
p_load(ggplot2, ggpubr, dplyr, survival, survminer, stringr,export)

# 外部验证之生存分析 ---------------------------------------------------------------
rm(list = ls())
cbioportal <- data.table::fread("data/other_dataset/survival——data.csv")

cbioportal |> str()
table(cbioportal$Status)
cbioportal <- cbioportal |> dplyr::mutate(
  group = as.factor(group),
  OS = ifelse(Status == "censored",0,1),
)
colnames(cbioportal)[6] <- "OS.time"
cln_xena1 <- cbioportal[,c(1,6:8)]
# 以下代码直接运行即可
data1 <- cln_xena1 #%>% filter(OS.time<=96)

fit<-survfit(Surv(OS.time,OS)~group, data = data1) 

surv_summary(fit) #查看生存率及其标准误
surv_median(fit = fit,combine = F) # 查看中位生存时间
pp <- ggsurvplot(fit, data=data1, linetype = 1, palette = c("jco"),size=1,
                 surv.scale = c("percent"),pval = TRUE,
                 legend.title = '', legend.labs=c("Mut","Wild"),
                 break.time.by =12,
                 xlim = c(0,72),risk.table = TRUE,
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
graph2ppt(file="output/plots/外部验证.pptx",append = T,)

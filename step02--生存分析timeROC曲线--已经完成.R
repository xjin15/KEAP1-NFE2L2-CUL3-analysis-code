
# time-ROC曲线 --------------------------------------------------------------------
p_load(timeROC,survival,survivalROC)
data(mayo)
rocdata_os <- dt2 %>% 
  select(sample, OS,OS.time,group,ajcc_T,ajcc_N,radiotherapy)

time_roc_res <- timeROC(
  T = rocdata_os$OS.time,
  delta = rocdata_os$OS,
  marker = rocdata_os$group,
  other_markers=as.matrix(rocdata_os[,c("ajcc_T","ajcc_N","radiotherapy")]),
  cause = 1,
  weighting="cox",
  times = c(1 * 12, 3 * 12, 5 * 12),
  ROC = TRUE,
  iid = TRUE
)

Mayo4.2 <-  survivalROC(Stime=rocdata_os$OS.time,
                     status=rocdata_os$OS,
                     marker = rocdata_os$group,
                     predict.time = 12,
                     method = "KM")
# time_roc_res <- timeROC(
#   T = mayo$time,
#   delta = mayo$censor,
#   marker = mayo$mayoscore5,
#   cause = 1,
#   weighting="marginal",
#   times = c(3 * 365, 5 * 365, 10 * 365),
#   ROC = TRUE,
#   iid = TRUE
# )
time_roc_res$AUC   # 查看AUC值
confint(time_roc_res, level = 0.95)$CI_AUC # 查看95%置信区间
# 这里的plot函数对应的是timeROC::plot.ipcwsurvivalROC函数）
plot(time_roc_res, time=1 * 12, col = "red", title = FALSE)  
plot(time_roc_res, time=3 * 12, add=TRUE, col="blue") 
plot(time_roc_res, time=5 * 12, add=TRUE, col="green") 
legend("bottomright",c("1 Years" ,"3 Years", "5 Years"),
       col=c("red", "blue", "green"), lty=1, lwd=2)

# 修改和美观time-ROC曲线图 --------------------------------------------------------

time_ROC_df <- data.frame(
  TP_3year = time_roc_res$TP[, 1],
  FP_3year = time_roc_res$FP[, 1],
  TP_5year = time_roc_res$TP[, 2],
  FP_5year = time_roc_res$FP[, 2],
  TP_10year = time_roc_res$TP[, 3],
  FP_10year = time_roc_res$FP[, 3]
)
library(ggplot2)
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_10year, y = TP_10year), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 10 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )






# 比较两个ROC曲线图 --------------------------------------------------------------
# 对mayoscore4做类似的分析

time_roc_res2 <- timeROC(
  T = mayo$time,
  delta = mayo$censor,
  marker = mayo$mayoscore4,
  cause = 1,
  weighting="marginal",
  times = c(3 * 365, 5 * 365, 10 * 365),
  ROC = TRUE,
  iid = TRUE
)
time_roc_res2$AUC

# 然后通过compare函数进行比较，
# 并输出矫正后的P值和相关系数矩阵
# 假设检验的原假设是两个AUC是相等的
compare(time_roc_res, time_roc_res2, adjusted = TRUE)

# 接着可通过plotAUCcurve函数绘制不同时间节点的AUC曲线及其置信区间
# 也可将多个ROC曲线的AUC值放在一起绘制（节点多一点，曲线会展示的更加细致一点）
plotAUCcurve(time_roc_res, conf.int=TRUE, col="red")
plotAUCcurve(time_roc_res2, conf.int=TRUE, col="blue", add=TRUE)
legend("bottomright",c("mayoscore5", "mayoscore4"), col = c("red","blue"), lty=1, lwd=2)


# ROC的最佳阈值 ----------------------------------------------------------------
# 对于上述timeROC的结果，
# 如3年ROC曲线的约登指数
# （因为TP代表的是True Positive fraction，即sensitivity；
# 而FP代表的是False Positive fraction，即1-specificity）：
mayo$mayoscore5[which.max(time_ROC_df$TP_3year - time_ROC_df$FP_3year)]
# 即对于mayoscore5这个marker而言，最佳阈值（cutoff）为6.27




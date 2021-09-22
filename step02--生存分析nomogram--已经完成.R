############### cox 分析和nomogram############

########################cox回归分析####################
#############1单因素cox分析

#########phe_f是完全整理好的临床数据！包括各种因素的数据情况
p_load(ggplot2, ggpubr, dplyr, survival, survminer, stringr,export)

phe_f <- phenotype_xena3 %>% 
  select(c("sample", "OS", "OS.time","PFS", "PFS.time", "group", 
                            "gender", "age", "age_median","age_65","race",  "smoke",
                            "smoke_group", "ajcc_T", "ajcc_N", "ajcc_M", "stage_group", "stage2",
                            "resection_site","radiotherapy"))

# 以下代码改写成within形式，注意必须用within，不能用with，也不能用attach，因为不是全局变量
# phe_f$age_group2 <- "<60"
# phe_f$age_group2[phe_f$age >= 60 &
#                  phe_f$age < 70  ] <- "60-70"
# phe_f$age_group2[phe_f$age > 70  ] <- ">70"
# phe_f <- within(data = phe_f, {age_group2 <- NA
#                               age_group2[age > 70] <- ">70"
#                               age_group2[age <= 70 & age >= 60] <- "60-70"
#                               age_group2[age < 60] <- "<60"  })
# 还有一种方法
# phe_f <- phe_f %>% 
#   mutate(age_group2= ifelse(age > 70, ">70", 
#                             ifelse(age < 60, "<60", 
#                                    ifelse(is.na(age), NA, "60-70") )))
str(phe_f)
# 把分类变量都转成因子

phe_f <- phe_f %>% 
  mutate(gender = factor(gender, levels = c("female","male")),
         group = factor(group, levels = c("wild","mut")),
         age_median = factor(age_median, levels = c("younger","older")),
         age_65 = factor(age_65, levels = c("<65","≥65")),
         race = factor(race, levels = c("other", "black", "white" ), ordered = F),
         smoke = factor(smoke, levels = c(1,3,4,5,2)),
         smoke_group = factor(smoke_group, levels = c("NO", "YES"), ordered = F),
         ajcc_T = as.factor(ajcc_T),
         ajcc_N = as.factor(ajcc_N),
         ajcc_M = factor(ajcc_M, levels = c("M0","M1","MX")),
         stage_group = as.factor(stage_group),
         stage2 = as.factor(stage2),
         resection_site = factor(resection_site, levels = c("Lower lobe","Middle lobe","Upper lobe", "other site")),
         age_group2= ifelse(age > 70, ">70", 
                            ifelse(age < 60, "<60", 
                                   ifelse(is.na(age), NA, "60-70") )),
         age_group2 = factor(age_group2, levels = c("<60", "60-70", ">70")),
         radiotherapy = factor(radiotherapy, levels = c("NO", "YES")),
         PFS = as.integer(PFS),
         OS = as.integer(OS)  ) %>% 
  select(sample, OS, OS.time, PFS, PFS.time, group, gender, age, age_median, age_65, age_group2, everything())
str(phe_f)
write.csv(phe_f, file = "outdata/phe_finish.csv")
######## data1 筛选PFS时间小于等于5年的进行分析
library(dplyr)
phe_f <- read.csv(file = "outdata/phe_finish.csv", row.names = 1)
data1 <- phe_f %>% filter(PFS.time <= 62)
ppcox <- coxph(Surv(PFS.time, PFS) ~ group, data = data1)
summary(ppcox)   ##不去除na值的p=0.016
############ data1 去除na值
data1 <- phe_f %>% filter(PFS.time <= 62)
data1 <- na.omit(data1)
ppcox <- coxph(Surv(PFS.time, PFS) ~ group, data = data1)
summary(ppcox)  ###去除na值得p=0.157
# 还是不去除

####################### 对年龄的单因素cox分析：连续OR分类？#############
{
  
  
  ### 当作连续变量
  agecox1 <- coxph(Surv(PFS.time,PFS==1)~age,data1)
  summary(agecox1)   # P==0.459
  
  # 年龄根据<60, 60-70, >70 3分类
  agecox2 <- coxph(Surv(PFS.time,PFS==1)~age_group2,data1)
  summary(agecox2)# likelihood p=0.5
  
  ### 年龄用ROC曲线，选择最佳截断值二分类。
  library(pROC) 
  #构建关于结局和年龄的函数
  aa<-roc(data1$PFS,data1$age) 
  ##roc曲线
  plot(aa,
       print.thres=TRUE,
       main="ROC曲线",
       col="#008600")
  ### 图上显示，最佳截断点是 65岁
  agecox3 <- coxph(Surv(PFS.time,PFS==1)~age_median, data1)
  summary(agecox3) # 0,195
  
  agecox4 <- coxph(Surv(time = OS.time, event = OS==1)~age_group2, data = data1 )
  summary(agecox4) # P=0.07
  
  ppcoxa <- coxph(Surv(PFS.time, PFS) ~ stage2, data = data1)
  a <- summary(ppcoxa)
  a
  ppcoxb <- coxph(Surv(PFS.time, PFS) ~ stage_group, data = data1)
  b <- summary(ppcoxb)
  b
}
############## 批量单因素cox分析#########################
phe_f <- phenotype_xena3[ c("sample", "OS", "OS.time","PFS", "PFS.time", "group", 
                            "gender", "age", "age_median","age_65","race",  "smoke",
                            "smoke_group", "ajcc_T", "ajcc_N", "ajcc_M", "stage_group", "stage2",
                            "resection_site","radiotherapy")]
library(mice)
md.pattern(phe_f)# 查看缺失值的分布
library(VIM)
aggr(phe_f)
complete.cases(phe_f) %>% table()
phe_f <- within(phe_f, {
  radiotherapy[is.na(radiotherapy)] <- "YES" #NO/Unknow
  age[is.na(age)] <- median(age,na.rm = T)
  smoke[is.na(smoke)] <- 7
  smoke_group[is.na(smoke_group)] <- "NO"
  ajcc_N[ajcc_N=="N2"|ajcc_N=="N3"] <- "N2+N3"
  ajcc_N[ajcc_N== "NX"] <- NA
})
complete.cases(phe_f) %>% table()


phe_f <- na.omit(phe_f)

phe_f <- phe_f %>% 
  mutate(gender = factor(gender, levels = c("female","male")),
         group = factor(group, levels = c("wild","mut")),
         age_median = factor(age_median, levels = c("younger","older")),
         age_65 = factor(age_65, levels = c("<65","≥65")),
         race = factor(race, levels = c("other", "black", "white" ), ordered = F),
         smoke = as.factor(smoke),
         smoke_group = as.factor(smoke_group),
         ajcc_T = as.factor(ajcc_T),
         ajcc_N = as.factor(ajcc_N),
         ajcc_M = factor(ajcc_M, levels = c("M0","M1","MX")),
         stage_group = as.factor(stage_group),
         stage2 = as.factor(stage2),
         resection_site = factor(resection_site, levels = c("Lower lobe","Middle lobe","Upper lobe", "other site")),
         age_group2= ifelse(age > 70, ">70", 
                            ifelse(age < 60, "<60", 
                                   ifelse(is.na(age), NA, "60-70") )),
         age_group2 = factor(age_group2, levels = c("<60", "60-70", ">70")),
         radiotherapy = as.factor(radiotherapy),#ifelse(radiotherapy=="YES", "YES","NO/Unknow"),
         PFS = as.integer(PFS),
         OS = as.integer(OS)  ) %>% 
  select(sample, OS, OS.time, PFS, PFS.time, group, gender, age, age_median, age_65, age_group2, everything()) # %>% 
  # mutate(radiotherapy = as.factor(radiotherapy))
str(phe_f)
data1 <- phe_f %>% filter(PFS.time <= 62)
data2 <- phe_f %>% filter(OS.time <= 62)
#1.构建函数
pfs_surv<- Surv(time = data1$PFS.time, event = data1$PFS==1)
data1$surv <- with(data1,pfs_surv)  
# 相当于 data1$surv1 <- pfs_surv 这两个是等价的
cox <- coxph(surv~group, data=data1)
summary(cox)
Unicox<- function(x){ 
  FML <- as.formula(paste0 ("surv~",x))
  cox <- coxph(FML,data=unidata)
  cox1 <- summary(cox)
  # beta <- round(cox1$coef[,1],3) # betaβ值就是coef偏回归系数
  HR <- round(cox1$coefficients[,2],2)  # HR 就是exp(coef)  
  PValue <- round(cox1$coefficients[,5],3)  # p值就是Pr(>|z|)
  xname <- rownames(cox1$conf.int)
  CI <- paste(round(cox1$conf.int[,3],2),round(cox1$conf.int[,4],2), sep = '-')  # 95%可信区间就是lower .95和 upper .95
  Unicox<- data.frame('Characteristics' = xname,
                      'Hazard Ratio' = HR,
                      'P Value' = PValue,
                      'CI95' = CI)
  return(Unicox)}      

#2.提取变量，构建数据框，批量生成表格
variable.names<- colnames(data1)[which(names(data1)=="group"):(ncol(data1)-1)] 
unidata <- data1
pfs_Univar <- lapply(variable.names, Unicox)

pfs_Univar <- plyr::ldply(pfs_Univar,data.frame)
pfs_Univar #查看表格结果
#3 一些修饰
pfs_Univar$P.Value[which(pfs_Univar$P.Value == 0 )] <- "<0.001"
pfs_Univar
# 挑选出pfs需要的因素group, age，gender，ajcc_T, ajcc_N, radiotherapy。
write.csv(pfs_Univar,"outdata/PFS单因素回归分析表.csv")
##用data2 做一遍 OS的分析
os_surv<- Surv(time = data2$OS.time, event = data2$OS==1)
data2$surv <- with(data2,os_surv)  
cox <- coxph(surv~group, data=data2)
summary(cox)

variable.names<- colnames(data2)[which(names(data2)=="group"):(ncol(data2)-1)] 
unidata <- data2
os_Univar <- lapply(variable.names, Unicox)
library(plyr)
os_Univar <- ldply(os_Univar,data.frame)
os_Univar$P.Value[which(os_Univar$P.Value == 0 )] <- "<0.001"
os_Univar 
write.csv(os_Univar,"outdata/OS单因素回归分析表.csv")

# 挑选出os需要的因素group, age，gender，ajcc_T, ajcc_N, radiotherapy, 。

#############多因素cox分析################


##########对PFS的cox分析

pfs_multicox <- coxph(Surv(PFS.time, PFS) ~ group+age+smoke_group+gender+ajcc_T+ajcc_N+radiotherapy, data = data1)
step(pfs_multicox)
pfs_multicox <- coxph(Surv(PFS.time, PFS) ~ group+age+ajcc_T+ajcc_N+radiotherapy, data = data1)
summary(pfs_multicox)  # p=0.01386
########提取cox多因素回归分析三线表##################
library(tableone)
library(broom)
#2 提取HR.P.95%CI
pfs_multi1<-ShowRegTable(pfs_multicox, 
                         exp=TRUE, 
                         digits=2, 
                         pDigits =3,
                         printToggle = TRUE, 
                         quote=FALSE, 
                         ciFun=confint)
#3 提取回归系数、统计量等                     
pfs_multi2<-tidy(pfs_multicox) #broom包
pfs_multi2
#4 将两次提取结果合并
pfs_multi<-cbind(pfs_multi1, pfs_multi2)
pfs_multi
write.csv(pfs_multi,file="outdata/pfs_multi.csv")

###森林图
ggforest(model = pfs_multicox, data = data1,main ="Forest of PFS multicox analysis")
graph2ppt(file="output/plots/forestplot_PFS")
graph2pdf(file="output/plots/pdf/forestplot_PFS")

###########对OS的cox分析

os_multicox <- coxph(Surv(OS.time, OS) ~ group+age+gender+ajcc_T+ajcc_N+radiotherapy, data = data2)
step(os_multicox)
os_multicox <- coxph(Surv(OS.time, OS) ~ group+ajcc_T+ajcc_N+radiotherapy, data = data2)
summary(os_multicox)  # p=0.007
ggforest(model = os_multicox, data = data2, main ="Forest of OS multicox Analysis")
graph2ppt(file="output/plots/forestplot_OS")
graph2pdf(file="output/plots/pdf/forestplot_OS")

#2 提取HR.P.95%CI
os_multi1<-ShowRegTable(os_multicox, 
                        exp=TRUE, 
                        digits=2, 
                        pDigits =3,
                        printToggle = TRUE, 
                        quote=FALSE, 
                        ciFun=confint)
#3 提取回归系数、统计量等                     
os_multi2<-tidy(os_multicox)
os_multi2
#4 将两次提取结果合并
os_multi<-cbind(os_multi1, os_multi2)
os_multi
write.csv(os_multi,file="outdata/os_multi.csv")




################ 用data1  PFS做nomogram图 使用了stage_group而不是TN分期#################
#需要加一个 R包
dt1 <- data1[-ncol(data1)]
dt1 <- na.omit(dt1)
uu <- dt1
y<-Surv(uu$PFS.time,uu$PFS==1)
kmfit<-survfit(y~1, data=y)
coxmodel <- coxph(y ~ group+age+gender+smoke_group+ajcc_T+ajcc_N+radiotherapy, data=uu)
step(coxmodel) ### 这一步是挑选合适的因素进行cox多因素分析，但其实我在这之前已经挑选过了
coxmodel <- coxph(y ~ group+age+ajcc_T+ajcc_N+radiotherapy, data=uu)
summary(coxmodel)
library(rms)
dd<-datadist(uu)
options(datadist="dd")

f <- cph( y ~ group+age+ajcc_T+ajcc_N+radiotherapy, x=T, y=T,surv=T, data=uu)
surv<-Survival(f)
nom1 <- nomogram(f, fun=list(function(x) surv(12, x), function(x) surv(36, x), function(x) surv(60, x)), lp=F,funlabel=c("1-year survival ", "3-year survival ", "5-year survival "), maxscale=10, fun.at=c(0.95, 0.9, 0.8, 0.7, 0.6, 0.5,0.4,0.3,0.2,0.1,0.05))

plot(nom1,cex.axis = 1.3,cex.var =1.5,xfrac=.4,lmgp=0.5,tcl=-0.5)
graph2ppt(file="output/plots/nomogram_pfs_舍弃", width=10,height=10)
graph2pdf(file="output/plots/nomogram_pfs_舍弃", width=10,height=10)

# 内部验证 
validate(f, method="boot", B=1000, dxy=T)
rcorrcens(Surv(PFS.time, PFS) ~ predict(f), data = uu) #C—index = 1-C==1-0.355=0.645

# 一致性检验
# 1年
f1 <- cph( y ~ group+age+ajcc_T+ajcc_N+radiotherapy, x=T, y=T,surv=T, data=uu,time.inc=12)
cal1 <- calibrate(f1, cmethod="KM", method="boot", u=12, m=100, B=1000)
{
  plot(cal1, lwd=2,cex.axis=1.5,tcl=-0.8, lty=1, errbar.col=c(rgb(0,118,192,maxColorValue=255)), xlab=" ",ylab = "",col=c(rgb(192,98,83,maxColorValue=255)))
  lines(cal1[,c("mean.predicted","KM")], type="b", lwd=3, col=c(rgb(192,98,83,maxColorValue=255)), pch=16)
  mtext("Nomogram predicted 1-year PFS ", side = 1, line = 2.5,cex=1.7)
  mtext("Actual 1-year PFS (proportion)  ", side = 2, line = 2.5,cex=1.7)
  box(lwd = 1)
  abline(0, 1, lty=3, lwd=2, col=c(rgb(0,118,192,maxColorValue=255)))
}
graph2ppt(file="output/plots/nomogram_pfs_cal1year")
# 3 年
f2<-cph( y ~ group+age+gender+stage_group+radiotherapy, x=T, y=T,surv=T, data=uu,time.inc=36)
cal2<-calibrate(f2, cmethod="KM", method="boot", u=36, m=80, B=1000)
{
  plot(cal2, lwd=2, cex.axis=1.5,tcl=-0.8,lty=1, errbar.col=c(rgb(0,118,192,maxColorValue=255)), xlab=" ",ylab = "",col=c(rgb(192,98,83,maxColorValue=255)))
  lines(cal2[,c("mean.predicted","KM")], type="b", lwd=3, col=c(rgb(192,98,83,maxColorValue=255)), pch=16)
  mtext("Nomogram predicted 3-year PFS ", side = 1, line = 2.5,cex=1.7)
  mtext("Actual 3-year PFS (proportion)  ", side = 2, line = 2.5,cex=1.7)
  box(lwd = 1)
  abline(0, 1, lty=3, lwd=2, col=c(rgb(0,118,192,maxColorValue=255)))
}
graph2ppt(file="output/plots/nomogram_pfs_cal3year")

# 5年  5年生存率做不了，没有这么多数据
# f3 <- cph( y ~ group+age+gender+stage_group+radiotherapy, x=T, y=T,surv=T, data=uu,time.inc=60)
# cal3 <- calibrate(f3, cmethod="KM", method="boot", u=60, m=100, B=1000)
# {
# plot(cal3, lwd=2,cex.axis=1.5,tcl=-0.8, lty=1, errbar.col=c(rgb(0,118,192,maxColorValue=255)), xlab=" ",ylab = "",col=c(rgb(192,98,83,maxColorValue=255)))
# lines(cal3[,c("mean.predicted","KM")], type="b", lwd=3, col=c(rgb(192,98,83,maxColorValue=255)), pch=16)
# mtext("Nomogram predicted 5-year PFS ", side = 1, line = 2.5,cex=1.7)
# mtext("Actual 5-year PFS (proportion)  ", side = 2, line = 2.5,cex=1.7)
# box(lwd = 1)
# abline(0, 1, lty=3, lwd=2, col=c(rgb(0,118,192,maxColorValue=255)))
# }


################ 用data2做nomogram图###########
#需要加一个 R包
dt2 <- data2[-ncol(data2)]
dt2 <- na.omit(dt2)
uu <- dt2
y<-Surv(uu$OS.time,uu$OS==1)
kmfit<-survfit(y~1, data=y)
coxmodel <- coxph(y ~ group+age+gender+ajcc_T+ajcc_N+smoke_group+radiotherapy, data=uu)
step(coxmodel) ### 这一步是挑选合适的因素进行cox多因素分析，但其实我在这之前已经挑选过了
coxmodel <- coxph(y ~ group+ajcc_T+ajcc_N+radiotherapy, data=uu)
library(rms)
dd<-datadist(uu)
options(datadist="dd")

f <- cph( y ~ group + ajcc_T + ajcc_N + radiotherapy, x=T, y=T,surv=T, data=uu)
surv<-Survival(f)
nom2 <- nomogram(f, fun=list(function(x) surv(12, x), function(x) surv(36, x), function(x) surv(60, x)), lp=F,funlabel=c("1-year survival ", "3-year survival ", "5-year survival "), maxscale=10, fun.at=c(0.95, 0.9, 0.8, 0.7, 0.6, 0.5,0.4,0.3,0.2,0.1,0.05))
plot(nom2,cex.axis = 1.3,cex.var =1.5,xfrac=.4,lmgp=0.5,tcl=-0.5)
graph2ppt(file="output/plots/nomogram_os_simple_TN", width=10,height=10)

# 内部验证 C—index = 1-C==1-0.354=0.646
validate(f, method="boot", B=1000, dxy=T)
rcorrcens(Surv(OS.time, OS) ~ predict(f), data = uu) #C—index = 1-C==1-0.328=0.672
Cindex_os_nomo <- summary(coxmodel) 
Cindex_os_nomo <- Cindex_os_nomo$concordance[1] ### C-index 也可以这样计算

# 一致性检验
# 1年
f1 <- cph( y ~ group+ajcc_T+ajcc_N+radiotherapy, x=T, y=T,surv=T, data=uu,time.inc=12)
cal1 <- calibrate(f1, cmethod="KM", method="boot", u=12, m=100, B=1000)
##### 调整u,m,B 的值可以获得好看的calibration图 见幕布学习笔记
{
  plot(cal1, lwd=2,cex.axis=1.5,tcl=-0.8, lty=1, errbar.col=c(rgb(0,118,192,maxColorValue=255)), xlab=" ",ylab = "",col=c(rgb(192,98,83,maxColorValue=255)))
  lines(cal1[,c("mean.predicted","KM")], type="b", lwd=3, col=c(rgb(192,98,83,maxColorValue=255)), pch=16)
  mtext("Nomogram predicted 1-year OS ", side = 1, line = 2.5,cex=1.7)
  mtext("Actual 1-year OS (proportion)  ", side = 2, line = 2.5,cex=1.7)
  box(lwd = 1)
  abline(0, 1, lty=3, lwd=2, col=c(rgb(0,118,192,maxColorValue=255)))
}
graph2ppt(file="output/plots/nomogram_os_cal1year")
graph2pdf(file="output/plots/nomogram_os_cal1year")
# 3 年
f2<-cph( y ~ group+ajcc_T+ajcc_N+radiotherapy, x=T, y=T,surv=T, data=uu,time.inc=36)
cal2<-calibrate(f2, cmethod="KM", method="boot", u=36, m=100, B=1000) 
{
  plot(cal2, lwd=2, cex.axis=1.5,tcl=-0.8,lty=1, errbar.col=c(rgb(0,118,192,maxColorValue=255)), xlab=" ",ylab = "",col=c(rgb(192,98,83,maxColorValue=255)))
  lines(cal2[,c("mean.predicted","KM")], type="b", lwd=3, col=c(rgb(192,98,83,maxColorValue=255)), pch=16)
  mtext("Nomogram predicted 3-year OS ", side = 1, line = 2.5,cex=1.7)
  mtext("Actual 3-year OS (proportion)  ", side = 2, line = 2.5,cex=1.7)
  box(lwd = 1)
  abline(0, 1, lty=3, lwd=2, col=c(rgb(0,118,192,maxColorValue=255)))
}
graph2ppt(file="output/plots/nomogram_os_cal3year")
graph2pdf(file="output/plots/nomogram_os_cal3year")

# 5年  5年生存率做不了，没有这么多数据
# f3 <- cph( y ~ group+age+gender+stage_group+radiotherapy, x=T, y=T,surv=T, data=uu,time.inc=60)
# cal3 <- calibrate(f3, cmethod="KM", method="boot", u=60, m=150, B=1000)
# {
#   plot(cal3, lwd=2,cex.axis=1.5,tcl=-0.8, lty=1, errbar.col=c(rgb(0,118,192,maxColorValue=255)), xlab=" ",ylab = "",col=c(rgb(192,98,83,maxColorValue=255)))
#   lines(cal3[,c("mean.predicted","KM")], type="b", lwd=3, col=c(rgb(192,98,83,maxColorValue=255)), pch=16)
#   mtext("Nomogram predicted 5-year overall survival ", side = 1, line = 2.5,cex=1.7)
#   mtext("Actual 5-year survival (proportion)  ", side = 2, line = 2.5,cex=1.7)
#   box(lwd = 1)
#   abline(0, 1, lty=3, lwd=2, col=c(rgb(0,118,192,maxColorValue=255)))
#   }
# 


beta <- coef(ppcox)
se <- sqrt(diag(vcov(ppcox)))
HR <- exp(beta)
HRse <- HR * se
summary(ppcox)

################## 做nomogram图#################
library(regplot)
library(survival)
library(survminer)
library(regplot)

obs <- data2[5,-22]
obs
nomo <- regplot(os_multicox, observation=obs, plots = c("density","boxes"),failtime =c(12,36,60), prfail = TRUE,
                boxcol="#ADD8E6",cexvars=1.2,cexscales=1.2,cexcats=1.0,droplines=TRUE, points = T,
                title = "Nomogram to predict OS in LUAD patients")
#res.cox表示模型,可以是广义线性模型（glm）,线性（lm）,生存分析（cox比例风险模型）
#observation指定某个患者各协变量的取值映射到相应的得分，并计算总得分
#failtime = c(12,36,60)计算其在1\3\5年的的累计事件发生概率
#本案例中lung是个生存类数据，status=1代表构建的是生存模型，因此prfail = T

graph2ppt(file="output/plots/nomogram_regplot",width=12, aspectr = 1.5)


nomo




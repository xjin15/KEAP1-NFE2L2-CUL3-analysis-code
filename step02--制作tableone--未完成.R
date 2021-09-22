#############R包tableone制作文章的tableone###############
rm(list = ls())
p_load(tableone, survival, broom,tidyverse,compareGroups)
phe_f <- read.csv("outdata/phe_finish.csv",row.names = 1,stringsAsFactors = F)
head(phe_f)
str(phe_f)



# 数据整理成因子 -----------------------------------------------------------------


phe_f$radiotherapy <-   phe_f$radiotherapy %>% replace_na("No/Unknown") %>% recode( NO = "No/Unknown")

####age 有10个NA值，不过在连续变量中不需要管  race 没有缺失值 stage_group  有6个NA
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
         age_group2 = factor(age_group2, levels = c("<60", "60-70", ">70")),
         radiotherapy = factor(radiotherapy, levels = c("No/Unknown", "YES")),
         PFS = as.integer(PFS),
         OS = as.integer(OS)  ) %>% 
  select(sample, OS, OS.time, PFS, PFS.time, group, gender, age, age_median, age_65, age_group2, everything())
str(phe_f)
# 连续变量正态性检验
shapiro.test(phe_f$age) # 不符合正态分布
shapiro.test(phe_f$PFS.time) # 不符合正态分布
shapiro.test(phe_f$OS.time) # 不符合正态分布
dput(names(phe_f))



# tableone的使用 -------------------------------------------------------------


myVars <- c( "gender","age","age_65", "race", "smoke_group", "ajcc_T", "ajcc_N", "ajcc_M" ,"stage_group",  
            "resection_site", "radiotherapy")
catVars <- c( "gender", "race", "age_65", "smoke_group", "ajcc_T", "ajcc_N", "ajcc_M", "stage_group",  
              "resection_site", "radiotherapy") 
nonvar <- c("age")
exactvars <- c("ajcc_T", "ajcc_N", "resection_site") #假如有T<5变量应使用Fisher精确检验


table<- CreateTableOne(vars = myVars,       #条件1
                       factorVars = catVars, #条件2
                       strata = "group", # 条件3
                       data = phe_f,  #源数据
                       addOverall = TRUE)  #增加overall列
table1 <- print(table,  #构建的table函数（包括条件1.2）
                showAllLevels=TRUE, #显示所有变量
                nonnormal = nonvar) #条件3


table1<- print(table, #构建的table函数（带条件1.2.3）
               nonnormal = nonvar,#条件4
               exact = exactvars,#条件5
               catDigits = 2,contDigits = 3,pDigits = 4, #附加条件
               showAllLevels=TRUE, #显示所有变量
               quote = FALSE, # 不显示引号
               noSpaces = TRUE, # #删除用于对齐的空格
               printToggle = TRUE) #展示输出结果
table1

write.csv(table1,file="outdata/tableone_5列.csv")
read.csv(file = "outdata/tableone_5列.csv")


########提取cox多因素回归分析三线表##################

#1 cox多因素回归模型
pfs_multicox <- coxph(Surv(PFS.time, PFS) ~ group+age+gender+ajcc_T+ajcc_N+radiotherapy, data = data1)
summary(pfs_multicox) 
#2 提取HR.P.95%CI
pfs_multi1<-ShowRegTable(pfs_multicox, 
                         exp=TRUE, 
                         digits=2, 
                         pDigits =3,
                         printToggle = TRUE, 
                         quote=FALSE, 
                         ciFun=confint)
#3 提取回归系数、统计量等                     
pfs_multi2<-tidy(pfs_multicox)
pfs_multi2
#4 将两次提取结果合并
pfs_multi<-cbind(pfs_multi1, pfs_multi2)
pfs_multi
write.csv(pfs_multi,file="outdata/pfs_multi.csv")


# comparegroups制作tableone -----------------------------------------------------------
rm(list = ls())
library(compareGroups) 
phetable <- descrTable(group ~ gender + age + age_65 + race + smoke_group + ajcc_T + ajcc_N + ajcc_M + stage_group + resection_site + radiotherapy,
           data = phe_f)
export2csv(phetable, file='outdata/table1bycompare.csv')
export2word(x = phetable,file = "table1bycompare.doc") ### 转成word也很不错




######### 下面是教程
data(predimed) # 加载数据集
View(predimed) # 预览数据集
str(predimed) # 查看数据集结构
# 分类变量都显示为因子，并且都添加了标签。
# 需要知道数据集中哪些变量是分类变量，将其编码为因子，并注意是不是有序分类变量；
# 给分类变量添加标签属性，默认情况下输出的基线特征表会包含变量标签。
descrTable( ~ ., data = predimed)
# ~ 的左边为分组变量或不填变量，不填变量则计算总研究人群的基线特征，并且不进行统计检验；
# ~ 的右边为基线特征表中需要统计分析的变量，如果没填变量仅出现一个.，则默认数据集的全部变量进行统计。
descrTable(group ~ ., data = predimed) #group 为分组变量
descrTable(group ~ age + sex + smoke + waist + hormo,  # 左边为分组变量，右边为基线表行变量
           data = predimed)  # 数据集 
descrTable(group ~ . - toevent - event - diab - p14, 
           data = predimed)  
# 除了选择部分变量进行统计分析外，我们还可以选择亚组人群进行分析，比如说只选取女性进行分析。
descrTable(group ~ age + smoke + waist + hormo + toevent + event + diab + p14, 
           data = predimed,                     
           subset = sex == "Female")   
# 除了选择亚组人群外，还可以在亚组人群基础上选取特定变量进行研究。
descrTable(group ~ age + sex + smoke + waist + hormo, 
           data = predimed, 
           selec = list(hormo = sex == "Female", waist = waist > 20)) 
# 基线特征表中的变量可以在公式中出现两次，比如说bmi：
descrTable(group ~ age + sex + bmi + bmi + waist + hormo, 
           data = predimed, 
           selec = list(bmi.1 = !is.na(hormo)))  
# 输出的基线表中会报告两次bmi的统计结果，第一个bmi表示所有患者的bmi结果，第二个bmi是报告hormo变量中排除缺失值时研究患者的bmi结果。
# 指定waist为非正态分布变量，则：
descrTable(group ~ age + smoke + waist + hormo, 
           data = predimed,               
           method = c(waist = 2)) 
# 绘制分层后的基线特征表，绘制分层基线特征表的函数为strataTable()函数。
restab <- descrTable(group ~ age + smoke + bmi + waist + hormo, 
                     data = predimed)  
strataTable(restab, "sex")
## 先绘制一个基线特征表
restab <- descrTable(group ~ age + smoke + bmi + waist + hormo, 
                     data = predimed)  
restab
export2csv(restab, file='table1.csv')

#############tableone 和comparegroups的教程############
######## 用tableone 整理基线表 ####


install.packages("tableone")
dt <- read.csv("data/pheno_xena.csv")
head(dt,10)
as_tibble(dt)
library(tableone)
dput(names(dt)) # dput函数可以查看数据的变量名

# step1 创建对象
tab1 <- CreateTableOne(
  vars = c("OS", "OS.time", "gender", "age", "race", 
           "smoke", "ajcc_T", "ajcc_N", "ajcc_M", "age_group", 
           "stage2", "stage3"), # 所有要逆袭额变量
  strata = "group",  # 指定分组
  data = dt,  # 指定数据集
  factorVars = c( "gender","smoke", "site", "ajcc_T", "ajcc_N", "ajcc_M",
                  "age_group", "stage2", "stage3")  # 指定分类变量
)
# step2 生成结果
print(tab1,
      nonnormal = c("OS.time"), # 指定不符合正态性检验的连续变量
      exact = 'smoke', # 需要用fisher确切概率法的变量
      showAllLevels = TRUE, # 展示多分类变量所有 level
      pDigits = 3,   # P值小数点后的位数
      quote = TRUE  #  粘贴进入excel，分列即可。
)
# 完毕

help(package = 'tableone') # 查看帮助文档


##### 第二种方法#####
# 初次使用请先安装
install.packages('compareGroups')
library(compareGroups)
# 安装后调用
# 查看一下数据集
head(dt,10)  
descrTable(group ~.,
           data = dt,
           sd.type = 2,
           method = 4,
           show.n = T)

#####  descrTable参数设置详细解释#####
# descrTable(formula, # 公式，左边为分组，右边为变量
#            data,  # 数据集
#            subset, # 按条件筛选子集
#            selec = NA, # 按条件筛选子集
#            method = 4, # 根据数据实际情况，自动选择统计方法
#            alpha = 0.05,  # 显著性水平
#            Q1 = 0.25, Q3 = 0.75,  # 默认输出p25和p75的分位数结果
#            show.n = T, # 显示样本量
#            show.ci = F, # 显示置信区间，默认是F
#            conf.level = 0.95, # 置信区间范围
#            type = 2, # 分类变量会显示频数和百分比
#            show.p.overall = TRUE, # 显示P值
#            digits.p = 3, # p值小数点位数
#            sd.type = 2, # 1位mean（sd），2位mean±sd 
# )
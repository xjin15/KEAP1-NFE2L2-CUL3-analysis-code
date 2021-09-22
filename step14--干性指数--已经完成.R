
# 干性指数 赵队--------------------------------------------------------------------
rm(list=ls())
p_load(tidyverse,export,ggpubr)
stem <- read.csv("data/StemnessScore_ZMN.csv") %>% select(-1)
allid <- read.csv("outdata/allid_493.txt") %>% .[ ,2 ]
stem$Sample <- str_sub(string = stem$Sample,start = 1,end = 16  )
stem <- stem[stem$Annotation == "TP",]
ids <- intersect(allid, stem$Sample)

stem_493 <- stem[match(ids, stem$Sample),]
stem_493$group <- "Mut"
stem_493$group[112:493] <- "Wild"

p_load(tidyverse)
stem <- read.csv("data/CELL_for_mRNAsi.csv") 
colnames(stem)[1] <- "Sample"
allid <- read.csv("outdata/allid_493.txt") %>% .[ ,2 ]
stem$Sample <- str_sub(string = stem$Sample,start = 1,end = 16  )
stem <- stem[stem$cancer.type == "LUAD" & stem$sample.type == 1,]
ids <- intersect(allid, stem$Sample)


stem_487 <- stem[match(ids, stem$Sample),]
stem_487$group <- "Mut"
stem_487$group[match(allid[112:493],stem_487$Sample)] <- "Wild"

# 画图 ----------------------------------------------------------------------
## 赵队干性指数的图
library(ggpubr)
ggviolin(x="group",y="stemnessScore", fill = "group",data = stem_493, size = 0.1,
              palette = c("#00CCFF","#FF3333"),add = "boxplot", 
              add.params = list(group = "group"))+
  xlab("")+
  ylab("")+
  # scale_y_continuous(limits =c(0,1.5),breaks = seq(0,1.5,0.3))+ 
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y = 1.2,size = 8)+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "RNAsi Plot")+
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))+
theme(legend.text = element_text(size = 14, face = "bold"), 
      legend.title = element_text(size = 14, face = "bold")) +
  labs(x = NULL, y = NULL)

library(export)
graph2ppt(file = "output/plots/肿瘤干细胞指数对比图",append = T,)

## TCGA干性指数图
library(ggpubr)
ggviolin(x="group",y="mRNAsi", fill = "group",data = stem_487, size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(group = "group"))+
  xlab("")+
  ylab("")+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y = 0.6,size = 8)+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "RNAsi Plot")+
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))+
  theme(legend.text = element_text(size = 14, face = "bold"), 
        legend.title = element_text(size = 14, face = "bold")) +
  labs(x = NULL, y = NULL)

library(export)
graph2ppt(file = "output/plots/肿瘤干细胞指数对比图",append = T,)

## TCGA干性指数图EREG.mRNAsi
library(ggpubr)
ggviolin(x="group",y="EREG.mRNAsi", fill = "group",data = stem_487, size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(group = "group"))+
  xlab("")+
  ylab("")+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y = 1,size = 8)+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "RNAsi Plot")+
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))+
  theme(legend.text = element_text(size = 14, face = "bold"), 
        legend.title = element_text(size = 14, face = "bold")) +
  labs(x = NULL, y = NULL)

library(export)
graph2ppt(file = "output/plots/肿瘤干细胞指数对比图",append = T,)


# mDNAsi ------------------------------------------------------------------

stem <- read.csv("data/CELL_for_mDNAsi.csv") 
colnames(stem)[1] <- "Sample"
allid <- read.csv("outdata/allid_493.txt") %>% .[ ,2 ]
stem$Sample <- str_sub(string = stem$Sample,start = 1,end = 16  )
stem <- stem[stem$cancer.type == "LUAD" & stem$sample.type == 1,]
ids <- intersect(allid, stem$Sample)

stem_434d <- stem[match(ids, stem$Sample),]
stem_434d$group <- "Mut"
stem_434d$group[match(allid[112:493],stem_434d$Sample)] <- "Wild"
### 画DNAsi图
library(ggpubr)
ggviolin(x="group",y="mDNAsi", fill = "group",data = stem_434d, size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(group = "group"))+
  xlab("")+
  ylab("")+
  # scale_y_continuous(limits =c(0,1.5),breaks = seq(0,1.5,0.3))+ 
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y = 0.4,size = 8)+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "DNAsi Plot")+
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))+
  theme(legend.text = element_text(size = 14, face = "bold"), 
        legend.title = element_text(size = 14, face = "bold")) +
  labs(x = NULL, y = NULL)

library(export)
graph2ppt(file = "output/plots/肿瘤干细胞指数对比图",append = T,)









# 别人的代码 -------------------------------------------------------------------

p + geom_violin(alpha = .25, draw_quantiles = c(0.25, 0.5, 0.75)) +
geom_jitter(alpha = .5, 
            height = 0, 
            width = .25) +
  aes(color = group, fill = group) +
  scale_color_manual(values = 
                       c("blue", "red")) +
  scale_fill_manual(values = 
                      c("blue", "red")) +
  theme_bw() +
  labs(title = "关注B站UP主--小猴子_爱你呀呀") +
  labs(subtitle = "一键三连哦") +
  labs(caption = "233333")



# 小提琴图完整代码教学 ----------------------------------------------------------------

rm(list=ls()) #清除工作区
library(ggplot2)



#生成模拟数据
Group <- rep(c("group1","group2"),each=200) #组别变量
Group <- factor(Group) #组别因子化
Attribute <- c(rep("Attribute_1",100),rep("Attribute_2",100),rep("Attribute_1",100),rep("Attribute_2",100)) #每个组别的两个属性
Attribute <- factor(Attribute) #属性因子化
value <- c(rnorm(100)+1,rnorm(100)+2,rnorm(100)+1.2,rnorm(100)+1.5) #随机赋值

Data <- data.frame(Group=Group,Attribute=Attribute,value=value) #生成数据框


#对数据进行统计的函数
#指定分组变量和求值变量后，可计算出不同分组变量(或分组变量间的组合)对应的求值变量的均值，标准差，标准误，置信区间ci

#汇总数据
#计算出计数，平均值，标准差，均值的标准误差和置信区间（默认为95％）
#data：一个数据框
#measurevar：包含要汇总的变量的列的名称
#groupvars：包含分组变量的列名称的向量
#na.rm：一个布尔值，表示是否忽略NA
## conf.interval：置信区间的百分比范围（默认为95％）

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


#依据分组对vale进行统计
Data_summary <- summarySE(Data, measurevar="value", groupvars=c("Group","Attribute"))

P1 <- ggplot(Data_summary,aes(x=Group, y=value, fill=Attribute)) + #“fill=”设置填充颜色依据Attribute指定
  geom_point(aes(x=Group, y=value),pch=19,position=position_dodge(0.9),size=2.5)+ #绘制均值为点图
  geom_bar(stat = "identity",position = "dodge",alpha = 0.7) + #绘制条形图
  
  #如果误差条想表示标准差：请设置 ymin = value-sd, ymax=value+sd
  #如果误差条想表示标准误：请设置 ymin = value-se, ymax=value+se
  geom_errorbar(aes(ymin = value-ci, ymax=value+ci), #误差条表示95%的置信区间
                width=0.1, #误差条末端短横线的宽度
                position=position_dodge(0.9), 
                color="black",
                alpha = 0.7,
                size=0.5) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("") #设置x轴和y轴的标题

P1
jpeg(file = "results_Value_1.jpg",width =1600,height = 2000,units = "px",res =300) #结果保存
print(P1)
dev.off()


P2<- ggplot(Data, aes(x=Group, y=value,fill=Attribute)) + 
  geom_violin(trim=FALSE,color="white") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("") #设置x轴和y轴的标题

P2

jpeg(file = "results_Value_2.jpg",width =1600,height = 2000,units = "px",res =300) #结果保存
print(P2)
dev.off()

P3 <- ggplot(Data, aes(x=Group, y=value,fill=Attribute)) + 
  geom_violin(trim=FALSE,color="white") + #绘制小提琴图
  geom_point(data = Data_summary,aes(x=Group, y=value),pch=19,position=position_dodge(0.9),size=1.5)+ #绘制均值为点图
  geom_errorbar(data = Data_summary,aes(ymin = value-ci, ymax=value+ci), #误差条表示95%的置信区间
                width=0.1, #误差条末端短横线的宽度
                position=position_dodge(0.9), 
                color="black",
                alpha = 0.7,
                size=0.5) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("") #设置x轴和y轴的标题

P3 

jpeg(file = "results_Value_3.jpg",width =1600,height = 2000,units = "px",res =300) #结果保存
print(P3)
dev.off()

GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}

P4 <- ggplot(data=Data, aes(x=Group, y=value,fill=Attribute)) + 
  geom_split_violin(trim=FALSE,color="white") + #绘制分半的小提琴图
  geom_point(data = Data_summary,aes(x=Group, y=value),pch=19,position=position_dodge(0.9),size=1.5)+ #绘制均值为点图
  geom_errorbar(data = Data_summary,aes(ymin = value-ci, ymax=value+ci), #误差条表示95%的置信区间
                width=0.1, #误差条末端短横线的宽度
                position=position_dodge(0.9), 
                color="black",
                alpha = 0.7,
                size=0.5) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("") #设置x轴和y轴的标题

P4

jpeg(file = "results_Value_4.jpg",width =1600,height = 2000,units = "px",res =300) #结果保存
print(P4)
dev.off()


p_load(tidyverse,ggplot2, limma, export)
# 免疫浸润比较 ------------------------------------------------------------------
allgsva<-read.csv("outdata/allgsva.csv",row.names = 1) %>% t()
allgsvamut <- allgsva[ ,1:111]
allgsvawild <- allgsva[ ,112:ncol(allgsva)]
immumut <- apply(allgsvamut, 1, mean)
immuwild <- apply(allgsvawild, 1, mean)
immu <- cbind(Mut = immumut, Wild = immuwild) 
head(immu)
immu <- t(immu)
cor <- cor(immu)
cor
corrplot::corrplot(cor)
#####量化CD8+T毒性指数和IFN-γ相关毒性
mrna <- fread("outdata/mrna_fpkm_id_finish.txt") %>% column_to_rownames( "V1")
mrna <- t(mrna)
cln_survival <- read.csv("outdata/cln_xena.csv")

mrna$group <- cln_survival$group[match(x = rownames(mrna), 
                                           table = cln_survival$sample)]
mrna <- reshape2::melt(mrna,id.vars=c("group"))
mrna$group <- as.factor(mrna$group)
mrna$value <- as.numeric(mrna$value)
compare_means(value ~ group, data = mrna)
p <- ggboxplot(mrna,x="variable",y="value",fill = "group",size = 0.1,
               palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("GSVA_ES")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()
p 
p+theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6)) 
library(export)

ggviolin(data = mrna,
         x = )
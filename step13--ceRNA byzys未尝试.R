###############十.2 绘制miRNA与基因调控关系的水平网络图################
#此处需要miRNA与基因相互调控的数据，网址为http://mirwalk.umm.uni-heidelberg.de/search_mirnas/
# Transform the adjacency matrix in a long format
miR_target_filter<-miR_target[miR_target$genesymbol %in% colnames(target_heatmap),] 
connect <- miR_target_filter[c(1,3,4)]
connect$mirnaid<-gsub("hsa-","",connect$mirnaid)
names(connect)[3]<-"value"
# Number of connection per person
c( as.character(connect$mirnaid), as.character(connect$genesymbol)) %>%
  as.tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> coauth
colnames(coauth) <- c("name", "n")
dim(coauth)

# Create a graph object with igraph
mygraph <- graph_from_data_frame(connect, vertices = coauth, directed = FALSE )

# Find community
com <- walktrap.community(mygraph)
max(com$membership)

#Reorder dataset(此处需要手动修改！找到对应的列)
#coauth<-rbind(coauth[13,],coauth[c(15:35),],coauth[c(1:12),],coauth[14,],coauth[c(36:43),])#出现多的
coauth<-rbind(coauth[7:27,],coauth[c(1:6),],coauth[c(28:38),])#下降多的

# Create a graph object with igraph
mygraph <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )

# prepare a vector of n color in the viridis scale
mycolor <- colormap(colormap=colormaps$jet, nshades=3)
mycolor <- sample(mycolor, length(mycolor))
coauth$color<-c(rep(mycolor[2],21),rep(mycolor[3],17))

# Make the graph
ggraph(mygraph, layout="linear") + 
  geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  geom_node_point(aes(size=n,color = coauth$color, fill = coauth$color), alpha=0.5) +
  scale_size_continuous(range=c(0.5,8)) +
  scale_color_manual(values=mycolor) +
  geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0.4,0), "null"),
    panel.spacing=unit(c(0,0,3.4,0), "null")
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2))

#最后去掉多余的一个miRNA做热图
miR_heapmap<-miR_heapmap[,colnames(miR_heapmap) %in% coauth$name]
pheatmap(miR_heapmap,cluster_rows = F, cluster_cols = F,treeheight_col =F)
           
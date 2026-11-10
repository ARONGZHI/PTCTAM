library(monocle)


monocledat <- subset(TAM,celltype2 %in% c(as.character(unique(TAM$celltype2[grepl("TAMs",TAM$celltype2)])),"Mono"))
monocledat <- subset(monocledat,celltype2 != "SELENOP+ TAMs")
monocledat$celltype2 <- as.factor(as.character(monocledat$celltype2))

exp <- as(as.matrix(monocledat@assays$RNA@counts),'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = monocledat@meta.data)
fData <- data.frame(gene_short_name = row.names(exp), row.names = row.names(exp))
fd <- new('AnnotatedDataFrame', data = fData)

cds <- newCellDataSet(exp,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily=VGAM::negbinomial.size(),
                      lowerDetectionLimit=1)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

#####

cds <- detectGenes(cds, min_expr = 0.1)

expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10)) #保存了数据集中至少10个细胞中表达的基因

pData(cds)$Cluster  <- as.factor(pData(cds)$celltype2)

cds <- cds[expressed_genes,]

deg.cluster <- FindAllMarkers(monocledat) 

express_genes <- deg.cluster %>% dplyr::filter(p_val_adj < 0.05) %>% group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)

cds <- setOrderingFilter(cds,ordering_genes=unique(express_genes$gene))
plot_ordering_genes(cds)

# 默认使用DDRTree方法进行降维
cds <- reduceDimension(cds,max_components = 2, reduction_method = 'DDRTree',
                       #residualModelFormulaStr = "~orig.ident"
)

#save(cds,file = "RAW/cds.rda")
#load("RAW/cds.rda")

# 根据关键基因的表达模式及其变化趋势，对细胞进行排序，并据此构建出细胞发育的轨迹。
cds <- orderCells(cds)

#网上关于修复这个报错的方法有很多，基本都要修改源代码，重新安装加载，
#都比较麻烦，这里教大家一种简单的方法，只需要一行代码，就能轻松解决，
#输入并运行：
#trace('project2MST', edit = T, where = asNamespace("monocle"))

plot_cell_trajectory(cds, color_by = "State",cell_size = 0.5)  +
  facet_wrap(~State)
plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 0.8) 

plot_cell_trajectory(cds, color_by = "Cluster",cell_size = 0.1) +
  facet_wrap(~Cluster) +
  ggsci::scale_color_lancet(alpha = 0.8)


spring <-cdsspring <- c("#f6bcfd","#8dd3c6","#caa2f4","#ff7d00","#ff6581","#b05a28",
                        "#a5da6b","#e578d6","#ffd2d8","#90e4cd","#84dce0",
                        "#fe65b3","#c2e699","#a4cde1","#f6f28f","#caa2f4","#b05a28")   


#########绘图
library(viridis)
Figure4A <- plot_cell_trajectory(cds,show_branch_points = F,
                     color_by = "Pseudotime",
                     cell_size = 0.5) +
  tidydr::theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_viridis_c(option = "viridis")
Figure4A
save_plot("Figure/Figure4A.pdf",Figure4A,base_width = 5,base_height =4)

Figure4B <- plot_cell_trajectory(cds, color_by = "Cluster",cell_name_size = 1,
                     cell_size = 0.5,show_branch_points = F,state_number_size = 1,cell_link_size = 0.6,
                     show_state_number = F,show_cell_names = F) +
  tidydr::theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_manual(values = alpha(spring,0.8)) +
  geom_curve(aes(x=-4.5,y=0.2,yend=0.2,xend=1.5),
             arrow=arrow(length=unit(0.03,"npc")),
             size=0.8,curvature=0,
             color="black") +
  geom_curve(aes(x=1.9,y=1,yend= 3,xend=1.5),
             arrow=arrow(length=unit(0.03,"npc")),
             size=0.8,curvature=0,
             color="black") +
  geom_curve(aes(x=3,y=-0.8,yend=-2,xend=4.5),
             arrow=arrow(length=unit(0.03,"npc")),
             size=0.8,curvature=0,
             color="black")
dev.off()
save_plot("Figure/Figure4B.pdf",Figure4B,base_width = 6,base_height =5)


df <- cds@phenoData@data
Figure4C <- ggplot(df,aes(x=Pseudotime,color=Cluster)) +
  geom_density(size = 1) +
  theme_bw(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),  # 设置白色背景
        panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),  # 设置边框
        legend.position = c(0.2,0.67),  # 设置图例位置
        panel.grid = element_blank()  # 去除网格线
  ) +
  scale_color_manual(values = alpha(spring,0.8)) 
#ggsci::scale_color_bmj() 


save_plot("Figure/Figure4C.pdf",Figure4C,base_width = 5,base_height =4)


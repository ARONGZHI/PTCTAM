library(cowplot)

TNK <- subset(sce_mer,celltype == "T/NK")
TNK <- NormalizeData(TNK,verbose = T) 
TNK <- FindVariableFeatures(TNK,selection.method = "vst", nfeatures = 2000)
TNK <- ScaleData(TNK,verbose = FALSE,vars.to.regress = "percent.mt")
TNK <- RunPCA(TNK,verbose = T,npcs = 50)
ElbowPlot(TNK,ndims = 50)

library(harmony)
#TNK <- RunHarmony(TNK,group.by.vars = c("orig.ident"), plot_convergence = TRUE)

TNK <- RunUMAP(TNK,dims = 1:15)
TNK <- FindNeighbors(TNK,dims = 1:15)
TNK <- FindClusters(TNK,resolution = 0.6)

marker <- c("CD3E","CD3D","CD4","CD8A","NCAM1","FCGR3A")
FeaturePlot(TNK,features = marker,reduction = "umap",label = T)
##########

VlnPlot(TNK,features = marker,
        combine = TRUE,raster=FALSE)
DimPlot(TNK,label = T)

TNK <-RenameIdents(TNK,"0" = "CD8+ T","1" = "CD8+ T","2" = "CD8+ T",
                   "3" = "CD4+ T","4" = "CD4+ T","5" = "CD8+ T",
                   "6" = "CD8+ T","7" = "CD4+ T","8" = "NK",
                   "9" ="CD8+ T","10" = "CD4+ T","11" = "DNT",
                   "12"= "NK","13" = "Unknown",'14' = "Unknown", 
                   "15" = "CD8+ T",'16' = "Unknown",'17' = "Unknown",
                   '18' = "Unknown")



TNK$celltype4 <- Idents(TNK)
TNK <- subset(TNK,celltype4 != "Unknown")
DimPlot(TNK,label = T)

p21 <- DimPlot(TNK,reduction = "umap",pt.size = 0.00001,label.size = 5, label = T,
        raster = F,) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        panel.background = element_blank()) 

save_plot("Figure/FigureS2A.pdf",p21,base_height = 5,base_width =5)


p22 <- DotPlot(TNK,features = marker,
        cols = c("#D7BDE2", "#E67E22"), 
        cluster.idents = TRUE, 
        scale = TRUE) +
  theme(
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    axis.title = element_blank(),
    axis.ticks.y = element_blank())

save_plot("Figure/FigureS2B.pdf",p22,base_height = 4,base_width =7)

metadat <- sce_mer@meta.data %>%
  rownames_to_column(var = "ID") %>%
  .[,c("ID","celltype2")] %>%
  mutate(celltype2 = as.character(celltype2)) %>%
  left_join(.,TNK@meta.data %>%
              rownames_to_column(var = "ID") %>%
              .[,c("ID","celltype4")] %>%
              mutate(celltype4 = as.character(celltype4)))

metadat$celltype4[is.na(metadat$celltype4)] <- metadat$celltype2[is.na(metadat$celltype4)]
metadat <- column_to_rownames(metadat,var = "ID")

sce_mer <- AddMetaData(sce_mer,metadat)
sce_mer <- subset(sce_mer,celltype4 != "T/NK")

table(sce_mer$celltype4)



samplecell <- table(sce_mer$orig.ident,sce_mer$celltype4) %>% as.data.frame()
prop_df <- sce_mer@meta.data %>%
  group_by(orig.ident, celltype4) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup()

prop_wide <- prop_df %>%
  select(orig.ident, celltype4, proportion) %>%
  pivot_wider(names_from =  celltype4, values_from = proportion, values_fill = 0)


cor.test(prop_wide$`MARCO+ TAMs`, 
                     prop_wide$`CD8+ T`, 
                     method = "spearman")

p23 <- ggplot(prop_wide, aes(x = `MARCO+ TAMs`, y = `CD8+ T`)) +
  geom_point(size = 4.5, color = "#E74C3C", alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray80", size = 1) +
  stat_cor(method = "spearman",
           label.x = 0.6 * max(prop_wide$`MARCO+ TAMs`),
           label.y = 0.9 * max(prop_wide$`CD8+ T`),
           size = 4,
           cor.coef.name = "R2") +  
  theme_classic(base_size = 14) +
  labs(x = "MARCO+ TAM Proportion (%)",
       y = "CD8+ T Cell Proportion (%)") +
  theme(axis.text = element_text(color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1))

save_plot("Figure/FigureS2C.pdf",p23,base_height = 4.5,base_width =4.5)


inflammatory_genes <- c(
  "IL1B", "IL6", "TNF", "IL12A", "IL18", "CXCL8", "CCL2", "CCL3", "CCL4", 
  "CCL5", "CXCL10", "IFNG", "TGFB1", "IL10", "PTGS2", "NOS2"
)


TAMs <- AddModuleScore(TAM,features = list(Inflammatory = inflammatory_genes))

TAM_sub <- subset(TAM,celltype2 %in% c(as.character(unique(TAM$celltype2[grepl("TAMs",TAM$celltype2)]))))

FigureS2D <- DotPlot(TAM_sub, cols = c("#D7BDE2", "#E67E22"),
                    features = c("IL1B", "IL6", "TNF", "IL12A", "CXCL8", "CCL2", "CCL5", "CXCL10"), 
                    group.by = "celltype2") +
  theme_bw(base_size = 14) +
  labs(x=NULL,y=NULL)


FigureS2D 
save_plot("Figure/FigureS2D.pdf",FigureS2D,base_height = 4,base_width =8)

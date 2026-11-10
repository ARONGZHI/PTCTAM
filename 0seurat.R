rm(list = ls())
gc()

library(GEOquery)
library(tidyverse)
library(Seurat)
library(future)
library(future.apply)
library(scDblFinder)
library(harmony)
library(survminer)
library(survival)
library(tidyverse)
library(tableone)
library(cowplot)

sessionInfo()

GSE193581 <- getGEO("GSE193581",AnnotGPL = F,getGPL = F,destdir = "RAW")

GSE193581_pheno <- pData(GSE193581$GSE193581_series_matrix.txt.gz) %>%
  .[,c("title","source_name_ch1","disease state:ch1","treatment:ch1")]
colnames(GSE193581_pheno) <- gsub(":ch1","",colnames(GSE193581_pheno))
colnames(GSE193581_pheno)

GSE193581_pheno <- GSE193581_pheno[grepl("thyroid cancer",GSE193581_pheno$`disease state`),] %>%
  .[grepl("Thyroid",.$source_name_ch1),]  %>%
  .[grepl("untreated",.$treatment),] 


plan("multicore",workers = 70)
options(future.globals.maxSize = 400000 * 1024^2)

sce_list <- lapply(GSE193581_pheno$title, function(i){

  filepa <- file.path("/mnt/data/file/HRZ/THCA/cellcountout",i,"outs/filtered_feature_bc_matrix")
  sce <- CreateSeuratObject(counts = Read10X(filepa),project = i,
                            min.cells = 3,min.features = 300)
  tmp <- as.SingleCellExperiment(sce)
  tmp <- scDblFinder(tmp)
  tmp <- as.Seurat(tmp)
  sce <- AddMetaData(sce, tmp@meta.data[,5:8])
  sce$type <- gsub("[0-9]","",i)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce,"^MT-")
  sce[["percent.rb"]] <- PercentageFeatureSet(sce,"^RP[SL]")
  sce[["percent.hb"]] <- PercentageFeatureSet(sce,"^Hb[^(P)]")
  return(sce)
})

sce_mer <- merge(sce_list[[1]],y= sce_list[-1],add.cell.ids = GSE193581_pheno$title,project = "THCA")

VlnPlot(sce_mer,group.by = "orig.ident",features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol =3,raster = F) +
  theme(legend.position = "none")

sce_mer <- subset(sce_mer, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA < 100000 & 
                    percent.mt < 20)

VlnPlot(sce_mer,group.by = "orig.ident",features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3,raster = F) +
  theme(legend.position = "none")

sce_mer <- NormalizeData(sce_mer,verbose = T) 
sce_mer <- FindVariableFeatures(sce_mer,selection.method = "vst", nfeatures = 2000)
sce_mer <- ScaleData(sce_mer,verbose = FALSE)
sce_mer <- RunPCA(sce_mer,verbose = T,npcs = 50)

ElbowPlot(sce_mer,ndims = 50)

p1 <- DimPlot(object = sce_mer, reduction = "pca", pt.size = .1, group.by = "orig.ident",raster = F)
p2 <- VlnPlot(object = sce_mer, features = "PC_1", group.by = "orig.ident", pt.size = .1,raster = F)
CombinePlots(plots=list(p1,p2))

library(harmony)
sce_mer <- RunHarmony(sce_mer,group.by.vars = c("orig.ident"), plot_convergence = TRUE)
p3 <- DimPlot(object = sce_mer, reduction = "harmony", pt.size = 0.1, group.by = "orig.ident",raster = F)
p4 <- VlnPlot(object = sce_mer, features = "harmony_1", group.by = "orig.ident", pt.size = .1,raster = F)
CombinePlots(plots=list(p3,p4))

sce_mer <- sce_mer %>% 
  RunUMAP(reduction = "harmony",dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:20)

sce_mer <- FindClusters(sce_mer,resolution = 0.08)

#save(sce_mer,file="RAW/sce_mer.rda")
load("RAW/sce_mer.rda")
DimPlot(sce_mer,reduction = "umap",pt.size = 0.0001, label = T,raster = F)

#9查看双细胞
DimPlot(sce_mer,reduction = "umap",group.by = "scDblFinder.class",pt.size = 0.0001, label = T,raster = F)

###上皮 2,5
VlnPlot(sce_mer,features = c("EPCAM","KRT18","KRT19","KRT8"),pt.size = 0)
VlnPlot(sce_mer,features = c("EPCAM","KRT18","KRT19","KRT8",
                             "KRT5","KRT17","TFF3","TPO","SLC26A4"),pt.size = 0)
table(sce_mer$RNA_snn_res.0.14,sce_mer$type)

##T/NK cell  0,6
FeaturePlot(sce_mer,features = c("CD3D","CD3E","TRAC","CD8B","CD3G","GZMK",
                                 'NKG7','GNLY',"KLRD1",'KLRF1'),reduction = "umap",label = T,raster = F)

VlnPlot(sce_mer,features = c("CD3D","CD3E","TRAC","CD8B","CD3G",
                             "GZMK",'NKG7','GNLY',"KLRD1",'KLRF1'),pt.size = 0)

###B cell 4
##FeaturePlot(sce_mer,features = c("MS4A1","JCHAIN","CD79A","MZB1"),reduction = "umap",raster = F)
FeaturePlot(sce_mer,features = c("MS4A1","CD79A"),reduction = "umap",label = T,raster = F)
VlnPlot(sce_mer,features = c("MS4A1","CD79A"),pt.size = 0)


###myeloid cell 1
FeaturePlot(sce_mer,features = c("CD74","CD14","CD68","LYZ","CSF1R"),reduction = "umap",label = T,raster = F)
VlnPlot(sce_mer,features = c("CD14","CD68","LYZ"),pt.size = 0)
##endothelial cell  7
FeaturePlot(sce_mer,features = c("VWF","PECAM1","PLVAP"),reduction = "umap",label = T,raster = F)
VlnPlot(sce_mer,features = c("VWF","PECAM1","PLVAP"),pt.size = 0)

####mesenchymal(fibroblasts)CAFs 3,8
FeaturePlot(sce_mer,features = c("PDGFRA","DCN","COL3A1"),reduction = "umap",label = T,raster = F)
VlnPlot(sce_mer,features = c("PDGFRA","DCN","COL3A1"),pt.size = 0)
VlnPlot(sce_mer,features = c("EPCAM","KRT18","KRT19","KRT8"),pt.size = 0)
VlnPlot(sce_mer,features = c("CD14","CD68","LYZ"),pt.size = 0)


#####
sce_mer <-RenameIdents(sce_mer,"0" = "T/NK",
                       "1" = "Myeloid",
                       "2" = "Epithelial",
                       "3" = "CAFs",
                       "4" = "B",
                       "5" = "Epithelial",
                       "6" = "T/NK",
                       "7" = "Endothelial",
                       '8' = 'CAFs',
                       "9"="Doublet",
                       "10" = "T/NK")

sce_mer$celltype <- Idents(sce_mer)
DimPlot(sce_mer)

sce_mer <- subset(sce_mer,celltype != "Doublet")

#################绘图
genes <- c("EPCAM","KRT18","KRT19","KRT8","CD3D","CD3E","GZMK",'NKG7',
           "MS4A1","CD79A","CD14","CD68","LYZ","VWF","PLVAP","DCN","COL3A1")

vln.df <- sce_mer@assays$RNA@data %>%
  t() %>%
  as.data.frame()%>%
  .[,genes] %>% 
  rownames_to_column("ID") %>% 
  mutate(cluster = as.character(sce_mer$seurat_clusters))%>%
  pivot_longer(cols = 2:(ncol(.)-1),
               names_to = "gene",
               values_to = "exp") %>% 
  mutate(gene = factor(gene,levels = genes))

library(cowplot)
spring <- c("#f6bcfd","#8dd3c6","#ffc512","#ff7d00","#ff6581","#ffa300",
            "#f8d90d","#a5da6b","#e578d6","#ffd2d8","#90e4cd","#84dce0",
            "#fe65b3","#c2e699","#a4cde1","#f6f28f","#caa2f4","#b05a28")               
use_colors <- c(
  `Epithelial` = "#a4cde1",
  `T/NK`= "#ffc512",
  CAFs = "#ff7d00",
  `B`= "#ffd2d8",
  #Epithelial = "#7A0177",
  Endothelial= "#c2e699",
  Myeloid = "#caa2f4"
  #Plasma = "#636363")
)
library(scales)           
show_col(spring,ncol=5)

vln.df$cluster <- factor(vln.df$cluster,levels = c(0,1,2,3,4,5,6,7,8,10))
p1 <- ggplot(vln.df, aes(gene,exp, fill = gene,color = gene)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0),
                     position="right",
                     labels = function(x)
                       c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(cluster), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none",
        panel.spacing = unit(0, "lines"),

        panel.background = element_rect(fill = NA, color = "black"),
        plot.margin = margin(6, 6, 0, 6, "pt"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.right = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = spring) +
  scale_color_manual(values = spring)
p1

df <- data.frame(gene = genes,
                 group = c("Epithelial","Epithelial","Epithelial","Epithelial",
                           "T/NK","T/NK","T/NK","T/NK","B","B","Myeloid","Myeloid",
                           "Myeloid","Endothelial","Endothelial","CAFs","CAFs")) %>%
  mutate(gene = factor(gene,levels = genes))

dat <- data.frame(x = 1:nrow(df)) %>% mutate(group = case_when(
  x >= 1 & x <= 4 ~ "Epithelial",
  x >= 5 & x <= 8 ~ "T/NK",
  x >= 9 & x <= 10 ~ "B",
  x >= 11 & x <= 13 ~ "Myeloid",
  x >= 14 & x <= 15 ~ "Endothelial",
  x >= 16 & x <= 17 ~ "CAFs"
))

label_pos <- dat %>%
  dplyr::group_by(group) %>%
  dplyr::summarize(midpoint = mean(x))

p2 <- ggplot(df, aes(x = gene, y = 1, fill = group)) +
  geom_tile() +
  geom_text(data = label_pos,
            aes(x = midpoint, y = 1, label = group),
            fontface = "bold",
            size = 4) +
  theme_bw(base_size = 14) +
  #scale_fill_manual(values = color) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none",
        panel.spacing = unit(0, "lines"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        plot.margin = margin(0, 0, 7, 0, "pt"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),    
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank())
p2
dfy <- data.frame(cluster = c(0:8,10),
                  group = c("T/NK","Myeloid","Epithelial","CAFs",
                            "B","Epithelial","T/NK","Endothelial",
                            "CAFs","T/NK")) %>%
  mutate(cluster = factor(cluster,level = rev(cluster))) 

p3 <- ggplot(dfy, aes(x = 1, y = cluster, fill = group)) +
  geom_tile(color = "grey50") +
  #theme_bw(base_size = 14) +
  scale_fill_manual(values = use_colors) +
  #facet_grid(rows = vars(cluster), scales = "free", switch = "x") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(position = "right",labels = c("T/NK","CAFs",
                                                 "Endothelial","T/NK","Epithelial",
                                                 "B","CAFs","Epithelial","Myeloid","T/NK"),
                   expand = c(0,0)) +
  theme_cowplot(font_size = 12,line_size = 0) +
  theme(legend.position = "none",
        panel.spacing = unit(0, "lines"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        #plot.margin = margin(-1, 9, 5, 5, "pt"),
        axis.text.y = element_text(size = 14),    
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank())
p3

p4 <- plot_grid(p1, p2, ncol = 1, rel_heights = c(0.6, 0.14), align = "hv", axis = "lr")
p4
save_plot("Figure/Figure1_2.pdf",p4,base_width = 7.4,base_height = 5)

#p4
p5 <- plot_grid(p1, p3, nrow = 1, rel_widths = c(0.6, 0.14),align = "h", axis = "rl")
p5
save_plot("Figure/Figure1_3.pdf",p5,base_width = 7,base_height = 4.2)

#####
DimPlot(sce_mer,reduction = "umap",pt.size = 0.00001,label.size = 5, label = T,
              raster = F,cols = alpha(use_colors,0.1)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        panel.background = element_blank()) 
####
library(plot1cell)
circ_data <- prepare_circlize_data(sce_mer, scale = 0.8 )

cluster_colors <- alpha(rand_color(length(levels(sce_mer))),0.05)

cluster_colors <- c("#DD3497","#FA9FB5","#bcbddc","#6bAEd6","#caa2f4",
                    "#c2e699","#f9b769","#fff500")


pdf("Figure/Figure1_1.pdf",height = 4.3,width = 4.3)
plot_circlize(circ_data,do.label = T, pt.size = 0.00001, 
              col.use = alpha(c("#ffd2d8","#ff7d00","#c2e699",
                          "#a4cde1","#caa2f4","#ffc512"),0.0005),    
              bg.color = 'white',
              kde2d.n = 200, repel = F, label.cex = 1.2)
dev.off()

#########
table(sce_mer$type)
sce_mer$celltype


sce_mer$cancermajor <- ifelse(sce_mer$type == "NORM","Normal","TC") %>%
  factor(.,levels=c("Normal","TC"))

p8 <- scRNAtoolVis::cellRatioPlot(object = sce_mer,fill.col = use_colors,
                                  sample.name = "cancermajor",
                                  celltype.name = "celltype")+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1))

p8
save_plot("Figure/Figure1_4.pdf",p8,base_width = 3.8,base_height = 4)

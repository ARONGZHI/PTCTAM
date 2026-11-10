library(cowplot)

TAM <- subset(sce_mer,celltype == "Myeloid")
TAM <- NormalizeData(TAM,verbose = T) 
TAM <- FindVariableFeatures(TAM,selection.method = "vst", nfeatures = 2000)
TAM <- ScaleData(TAM,verbose = FALSE,vars.to.regress = "percent.mt")
TAM <- RunPCA(TAM,verbose = T,npcs = 50)
ElbowPlot(TAM,ndims = 50)
library(harmony)
TAM <- RunHarmony(TAM,group.by.vars = c("orig.ident"), plot_convergence = TRUE)
TAM <- TAM %>% 
  RunUMAP(reduction = "harmony",dims = 1:10) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:10)
TAM <- FindClusters(TAM,resolution = 0.45)
#load("RAW/TAM.rda")
##########

TAMMarkers <- FindAllMarkers(TAM, only.pos = TRUE,min.pct = 0.25,test.use = "t") 

TopTAMMarkers <- TAMMarkers %>% dplyr::filter(p_val_adj < 0.05) %>% group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC)

#######################
FeaturePlot(TAM,features = c("S100A9","S100A8","FCGR3A","MMP19","IL1B",
                            "C1QA","C1QC","C1QB","DAB2","SEPP1",
                            'CD1E','CD1C'),reduction = "umap",label = T)
##########DC3 9
FeaturePlot(TAM,features = c("FSCN1","LAMP3",'CCL19',"CCR7"),reduction = "umap",label = T)
VlnPlot(TAM,features = c("LAMP3",'CCL19',"CCR7"),pt.size = 0)
### DC1 TAMMarkers###  10
FeaturePlot(TAM,features = c("CLEC9A","XCR1"),reduction = "umap",label = T)
VlnPlot(TAM,features =  c("CLEC9A","XCR1"),pt.size = 0)
VlnPlot(TAM,features =  c("FCGR3A","CD14"),pt.size = 0)
### DC2 1,4
FeaturePlot(TAM,features = c('CD1E','CD1C',"CD1A","FCER1A"),reduction = "umap",label = T)
VlnPlot(TAM,features =  c('CD1E','CD1C'),pt.size = 0)
##pDC 11
FeaturePlot(TAM,features = c('LILRA4',"IRF7",'PACSIN1',"TCF4"),reduction = "umap",label = T)
VlnPlot(TAM,features =  c('LILRA4','PACSIN1'),pt.size = 0)

DimPlot(TAM,label = T)

##########
VlnPlot(TAM,features = 	 c("S100A9","S100A8","FCGR3A","MMP19","IL1B",
                           "C1QA","C1QC","C1QB","DAB2","SEPP1",
                           'CD1E','CD1C'),pt.size = 0)

FeaturePlot(TAM,features = TopTAMMarkers$gene[TopTAMMarkers$cluster==5],reduction = "umap",label = T)


FeaturePlot(TAM,features = c("S100A9","S100A8","FCGR3A","MMP19","IL1B",
                           "C1QA","C1QC","C1QB","DAB2","SEPP1",
                           'CD1E','CD1C'),reduction = "umap",label = T)

TAM <-RenameIdents(TAM,"0" = "FCGBP+ TAMs","1" = "cDC2","2" = "Mono",
                   "3" = "SELENOP+ TAMs","4" = "cDC2",
                   "5" = "IFITM2+ TAMs","6" = "MARCO+ TAMs",
                   "7" = "CXCL10+ TAMs","8" = "FTL+ TAMs","9"= "Unknow",
                   "10" = "cDC3",'11' = "cDC1", "12" = "pDC")

TAM$celltype2 <- Idents(TAM)

TAM <- subset(TAM,celltype2 != "Unknow")


p21 <- DimPlot(TAM,reduction='umap',group.by="celltype2", #split.by = 'group',
               label=T,#cols=alpha(spring,0.8),
               pt.size = 0.0001,label.size = 4.5) + 
  tidydr::theme_dr() + 
  labs(title = '') +
  theme(panel.grid=element_blank(),
        legend.text = element_text(size=10))
p21
cowplot::save_plot("Figure/Figure2_1.pdf",p21,base_width = 6,base_height = 4)


###################无用
tamgenes <- c("LAMP3","CCR7","CLEC9A","XCR1",
              'CD1E','CD1C','FCER1A',
              "C1QA","C1QC","C1QB",
              'LILRA4','PACSIN1',
              'FCGBP',"NCKAP5","KCNQ3",
              "FOLR2","SELENOP","PLTP",
              "FCN1","VCAN","THBS1",
              "PREX1","IFITM2",'FN1',
              'MARCO',"CCL18","GPNMB",
              'CXCL10',"APOBEC3A","RSAD2",
              "FTL","CTSD","IFI27")

##
p22 <- DotPlot(TAM,features=tamgenes,group.by ="celltype2",cluster.idents = T,scale = T)+
  RotatedAxis() +
  theme(
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45),
    axis.ticks.y = element_blank())
p22

############################################
################

###########
Celldata <- data.frame(celltype = TAM$celltype2,cluster = as.character(TAM$seurat_clusters)) %>%
  dplyr::filter(!duplicated(cluster))


TAMmarkernew <- left_join(TAMMarkers,Celldata) %>%
  dplyr::filter(!is.na(celltype))

TAMMarkerssig <- TAMmarkernew[TAMmarkernew$p_val_adj < 0.05,] %>% 
  dplyr::select(celltype,gene,pct.1,pct.2,avg_log2FC,p_val,p_val_adj) 

#write.csv(TAMMarkerssig,file = "Figure/TAMMarkers.csv",row.names = F,quote = F)
########
###
df <- TAMmarkernew[grepl("TAMs",TAMmarkernew$celltype) ,]
df$point <- ifelse(df$p_val<0.05, "color1", "color2")

df_tile <- data.frame(x=unique(df$celltype), y=0)
df_bar <- df %>% dplyr::group_by(celltype)%>% dplyr::summarise(max = max(avg_log2FC))
df_label <- df %>% group_by(celltype) %>% top_n(n=10,wt=avg_log2FC)

point_color <- c("#cc3333","#0099cc")

p23 <- ggplot() + 
  geom_col(data = df_bar, aes(x = celltype,y = max),fill = "#dcdcdc",alpha=0.5) + 
  geom_tile(data = df_tile, aes(x = x,y = 0,fill = x ), height = 0.5, show.legend = F) +
  scale_fill_manual(values = spring)+
  geom_jitter(data = df, aes(x = celltype, y = avg_log2FC, color = point), size = 0.7, width =0.4)+
  scale_color_manual(values = point_color)+
  ggrepel::geom_text_repel(data=df_label,aes(x=celltype, y=avg_log2FC, label=gene),force = 0.5,size = 5) +
  geom_text(data = df_tile, aes(x=x,y=0,label=x),size =6) +
  labs(x=NULL,y="Average log2FoldChange") +
  theme_classic() +
  theme(axis.line.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(), #不展示x轴的刻度、轴线等各种东西
        axis.text.y = element_text(size = 13),   #axis.text.y更改纵坐标/y轴刻度字体大小
        axis.title.y = element_text(size = 16),   #axis.text.y更改纵坐标/y轴标题字体大小
        legend.position="none"   ) 

save_plot("Figure/Figure2_3.pdf",p23,base_width = 14,base_height = 5)


####################

p24 <- scRNAtoolVis::cellRatioPlot(object = TAM,
                            sample.name = "cancermajor",
                            celltype.name = "celltype2")+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1))

p24
save_plot("Figure/Figure2_4.pdf",p24,base_width = 4,base_height = 4)

########
Macmarker <- FindAllMarkers(TAMsubmarc, only.pos = F,min.pct = 0.25) 

Macmarker_deg <- split(Macmarker,f = Macmarker$cluster)

Macmarker_deg <- lapply(Macmarker_deg, function(i){
  res <- i %>%
    dplyr::filter(p_val_adj < 0.05)
  return(res)
})

library(fgsea)
hallmark.list <- clusterProfiler::read.gmt("../Supp/h.all.v2024.1.Hs.symbols.gmt") %>%
  split(.$term) %>% lapply( "[[", 2)


gseares <- lapply(Macmarker_deg,function(i){
  res <- i
  id <- res$avg_log2FC
  names(id) <- res$gene
  id <- id[order(id,decreasing = T)]
  
  gsea_re <- fgsea(pathways = hallmark.list,
                   stats = id,
                   minSize=1,
                   maxSize=1000,
                   nperm = 1000)
  return(gsea_re)
})

gsepval <- lapply(names(gseares), function(i) {
  retmp <- gseares[[i]] %>% dplyr::filter(pval < 0.05) %>%
    mutate(cluster = i)
})

gsepval <- do.call(rbind,gsepval) %>%
  mutate(pathway =gsub("HALLMARK_","",pathway)) %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = pval)

p3 <- ggplot(gsepval,aes(x=-log10(pval),y=pathway,fill=cluster)) +
  geom_bar(stat = "identity") +
  facet_wrap(~cluster,scales = "free_y",ncol = 1,strip.position = "right") +
  scale_fill_manual(values =spring) +
  theme_cowplot(font_size = 13) +
  labs(y=NULL,x="-Log10(P Value)") +
  theme(legend.position = "none") 
p3
save_plot("Figure/Figure2_6.pdf",p3,base_width = 6,base_height = 10)


library(ggpubr)

samplecell <- table(TAM$orig.ident,TAM$celltype2) %>% as.data.frame()
bar_per <- samplecell %>%
  dplyr::group_by(Var1) %>%
  dplyr::mutate(sum(Freq)) %>%
  mutate(percent = Freq /`sum(Freq)`) %>%
  dplyr::filter(Freq != 0)

colnames(bar_per)[1] <- "ID"
bar_per <- bar_per %>%
  mutate(cancertype = gsub("[0-9]","",ID),
         percent = percent*100) %>%
  mutate(cancertype2 = ifelse(cancertype=="NORM","Normal","TC"),
         cancertype = ifelse(cancertype=="NORM","Normal",cancertype)) %>%
  mutate(cancertype = factor(cancertype,levels = c("Normal","PTC","ATC")))


p29 <- ggplot(bar_per[grepl("TAMs",bar_per$Var2),],aes(x= Var2,y=percent,fill = cancertype2)) +
  geom_boxplot(outliers = F) +
  scale_y_continuous(limits = c(0,40)) +
  labs(y = "Cell Ratio") +
  ggsci::scale_fill_aaas(alpha = 0.8) +
  stat_compare_means(label = "p.signif",label.y = 38,size=5,method = "t.test",
                     color = c("black","red",'black',
                               'red',"black","black")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "top",legend.title = element_blank())


p29
save_plot("Figure/Figure2_5.pdf",p29,base_width = 5,base_height = 4)


library(tidyverse)
library(Seurat)
library(parallel)
library(cowplot)
options(future.globals.maxSize = 400000 * 1024^2)

gtf <- rtracklayer::import("RAW/Homo_sapiens.GRCh38.110.gtf.gz") %>%
  as.data.frame() %>%
  .[,c("gene_id","type","gene_name","gene_biotype")] %>%
  dplyr::filter(type =="gene")  %>%
  dplyr::filter(!is.na(gene_name))

survdat <- data.table::fread("RAW/TCGA-THCA.survival.tsv.gz") %>%
  as.data.frame() %>%
  dplyr::filter(!duplicated(sample)) %>%
  dplyr::filter(!is.na(OS)) %>%
  dplyr::filter(!is.na(OS.time)) 

rownames(survdat) <- NULL

bulkexp <- data.table::fread("RAW/TCGA-THCA.star_tpm.tsv.gz") %>%
  as.data.frame()

Normalsample <- colnames(bulkexp)[as.numeric(substr(colnames(bulkexp),14,15)) > 10] %>% 
  .[!is.na(.)] 

Primarysample <- colnames(bulkexp)[substr(colnames(bulkexp),14,15) == "01"] %>%
  .[!is.na(.)] %>% 
  .[substr(.,1,15) %in% survdat$sample] 


load("RAW/BisqueRNA.rda")
bioque_estimates <- res$bulk.props %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  mutate(sample = substr(sample, 1,15)) %>%
  inner_join(survdat[,c(1:3)],.)

dat <- bioque_estimates

dat$Group <- ifelse(dat$`MARCO+ TAMs` > median(dat$`MARCO+ TAMs`),"High MARCO+ TAMs","Low MARCO+ TAMs") %>%
  factor(.,levels=c("Low MARCO+ TAMs","High MARCO+ TAMs"))

#########差异性分析

bulkexp <- data.table::fread("RAW/TCGA-THCA.star_tpm.tsv.gz") %>%
  as.data.frame() %>%
  .[,c("Ensembl_ID",Normalsample,Primarysample)] %>%
  mutate(`Ensembl_ID` = str_split(`Ensembl_ID`,pattern = "\\.",simplify = T)[,1]) %>%
  dplyr::rename("gene_id" =`Ensembl_ID`) %>%
  inner_join(gtf[,c(1,3,4)],.) %>%
  dplyr::filter(!duplicated(gene_name)) %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  dplyr::select(-"gene_id") %>%
  dplyr::select(-"gene_biotype") %>%
  column_to_rownames("gene_name") %>%
  2^.-1 
colnames(bulkexp) <- substr(colnames(bulkexp),1,15)

exp <- bulkexp[,dat$sample] 
exp <- exp[apply(exp, 1, function(x) sum(x > 1) > ncol(exp)*0.25), ]

###############
group <- data.frame(barcode = dat$sample,
                    condition = ifelse(dat$Group =="High MARCO+ TAMs","High","Low") %>%
                      factor(.,levels=c("High","Low")))

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = round(exp),
                              colData = group,
                              design = ~condition)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "High", "Low"))
#save(res,file="RAW/res.rda")
#load("RAW/res.rda")
#write.csv(res,file = "Figure/diffres.csv",quote = F,row.names = T)

sel_res <- res %>%
  as.data.frame() %>%
  mutate(class = ifelse(padj < 0.05 & log2FoldChange > 1,"up",
                        ifelse(padj < 0.05 & log2FoldChange < -1,"down","ns"))) 
table(sel_res$class)

p51 <- ggplot(sel_res,aes(y=-log10(padj),x=log2FoldChange)) +
  geom_point(mapping=aes(color=class),alpha = 0.6, shape = 16, size = 2.5) +
  geom_vline(xintercept = c(-0.8,0.8), linetype ="dashed") +
  geom_hline(yintercept = -log10(0.05), linetype ="dashed") +
  scale_color_manual(values = c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey"),
                       labels = c('down 156', 'ns 13633', "up 41")) +
  labs(x = "log2(fold change)", y = "-log10(adjusted P-value)", colour = "Expression change") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold", color = "black", size = 16),
        axis.text = element_text(color = "black", size = 16, face = "bold"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", color = "black", size = 14),
        legend.spacing.x = unit(0, "cm"),
        legend.position = c(0.2, 0.89)) # 设置图例位置为右上角
p51
save_plot("Figure/FigureS5_1.pdf",p51,base_height = 5.5,base_width = 4.5)



sel_res2 <- sel_res %>% dplyr::filter(padj < 0.05 & abs(log2FoldChange) >1)


#####
datres <- data.table::fread("RAW/TCGA-THCA.star_tpm.tsv.gz") %>%
  as.data.frame() %>%
  .[,c("Ensembl_ID",Primarysample)] %>%
  mutate(`Ensembl_ID` = str_split(`Ensembl_ID`,pattern = "\\.",simplify = T)[,1]) %>%
  dplyr::rename("gene_id" =`Ensembl_ID`) %>%
  inner_join(gtf[,c(1,3)],.) %>%
  dplyr::filter(!duplicated(gene_name)) %>%
  dplyr::select(-"gene_id") %>%
  column_to_rownames("gene_name") 

datres <- datres[rownames(sel_res2),] %>%
  t() %>%
  as.data.frame( ) %>%
  rownames_to_column(var = "sample") %>%
  inner_join(survdat[,c(1:3)],.) 



genecox <- lapply(colnames(datres)[4:ncol(datres)], function(i){
  x <- datres[,i]
  cx <- survival::coxph(survival::Surv(datres$OS.time/365, datres$OS ==1)~x)
  cx <- summary(cx)
  cHR_lower <- signif(cx$conf.int[,"lower .95"], 3)
  cHR_upper <- signif(cx$conf.int[,"upper .95"],3)
  cHR <- cx$conf.int[,"exp(coef)"]
  cpValue <- cx$coefficients[,"Pr(>|z|)"]
  
  res <- data.frame(cHR_lower = cHR_lower, cHR_upper =cHR_upper,
                    CoxHR = cHR,CoxpValue = cpValue)
  return(res)
})


names(genecox) <-  colnames(datres)[4:ncol(datres)]

genecox <- do.call(rbind,genecox) %>% 
  mutate(CoxHR = round(CoxHR,2))

write.csv(genecox,file = "Figure/genecox.csv",quote = F,row.names = T)

genecoxsig  <- genecox  %>%
  dplyr::filter(CoxpValue < 0.05)

datres2 <- datres %>%
  dplyr::rename("ID" = sample) %>%
  .[,c("ID",'OS.time','OS',rownames(genecoxsig))]

save(datres2,file = "RAW/surda.rda")

#load("RAW/surda.rda")

set.seed(11)
library(caret)
library(Mime1)
library(glmnet)
library(aplot)

trainIndex <- createDataPartition(datres2$OS, p = 0.7, # 划分比例
                                  list = FALSE)

datTrain <- datres2[ trainIndex,] # 训练集
table(datTrain$OS)
datTest  <- datres2[-trainIndex,] 
table(datTest $OS)
list_Data <- list(Train = datTrain,Test = datTest)
save(list_Data,file = "sur_datalist.rda")
#load("sur_datalist.rda")
list_Data  <- lapply(list_Data,function(x){
  colnames(x) = gsub('-','.',colnames(x))
  return(x)})

list_Data <- lapply(list_Data,function(x){
  x[,-c(1:3)] <- apply(x[,-c(1:3)],2,as.numeric)
  return(x)})

list_Data <- lapply(list_Data,function(x){
  x[,c(2:3)] <- apply(x[,c(2:3)],2,as.numeric)
  return(x)})

candidate_genes <- colnames(datTrain)[4:36]
candidate_genes <- gsub('-','.',candidate_genes)

library(survival)

set.seed(12)
x1 <- as.matrix(list_Data$Train[,-c(1:3)])
x2 <- as.matrix(Surv(list_Data$Train$OS.time,list_Data$Train$OS))

fit <- glmnet(x1, x2, family = "cox",alpha=1)

pdf("Figure/FigureS5_2.pdf",height = 4,width = 3.5)
plot(fit,
     xvar='lambda',#'norm'是默认值，为L1 Norm；'lambda'：log-lambda;'dev':模型解释的%deviance
     label=T)
dev.off()

cv.fit <- cv.glmnet(x1, x2,
                 nfold = 10, #10-fold cross-validation
                 family = "cox", alpha = 1,
                 type.measure = "class")
pdf("Figure/FigureS5_3.pdf",height = 4,width = 3.5)
plot(cv.fit,
     label=T)
dev.off()

myCoefs <- coef(cv.fit, s = "lambda.min")


Lassogenes <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]

#################
source("R/ML.Dev.Prog.Sig.R")

MLres <- ML.Dev.Prog.Sig(train_data = list_Data$Train,
                         list_train_vali_Data = list_Data,
                         unicox.filter.for.candi = F,
                         unicox_p_cutoff = 0.05,
                         candidate_genes = Lassogenes,
                         direction_for_stepcox = "both",
                         mode = 'all',seed = 11,nodesize = 5)
#save(MLres,file = 'RAW/MLres.rda')
cindex <- MLres[["Cindex.res"]]
#load("RAW/MLres.rda")
library(c)
ML <- cindex %>%
  mutate(Cindex = round(Cindex,3))
write.csv(ML,file = "Figure/ML.txt",row.names = F,quote = F)

dim(cindex)[1]/2

gc()
#########
set.seed(11)
stepfit <- step(coxph(Surv(OS.time,OS)~., 
                      list_Data$Train[,c("OS.time","OS",Lassogenes)]),
                direction = "both")
coef(stepfit)


##top50
cindex <- cindex[1:60,]

cindex$Cindex <- as.numeric(sprintf("%.3f", cindex$Cindex))
mean_cindex <- aggregate(
  x = cindex$Cindex,
  by = list(cindex$Model),
  FUN = mean) %>%
  dplyr::rename(Model = Group.1,
                mean = x) %>%
  mutate(mean = as.numeric(sprintf("%.3f", mean))) %>%
  .[order(.$mean,decreasing = F),] %>%
  mutate(Value = "Mean C-index in all cohorts") %>%
  as.data.frame()


labels <- as.data.frame(unique(cindex$ID))
colnames(labels) <- "Cohort"

labels$Cohort <- factor(labels$Cohort, levels = c("Train", "Test"))
cindex$ID <- factor(cindex$ID, levels = c("Train", "Test"))



cindex$Model <- factor(cindex$Model, levels = mean_cindex$Model) 

mean_cindex$Model <- factor(mean_cindex$Model, levels = mean_cindex$Model)

######



p1 <- ggplot(cindex, aes(x = ID, y = Model)) +
  geom_tile(aes(fill = Cindex), color = "white", size = 0.5) +
  geom_text(aes(label = Cindex), vjust = 0.5, color = "black", size = 3) +
  scale_fill_viridis_c(alpha = 0.7,name = "C-index") +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_blank()
  )
p1

p2 <- ggplot(labels, aes(Cohort, y = 1)) +
  geom_tile(aes(fill = Cohort), color = "white", size = 0.5) +
  ggsci::scale_fill_aaas(alpha = 0.7) +
  #scale_fill_manual(name = "Cohort") +
  theme_void()
p2

p3 <- ggplot(data = mean_cindex, aes(x = Model, y = mean, fill = Value)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = alpha(ggsci::pal_aaas()(6)[3],0.7)) +
  geom_text(aes(label = mean), position = position_dodge(width = 0.9), 
            vjust = 0.5, hjust = 1.2, size =3) +
  theme(
    axis.title = element_text(size = 14),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    # panel.border = element_rect(colour = "black", fill=NA,size = 0.2),
    panel.grid = element_blank(),
    # legend.position = "",
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  ) +
  # scale_y_continuous(#position = "right",
  #   expand = c(0, 0.001))+
  labs(y = "") +
  coord_flip()

p3
p52 <- p1 %>%
  insert_top(p2, height = 0.03) %>%
  insert_right(p3, width = 0.3)

p52
save_plot("Figure/Figure5_2.pdf",p52,base_height = 8,base_width = 7)


####
RS <- MLres$riskscore$`StepCox[both]`

library(timeROC)
ROCTrain <-timeROC(T=RS$Train$OS.time/365,
                   delta=RS$Train$OS,
                   marker=RS$Train$RS,#计算ROC的biomaker，默认是marker值越大，事件越可能发生；反之的话，前面加-号。
                   cause=1, #所关心的事件结局
                   weighting="cox",#COX模型
                   times=c(3,5),#想计算的ROC曲线的时间节点
                   iid=F)

ROCTest <-timeROC(T=RS$Test$OS.time/365,
                  delta=RS$Test$OS,
                  marker=RS$Test$RS,#计算ROC的biomaker，默认是marker值越大，事件越可能发生；反之的话，前面加-号。
                  cause=1, #所关心的事件结局
                  weighting="cox",#COX模型
                  times=c(3,5),#想计算的ROC曲线的时间节点
                  iid=F)


ROCdata <- data.frame(Year = c("3 years-ROC","5 years-ROC","3 years-ROC","5 years-ROC"),
                      Value = c(ROCTrain$AUC[1],ROCTrain$AUC[2],ROCTest$AUC[1],ROCTest$AUC[2]),
                      Cohort =  c("Train","Train","Test","Test")) %>%
  mutate(Value = round(Value,3))



ROCdata <- data.frame(Year = c("5 years-ROC","5 years-ROC","5 years-ROC"),
                      Value = c(ROCTrain$AUC[2],0.68,0.878),
                      Cohort =  c("MTAMRS signature","Cuproptosis-Related signature (Tian)",
                                  "Lipid Metabolism-Related signature (Liao)")) %>%
  mutate(Value = round(Value,3))


p53 <- ggplot(data = ROCdata, aes(x = Cohort, y = Value, fill = Cohort)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(label = Value),
            position = position_dodge(width = 0.9), 
            vjust = 1.5, hjust = 0.5, size =4) +
  #coord_flip() +
  #scale_fill_viridis_d() +
  ggsci::scale_fill_aaas(alpha = 0.6) +
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1)) +
  theme(legend.position = "none",
        axis.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
    # legend.position = "",
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white")
  ) +
  # scale_y_continuous(position = "right")+
  labs(y = "AUC", x = "", title = "") +
  facet_grid(. ~ Year) +
  coord_flip()
p53
save_plot("Figure/Figure5_3.pdf",p53,base_height = 3,base_width = 7)


#########
RS$Train$Group <- ifelse(RS$Train$RS > median(RS$Train$RS), 
                         "High Risk Score", "Low Risk Score")

sfit1 <- survfit(Surv(OS.time/365, OS) ~ Group, data = RS$Train)

library(survminer)

p54 <- ggsurvplot(sfit1,palette = c("#FF8C00","#1E90FF"),
                          legend.title = element_blank(),
                          legend.labs = c("High Risk Score","Low Risk Score"), 
                          pval = T,risk.table = F,
                          legend.lab = c("High Risk Score","Low Risk Score"),
                          ggtheme = theme_bw(),
                          xlab = "Time in years",ylab = "Overall Survival Probility")

p54  <- p54$plot + scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(title = "Train Cohort") +
  theme_cowplot() +
  theme(legend.position = c(0.5,0.6),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=12))
p54

save_plot("Figure/Figure5_4.pdf",p54,base_width = 4,base_height =4)


RS$Test$Group <- ifelse(RS$Test$RS > median(RS$Train$RS), 
                         "High Risk Score", "Low Risk Score")

sfit2 <- survfit(Surv(OS.time/365, OS) ~ Group, data = RS$Test)

p55 <- ggsurvplot(sfit2,palette = c("#FF8C00","#1E90FF"),
                  legend.title = element_blank(),
                  legend.labs = c("High Risk Score","Low Risk Score"), 
                  pval = T,risk.table = F,
                  legend.lab = c("High Risk Score","Low Risk Score"),
                  ggtheme = theme_bw(),
                  xlab = "Time in years",ylab = "Overall Survival Probility")

p55  <- p55$plot + scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(title = "Test Cohort") +
  theme_cowplot() +
  theme(legend.position = c(0.5,0.6),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=12))
p55

save_plot("Figure/Figure5_5.pdf",p55,base_width = 4,base_height =4)

datatra <- cbind(ROCTrain$FP %>% as.data.frame() %>%
                   dplyr::rename(`3 years-ROC: 0.901` = `t=3`,
                                 `5 years-ROC: 0.936` = `t=5`) %>% 
                   data.table::melt() %>%
                   dplyr::rename(FP = value),
                 ROCTrain$TP %>% as.data.frame() %>%
                   dplyr::rename(`3 years-ROC: 0.90` = `t=3`,
                                 `5 years-ROC: 0.936` = `t=5`) %>% 
                   data.table::melt() %>%
                   dplyr::rename(TP = value)) %>% .[,c(1,2,4)]

p56 <- ggplot(datatra,aes(FP,TP,col=variable))+
  geom_line() +
  labs(title = "Train Cohort") +
  theme_bw() +
  labs(x="False Positive Rate",y="True Positive Rate") +
  geom_abline(intercept = 0, slope = 1,color="grey50") +
  ggsci::scale_color_aaas() +
  theme(legend.title = element_blank(),
        legend.position = c(0.7,0.2),
        panel.grid = element_blank())
p56
save_plot("Figure/Figure5_6.pdf",p56,base_width = 4.2,base_height =4.2)


round(ROCTest$AUC[2],3)
datatratest <- cbind(ROCTest$FP %>% as.data.frame() %>%
                       dplyr::rename(`3 years-ROC: 0.944` = `t=3`,
                                     `5 years-ROC: 0.929` = `t=5`) %>% 
                       data.table::melt() %>%
                       dplyr::rename(FP = value),
                     ROCTest$TP %>% as.data.frame() %>%
                       dplyr::rename(`3 years-ROC: 0.944` = `t=3`,
                                     `5 years-ROC: 0.929` = `t=5`) %>% 
                       data.table::melt() %>%
                       dplyr::rename(TP = value)) %>% .[,c(1,2,4)]

p57 <- ggplot(datatratest,aes(FP,TP,col=variable))+
  geom_line() +
  labs(title = "Test Cohort") +
  labs(x="False Positive Rate",y="True Positive Rate") +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1,color="grey50") +
  ggsci::scale_color_aaas() +
  theme(legend.title = element_blank(),
        legend.position = c(0.7,0.2),
        panel.grid = element_blank())

save_plot("Figure/Figure5_7.pdf",p57,base_width = 4.2,base_height =4.2)


##################end
##
RSall <- rbind(RS$Train,RS$Test)

multidat <- data.table::fread("RAW/TCGA-THCA.clinical.tsv.gz") %>%
  as.data.frame() %>%
  dplyr::filter(!duplicated(sample)) %>%
  dplyr::filter(sample %in% Primarysample) %>%
  .[,c("sample","age_at_diagnosis.diagnoses","gender.demographic",
       "ajcc_pathologic_t.diagnoses","ajcc_pathologic_n.diagnoses",
       "ajcc_pathologic_m.diagnoses","ajcc_pathologic_stage.diagnoses")] %>%
  mutate(`age_at_diagnosis.diagnoses` = `age_at_diagnosis.diagnoses`/365,
         `ajcc_pathologic_stage.diagnoses` = gsub('Stage ',"",`ajcc_pathologic_stage.diagnoses`),
         `ajcc_pathologic_t.diagnoses` = gsub('[a|b]',"",`ajcc_pathologic_t.diagnoses`),
         `ajcc_pathologic_n.diagnoses` = gsub('[a|b]',"",`ajcc_pathologic_n.diagnoses`),
         `ajcc_pathologic_m.diagnoses` = gsub('[a|b]',"",`ajcc_pathologic_m.diagnoses`),
         `ajcc_pathologic_stage.diagnoses` = gsub("[A|C]","",`ajcc_pathologic_stage.diagnoses`)) %>%
  dplyr::rename("age" = `age_at_diagnosis.diagnoses`,
                "gender"  = gender.demographic,
                "pathologic_stage" = `ajcc_pathologic_stage.diagnoses`,
                "pathologic_t" = `ajcc_pathologic_t.diagnoses`,
                "pathologic_n" = `ajcc_pathologic_n.diagnoses`,
                "pathologic_m" = `ajcc_pathologic_m.diagnoses`) %>%
  inner_join( RSall %>% dplyr::rename(sample = ID),.) 

multidat$pathologic_stage[multidat$pathologic_stage =="IVA"] <- "IV"
multidat$pathologic_stage[multidat$pathologic_stag2 =="IVC"] <- "IV"
multidat$pathologic_t[multidat$pathologic_t =="TX"] <- NA
multidat$pathologic_n[multidat$pathologic_n =="NX"] <- NA
multidat$pathologic_m[multidat$pathologic_m =="MX"] <- NA

multidat[c(6:11)] <- lapply(multidat[,c(6:11)],as.factor) 
multidat[c(6:11)] <- lapply(multidat[,c(6:11)],as.numeric) 


multifit <- coxph(Surv(OS.time/365,OS)~age + gender  + pathologic_t + pathologic_n + pathologic_m +
                    pathologic_stage + RS,
                  data = multidat)

summary(multifit)
multifitres <- data.frame(Feature = rownames(summary(multifit)$coefficients),
                          HR = signif(as.matrix(summary(multifit)$coefficients)[ , 2], 2),
                          low = signif(summary(multifit)$conf.int[ , 3], 5),
                          high = signif(summary(multifit)$conf.int[ , 4], 5),
                          pvalue = signif(as.matrix(summary(multifit)$coefficients)[ , 5], 2)) %>%
  mutate(Index = paste(rep(" ",12),collapse = " "),"P value" = ifelse(pvalue < 0.05,paste0(pvalue,"*"),pvalue),
         `HR (95%CI)` = paste0(HR," (",low,
                     " - ",high,")"))
multifitres$Feature <- c("Age","Gender","Tumor size",
                         "Node","Metastasis","Stage","Risk score")

multifitres
colnames(multifitres) <- c("Characteristic","HR","lor","uor","P","","P Value","HR (95%CI)")

tm <- forest_theme(ci_col = "#008B45FF",refline_col = "#EE0000FF")

pdf(file = "Figure/FigureS5_4.pdf",width = 6,height = 4)

p1 <- forest(multifitres[,c(1,8,6,7)],est = multifitres$HR,
       lower = multifitres$lor,
       upper = multifitres$uor,ci_column = 3,xlim = c(-1,6),
       ticks_at = c(-0.5,1,5),ref_line = 1,theme = tm)
edit_plot(p1,row = which(multifitres$P <0.05),
          gp = grid::gpar(col = "red", fontface = "bold"))
dev.off()


p1


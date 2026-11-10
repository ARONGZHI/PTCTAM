library(tableone)
####
#download.file("https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-THCA.survival.tsv.gz",
#              destfile = "RAW/TCGA-THCA.survival.tsv.gz")

#download.file("https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-THCA.clinical.tsv.gz",
#              destfile = "RAW/TCGA-THCA.clinical.tsv.gz")

#download.file("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FTHCA_survival.txt",
#              destfile = "RAW/TCGA-THCA.survival.txt")

#download.file("https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-THCA.star_tpm.tsv.gz",
#              destfile = "RAW/TCGA-THCA.star_tpm.tsv.gz")

#download.file("https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-THCA.star_counts.tsv.gz",
#              destfile = "RAW/TCGA-THCA.star_counts.tsv.gz")


gtf <- rtracklayer::import("RAW/Homo_sapiens.GRCh38.110.gtf.gz") %>%
  as.data.frame() %>%
  .[,c("gene_id","type","gene_name","gene_biotype")] %>%
  dplyr::filter(type =="gene")  %>%
  dplyr::filter(!is.na(gene_name))

survdat <- data.table::fread("RAW/TCGA-THCA.survival.tsv.gz") %>%
  as.data.frame() %>%
  dplyr::filter(!duplicated(sample))

rownames(survdat) <- NULL

bulkexp <- data.table::fread("RAW/TCGA-THCA.star_tpm.tsv.gz") %>%
  as.data.frame()

Normalsample <- colnames(bulkexp)[as.numeric(substr(colnames(bulkexp),14,15)) > 10] %>% 
  .[!is.na(.)] 

Primarysample <- colnames(bulkexp)[substr(colnames(bulkexp),14,15) == "01"] %>%
  .[!is.na(.)] %>% 
  .[. %in% survdat$sample]

bulkexp <- bulkexp[,c("Ensembl_ID",Normalsample,Primarysample)] %>%
  mutate(`Ensembl_ID` = str_split(`Ensembl_ID`,pattern = "\\.",simplify = T)[,1]) %>%
  dplyr::rename("gene_id" =`Ensembl_ID`) %>%
  inner_join(gtf[,c(1,3)],.) %>%
  dplyr::filter(!duplicated(gene_name)) %>%
  dplyr::select(-"gene_id") %>%
  column_to_rownames("gene_name") %>%
  2^.-1 

bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(bulkexp[,c(Normalsample,Primarysample)]))

metadat <- sce_mer@meta.data %>%
  rownames_to_column(var = "ID") %>%
  #.[,c("ID","celltype")] %>%
  mutate(celltype = as.character(celltype)) %>%
  left_join(.,TAM@meta.data %>%
               rownames_to_column(var = "ID") %>%
               .[,c("ID","celltype2")] %>%
            mutate(celltype2 = as.character(celltype2)))

metadat$celltype2[is.na(metadat$celltype2)] <- metadat$celltype[is.na(metadat$celltype2)]
metadat <- column_to_rownames(metadat,var = "ID")

sce_mer <- AddMetaData(sce_mer,metadat)
sce_mer <- subset(sce_mer,celltype2 != "Myeloid")

tamcounts_matrix <-  as.matrix(sce_mer@assays$RNA@counts)

sample_ids <- colnames(tamcounts_matrix) 
individual_labels <- sce_mer@meta.data$orig.ident
celltype_labels <- as.character(sce_mer$celltype2)

sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=sample_ids,
                       SubjectName=individual_labels,
                       cellType=celltype_labels)

sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.pdata <- new("AnnotatedDataFrame",data=sc.pheno,varMetadata=sc.meta)

sc.eset <- Biobase::ExpressionSet(assayData=tamcounts_matrix,phenoData=sc.pdata)

seed <- 123

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset,markers=NULL, use.overlap=F)
#save(res,file = "RAW/BisqueRNA.rda")
#load(file = "RAW/BisqueRNA.rda")

bioque_estimates <- res$bulk.props %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  inner_join(survdat[,c(1:3)],.) 

pheno <- data.table::fread("RAW/TCGA-THCA.clinical.tsv.gz") %>%
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
         `ajcc_pathologic_m.diagnoses` = gsub('[a|b]',"",`ajcc_pathologic_m.diagnoses`)) %>%
  dplyr::rename("age" = `age_at_diagnosis.diagnoses`,
              "gender"  = gender.demographic,
         "pathologic_stage" = `ajcc_pathologic_stage.diagnoses`,
         "pathologic_t" = `ajcc_pathologic_t.diagnoses`,
         "pathologic_n" = `ajcc_pathologic_n.diagnoses`,
         "pathologic_m" = `ajcc_pathologic_m.diagnoses`) %>%
  inner_join(survdat[,c(1:3)],.) %>%
  mutate(Survival = ifelse(OS== 1,"Dead","Alive")) %>%
  inner_join(.,res$bulk.props %>%
               as.data.frame() %>%
               t() %>%
               as.data.frame() %>%
               rownames_to_column(var = "sample")) %>%
  mutate(Group = ifelse(`MARCO+ TAMs` > median(`MARCO+ TAMs`),"High MARCO+ TAMs","Low MARCO+ TAMs") %>%
           factor(.,levels=c("Low MARCO+ TAMs","High MARCO+ TAMs")),
         pathologic_stage = gsub("[A|C]","",pathologic_stage)) 

pheno[pheno==""] <- NA

Factors <- c("age","OS.time","gender", "pathologic_t", "pathologic_n", "pathologic_m",
             "pathologic_stage","Survival")
table1 <- CreateTableOne(vars = Factors,strata = "Group", data = pheno)
dat <- print(table1, showAllLevels = TRUE, formatOptions = list(big.mark = ","))
write.table(dat, "Figure/table.txt",sep = "\t",quote = F)

#########
###########
#################
datkm <- bioque_estimates

datkm$`MARCO+ TAMs` <- ifelse(datkm$`MARCO+ TAMs` > median(datkm$`MARCO+ TAMs`),"High MARCO+ TAMs","Low MARCO+ TAMs") %>%
  factor(.,levels=c("Low MARCO+ TAMs","High MARCO+ TAMs"))

sfit <- survfit(Surv(datkm$OS.time/365, datkm$OS) ~ datkm$`MARCO+ TAMs`,data = datkm)
p2_10 <-  ggsurvplot(sfit,palette = c("#1E90FF","#FF8C00"),
                     legend.title = element_blank(),
                     legend.labs = c("Low MARCO+ TAMs","High MARCO+ TAMs"), 
                     pval = T,risk.table = F,
                     legend = c(0.7,0.6),
                     legend.lab = c("Low MARCO+ TAMs","High MARCO+ TAMs"),
                     ggtheme = theme_bw(),
                     xlab = "Time in years",ylab = "Overall Survival Probility") 

p2_10 <- p2_10$plot + scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.background = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=12))
p2_10

save_plot("Figure/Figure2_7.pdf",p2_10,base_width = 4.5,base_height =3.2)


datkm$`SELENOP+ TAMs` <- ifelse(datkm$`SELENOP+ TAMs` > median(datkm$`SELENOP+ TAMs`),"High SELENOP+ TAMs","Low SELENOP+ TAMs") %>%
  factor(.,levels=c("Low SELENOP+ TAMs","High SELENOP+ TAMs"))
sfit2 <- survfit(Surv(datkm$OS.time/365, datkm$OS==1) ~ datkm$`SELENOP+ TAMs`,data = datkm)
p2_11 <-  ggsurvplot(sfit2,palette = c("#1E90FF","#FF8C00"),
                     legend.title = element_blank(),
                     legend.labs = c("Low SELENOP+ TAMs","High SELENOP+ TAMs"), 
                     pval = T,risk.table = F,
                     legend = c(0.7,0.6),
                     legend.lab = c("Low SELENOP+ TAMs","High SELENOP+ TAMs"),
                     ggtheme = theme_bw(),
                     xlab = "Time in years",ylab = "Overall Survival Probility") 

p2_11 <- p2_11$plot + scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.background = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=12))

p2_11
save_plot("Figure/Figure2_8.pdf",p2_11,base_width = 4.5,base_height =3.2)




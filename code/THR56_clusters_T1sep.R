

# Clean the working directory
rm(list = ls())

# Load necessary packages
require(Biobase)
require(limma)
require(pROC)
require(caret)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(plotROC)
library(mltools)
library(xtable)
library(dplyr)
library(precrec)
library(patchwork)
library(survminer)
library(survival)
library(tidyverse)
library(pheatmap)
library(glmnet)
library(readxl)

#################
# load the THR signatures
THR_signature <- readxl::read_xlsx("./data/THR Signatures_sep23.xlsx")

# get the THR50 signature
THR_50 <- THR_signature$`THR-50.1`[!is.na(THR_signature$`THR-50.1`)]

THR_50 <- gsub('-', '', THR_50)

################
# load the ET signatures
################
ET_signatures <- read_xlsx('objs/ET-9 Selection Steps.xlsx')
ET60 <- ET_signatures$`ET-60`[!is.na(ET_signatures$`ET-60`)]
ET9 <- ET_signatures$`ET-9`[!is.na(ET_signatures$`ET-9`)]

################
# fix gene names
################
ET60[ET60=='ADGRG1'] <- 'GPR56'
ET60[ET60=='DENND6B'] <- 'FAM116B'
ET60[ET60=='TENT5B'] <- 'FAM46B'
ET60[ET60=='HSA011916'] <- 'CTDNEP1'

ET9[ET9=='ADGRG1'] <- 'GPR56'
ET9[ET9=='DENND6B'] <- 'FAM116B'
ET9[ET9=='TENT5B'] <- 'FAM46B'
ET9[ET9=='HSA011916'] <- 'CTDNEP1'

################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')

# fix gene names
rownames(Expr_metabric_refAll)[grep('^ZNF652', rownames(Expr_metabric_refAll))]

# filter the THR signatures to include only the genes present in the expr matrices
THR_50_fil <- THR_50[THR_50 %in% rownames(Expr_metabric_refAll)]

THR_56 <- c(THR_50_fil, 'CDC20', 'LMNB2', 'KIF2C', 'FAM64A', 'KIF4A', 'TPX2')

# same for ET60 and ET9
ET60_fil <- ET60[ET60 %in% rownames(Expr_metabric_refAll)]
ET9_fil <- ET9[ET9 %in% rownames(Expr_metabric_refAll)]


#############################################################################################
#############################################################################################
## heatmap (THR 50)
Expr_metabric_refAll_heatmap <- Expr_metabric_refAll[THR_56, ] 

# Create annotation for columns/samples based on some clinical variables:
Pheno_metabric_forHeatmap <- Pheno_metabric
rownames(Pheno_metabric_forHeatmap) <- NULL

AnnAll_metabric <- Pheno_metabric_forHeatmap %>% 
  as.data.frame() %>%
  dplyr::select(Sample.ID, Pam50...Claudin.low.subtype, X3.Gene.classifier.subtype, HER2.Status, PR.Status, ER.status.measured.by.IHC, Neoplasm.Histologic.Grade) %>%
  column_to_rownames(var = "Sample.ID") %>%
  filter(Pam50...Claudin.low.subtype %in% c('Basal', 'claudin-low', 'Her2', 'LumA', 'LumB')) %>%
  dplyr::mutate(X3.Gene.classifier.subtype = as.factor(X3.Gene.classifier.subtype),
                ER.status.measured.by.IHC = as.factor(ER.status.measured.by.IHC),
                Pam50...Claudin.low.subtype = as.factor(Pam50...Claudin.low.subtype),
                HER2.Status = as.factor(HER2.Status), 
                PR.Status = as.factor(PR.Status),
                #Overall.Survival.Status = as.factor(Overall.Survival.Status),
                #Relapse.Free.Status = as.factor(Relapse.Free.Status), 
                #Tumor.Stage = as.factor(Tumor.Stage), 
                Neoplasm.Histologic.Grade = as.factor(Neoplasm.Histologic.Grade))


# filter and transpose the expression matrix
Expr_metabric_refAll_heatmap <- Expr_metabric_refAll_heatmap[, rownames(AnnAll_metabric)]
Expr_metabric_refAll_heatmap_t <- t(Expr_metabric_refAll_heatmap)

# filter pheno (above we remove normal and NC)
# save all pheno_metabric for future use
#Pheno_metabric2 <- Pheno_metabric
Pheno_metabric <- Pheno_metabric[rownames(AnnAll_metabric), ]

# colors
ann_colors = list()
ann_colors$Pam50...Claudin.low.subtype <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(5)
names(ann_colors$Pam50...Claudin.low.subtype) <- levels(AnnAll_metabric$Pam50...Claudin.low.subtype)

ann_colors$ER.status.measured.by.IHC <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(2)
names(ann_colors$ER.status.measured.by.IHC) <- levels(AnnAll_metabric$ER.status.measured.by.IHC)

ann_colors$X3.Gene.classifier.subtype <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(4)
names(ann_colors$X3.Gene.classifier.subtype) <- levels(AnnAll_metabric$X3.Gene.classifier.subtype)

ann_colors$Neoplasm.Histologic.Grade <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(3)
names(ann_colors$Neoplasm.Histologic.Grade) <- levels(AnnAll_metabric$Neoplasm.Histologic.Grade)


breaksList = seq(-4, 4, by = 1)
ColPal <- colorRampPalette(colors = rev(brewer.pal(11,"RdYlBu")))(20)
ColPal2 <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(20))


#######################################################
# get the 5 groups
heat_metabric <- pheatmap(Expr_metabric_refAll_heatmap, 
                          scale = "none",
                          #color = rev(heat.colors(20)),
                          color =ColPal,
                          annotation_colors = ann_colors,
                          cluster_cols = T, 
                          cluster_rows = T, 
                          clustering_distance_cols = 'correlation',
                          clustering_distance_rows = 'correlation',
                          clustering_method = 'ward.D',
                          show_colnames = F,
                          show_rownames = T,
                          annotation_col = AnnAll_metabric,
                          annotation_names_col = T,
                          #annotation_row = AnnAll_metabric,
                          annotation_names_row = T,
                          fontsize = 7,
                          #fontsize_col = 3,
                          fontsize_row = 10,
                          silent = TRUE,
                          cex = 1,
                          cutree_cols = 5,
                          cutree_rows = 5,
                          breaks = seq(-1, 1, by = 0.1),
                          main = "")

clusters_metabric <- as.data.frame(cbind(t(Expr_metabric_refAll_heatmap), 
                                         'THR clusters' = cutree(heat_metabric$tree_col, 
                                                                 k = 5)))

table(clusters_metabric$`THR clusters`)

# add the cluster info to the phenotype table
all(rownames(clusters_metabric) == rownames(Pheno_metabric))
Pheno_metabric$`THR clusters` <- clusters_metabric$`THR clusters`
#Pheno_metabric2$`THR clusters` <- clusters_metabric$`THR clusters`

# add the cluster info to the Ann dataframe and re-plot the heatmap
all(rownames(clusters_metabric) == rownames(AnnAll_metabric))
AnnAll_metabric$`THR clusters` <- as.factor(paste0('c', clusters_metabric$`THR clusters`))
table(AnnAll_metabric$`THR clusters`)

# re-order the annotation dataframe then the expression matrix by cluster
#AnnAll_metabric <- AnnAll_metabric[order(AnnAll_metabric$cluster, decreasing = FALSE), ]
#Expr_metabric_refAll_heatmap <- Expr_metabric_refAll_heatmap[, rownames(AnnAll_metabric)]


ann_colors$`THR clusters` <- c('#66A61E', "#7570B3", "#1B9E77" , '#D95F02', "#E7298A") 
levels(AnnAll_metabric$`THR clusters`) <- c('E1', 'T1', 'E3', 'E4', 'E2')
names(ann_colors$`THR clusters`) <- levels(AnnAll_metabric$`THR clusters`)

# CDC20 is duplicated, remove one
Expr_metabric_refAll_heatmap2 <- Expr_metabric_refAll_heatmap[!duplicated(rownames(Expr_metabric_refAll_heatmap)), ]

# merge E2 and E4
# Updating the 'THR clusters' column
clusters_metabric$`THR clusters` <- as.character(clusters_metabric$`THR clusters`)
clusters_metabric$`THR clusters`[clusters_metabric$`THR clusters` == "E4"] <- "E2"

AnnAll_metabric$`THR clusters` <- as.character(AnnAll_metabric$`THR clusters`)
AnnAll_metabric$`THR clusters`[AnnAll_metabric$`THR clusters` == "E4"] <- "E2"

# Updating the color scheme
ann_colors$`THR clusters` <- c('#66A61E', "#7570B3", "#1B9E77" , '#E7298A')
names(ann_colors$`THR clusters`) <- c('E1', 'T1', 'E3', 'E2')

# fix the cluster names in the pheno table
table(Pheno_metabric$`THR clusters`)
Pheno_metabric$`THR clusters` <- as.factor(Pheno_metabric$`THR clusters`)
#Pheno_metabric2$`THR clusters` <- as.factor(Pheno_metabric2$`THR clusters`)
levels(Pheno_metabric$`THR clusters`) <- c('E1', 'T1', 'E3', 'E2', 'E2')
#levels(Pheno_metabric2$`THR clusters`) <- c('E1', 'T1', 'E3', 'E2', 'E2')
table(Pheno_metabric$`THR clusters`)


#############################################################################################################
# get cluster 3 with crossing curves
T1_pheno <- Pheno_metabric[Pheno_metabric$`THR clusters` == 'T1', ]
T1_expr <- Expr_metabric_refAll[, rownames(T1_pheno)]

all(rownames(T1_pheno) == colnames(T1_expr))

######################################
## divide based on survival (RFS)
summary(T1_pheno$Relapse.Free.Status..Months.)

summary(T1_pheno$Relapse.Free.Status..Months. >= 50)

T1_pheno$T1_rfs_binary <- ifelse(T1_pheno$Relapse.Free.Status..Months. >= 50, 'longSurv', 'shortSurv')
table(T1_pheno$T1_rfs_binary)


#################################################################################################
## training
#################################################################################################

### combine in 1 dataset: Training
RFS_T1 <- as.factor(T1_pheno$T1_rfs_binary)
Data_T1 <- as.data.frame(cbind(t(T1_expr), RFS_T1))
Data_T1$RFS_T1 <- as.factor(Data_T1$RFS_T1)
levels(Data_T1$RFS_T1) <- c('longSurv', 'shortSurv')
table(Data_T1$RFS_T1)

################################################################
# separate T1 using THR50 derived I20 genes
################################################################
# load the THR50-derived i20 genes
THR50_i20 <- readxl::read_xlsx("./figures/c3_DE_THR50_RFS/THR50_c3_longVSshortSurv_DE.xlsx")$gene
THR70_i20 <- readxl::read_xlsx("./figures/T1_DE_THR70_RFS/THR70_T1_longVSshortSurv_DE.xlsx")$gene

# genes in common () THR50 I20 and THR70 I20 (i6)
I20common <- intersect(THR50_i20, THR70_i20)
I20common

##
################################################################
# THR50 i20 model
################################################################
i20model <- glm(as.formula((paste("RFS_T1 ~", paste(THR50_i20, collapse = "+")))), data = Data_T1, family = "binomial")
summary(i20model)
#save(i20model, file = './objs/THR50i20model_onTHR56clusters.rda')

################################################################
# i6 (common between THR50i20 and THR70i20) + HER2 model
# if i6+HER2- >> T1a, otherwise, T1b
################################################################
Data_T1_mod <- as.data.frame(cbind(t(T1_expr), T1_pheno))

colnames(T1_pheno)
table(T1_pheno$HER2.Status)

# Compute the mean expression of c3_gns genes for each sample
Data_T1_mod$EnrichmentScore_i6 <- rowMeans(Data_T1_mod[, I20common])

# Define the threshold as the median enrichment score
thr <- quantile(Data_T1_mod$EnrichmentScore_i6, 0.9)

# Create a new variable for subcluster assignment
Data_T1_mod$predClasses_i6_HER2 <- ifelse(Data_T1_mod$EnrichmentScore_i6 > thr & Data_T1_mod$HER2.Status == 'Negative', "T1a", "T1b")

table(Data_T1_mod$predClasses_i6_HER2)

################################################################
# THR50i20 + HER2 model
# if i20+HER2- >> T1a, otherwise, T1b
################################################################
# Compute the mean expression of c3_gns genes for each sample
Data_T1_mod$EnrichmentScore_THR50i20 <- rowMeans(Data_T1_mod[, THR50_i20])

# Define the threshold as the median enrichment score
thr_THR50i20 <- quantile(Data_T1_mod$EnrichmentScore_THR50i20, 0.9)

# Create a new variable for subcluster assignment
Data_T1_mod$predClasses_THR50i20_HER2 <- ifelse(Data_T1_mod$EnrichmentScore_THR50i20 > thr_THR50i20 & Data_T1_mod$HER2.Status == 'Negative', "T1a", "T1b")

table(Data_T1_mod$predClasses_THR50i20_HER2)

################################################################
# THR50i20 + HER2 model
# 4 classes: 
#i20 plus + HER2 minus
#i20 plus + HER2 plus
#i20 minus + HER2 plus
#i20 plus + HER2 minus
################################################################
# Compute the mean expression of c3_gns genes for each sample
Data_T1_mod$EnrichmentScore_THR50i20 <- rowMeans(Data_T1_mod[, THR50_i20])

# Define the threshold as the median enrichment score
thr_THR50i20 <- quantile(Data_T1_mod$EnrichmentScore_THR50i20, 0.9)

# Create a new variable for subcluster assignment
# Create a new variable for subcluster assignment
Data_T1_mod$predClasses_THR50i20_HER2_4classes <- ifelse(Data_T1_mod$EnrichmentScore_THR50i20 > thr_THR50i20 & Data_T1_mod$HER2.Status == 'Negative', "I20+HER2-",
                                                ifelse(Data_T1_mod$EnrichmentScore_THR50i20 > thr_THR50i20 & Data_T1_mod$HER2.Status == 'Positive', "I20+HER2+",
                                                       ifelse(Data_T1_mod$EnrichmentScore_THR50i20 <= thr_THR50i20 & Data_T1_mod$HER2.Status == 'Positive', "I20-HER2+", "I20-HER2-")))

table(Data_T1_mod$predClasses_THR50i20_HER2_4classes)

################################################################
# THR50-i20 + HER2 logreg model
################################################################

THR50i20_HER2_logreg_model <- glm(as.formula((paste("RFS_T1 ~", paste(c(THR50_i20, 'HER2.Status'), collapse = "+")))), data = Data_T1_mod, family = "binomial")
summary(THR50i20_HER2_logreg_model)
#save(THR50i20_HER2_logreg_model, file = './objs/THR50i20_HER2_logreg_model_onTHR56clusters.rda')

################################################################
# THR70 i20 model
################################################################
# fix HLA-DOB
THR70_i20[THR70_i20 == 'HLA-DOB'] <- 'HLA_DOB'
colnames(Data_T1)[colnames(Data_T1) == 'HLA-DOB'] <- 'HLA_DOB'


THR70i20model <- glm(as.formula((paste("RFS_T1 ~", paste(THR70_i20, collapse = "+")))), data = Data_T1, family = "binomial")
summary(THR70i20model)
#save(THR70i20model, file = './objs/THR70i20model_onTHR56clusters.rda')

################################################################
# ET60 model
################################################################
ET60model <- glm(as.formula((paste("RFS_T1 ~", paste(ET60_fil, collapse = "+")))), data = Data_T1, family = "binomial")
summary(ET60model)
#save(ET60model, file = './objs/ET60model_onTHR56clusters.rda')

################################################################
# ET9 model
################################################################
ET9model <- glm(as.formula((paste("RFS_T1 ~", paste(ET9_fil, collapse = "+")))), data = Data_T1, family = "binomial")
summary(ET9model)
#save(ET9model, file = './objs/ET9model_onTHR56clusters.rda')

################################################################
# CCP model
################################################################
CCP_genes <- c('FOXM1',	'ASPM',	'TK1',	'PRC1',
               'CDC20',	'BUB1B',	'PBK',	'DTL',
               'CDKN3',	'RRM2',	'ASF1B',	'CEP55',
               'CDC2',	'DLGAP5',	'C18orf24',	'RAD51',
               'KIF11',	'BIRC5',	'RAD54L',	'CENPM',
               'KIAA0101',	'KIF20A',	'PTTG1',	'CDCA8',
               'NUSAP1',	'PLK1',	'CDCA3',	'ORC6L',
               'CENPF',	'TOP2A',	'MCM10')

CCP_genes <- CCP_genes[CCP_genes %in% rownames(Expr_metabric_refAll)]

CCPmodel <- glm(as.formula((paste("RFS_T1 ~", paste(CCP_genes, collapse = "+")))), data = Data_T1, family = "binomial")
summary(CCPmodel)
#save(CCPmodel, file = './objs/CCPmodel_onTHR56clusters.rda')

############################################################################
# Make predictions

THR56_T1_prob_i20model <- i20model %>% predict(Data_T1 , type = "response")
THR56_T1_prob_THR50i20_HER2_logreg_model <- THR50i20_HER2_logreg_model %>% predict(Data_T1_mod , type = "response")
THR56_T1_prob_THR70i20model <- THR70i20model %>% predict(Data_T1 , type = "response")
THR56_T1_prob_ET60model <- ET60model %>% predict(Data_T1 , type = "response")
THR56_T1_prob_ET9model <- ET9model %>% predict(Data_T1 , type = "response")
THR56_T1_prob_CCPmodel <- CCPmodel %>% predict(Data_T1 , type = "response")

### Threshold
thr_THR56_T1_i20model <- coords(roc(RFS_T1, THR56_T1_prob_i20model, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR56_T1_THR50i20_HER2_logreg_model <- coords(roc(RFS_T1, THR56_T1_prob_THR50i20_HER2_logreg_model, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR56_T1_THR70i20model <- coords(roc(RFS_T1, THR56_T1_prob_THR70i20model, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR56_T1_ET60model <- coords(roc(RFS_T1, THR56_T1_prob_ET60model, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR56_T1_ET9model <- coords(roc(RFS_T1, THR56_T1_prob_ET9model, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR56_T1_CCPmodel <- coords(roc(RFS_T1, THR56_T1_prob_CCPmodel, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]


### ROC Curve
ROC_THR56_T1_i20model <- roc(RFS_T1, THR56_T1_prob_i20model, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_THR56_T1_i20model

ROC_THR56_T1_THR50i20_HER2_logreg_model <- roc(RFS_T1, THR56_T1_prob_THR50i20_HER2_logreg_model, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_THR56_T1_THR50i20_HER2_logreg_model

ROC_THR56_T1_THR70i20model <- roc(RFS_T1, THR56_T1_prob_THR70i20model, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_THR56_T1_THR70i20model

ROC_THR56_T1_ET60model <- roc(RFS_T1, THR56_T1_prob_ET60model, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_THR56_T1_ET60model

ROC_THR56_T1_ET9model <- roc(RFS_T1, THR56_T1_prob_ET9model, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_THR56_T1_ET9model

ROC_THR56_T1_CCPmodel <- roc(RFS_T1, THR56_T1_prob_CCPmodel, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_THR56_T1_CCPmodel

### Get predictions based on best threshold from ROC curve
predClasses_THR56_T1_i20model <- ifelse(THR56_T1_prob_i20model >= thr_THR56_T1_i20model$threshold, "longSurv", "shortSurv")
table(predClasses_THR56_T1_i20model)
predClasses_THR56_T1_i20model <- factor(predClasses_THR56_T1_i20model, levels = c('longSurv', 'shortSurv'))

### Get predictions based on best threshold from ROC curve
predClasses_THR56_T1_THR50i20_HER2_logreg_model <- ifelse(THR56_T1_prob_THR50i20_HER2_logreg_model >= thr_THR56_T1_THR50i20_HER2_logreg_model$threshold, "longSurv", "shortSurv")
table(predClasses_THR56_T1_THR50i20_HER2_logreg_model)
predClasses_THR56_T1_THR50i20_HER2_logreg_model <- factor(predClasses_THR56_T1_THR50i20_HER2_logreg_model, levels = c('longSurv', 'shortSurv'))

### Get predictions based on best threshold from ROC curve: THR70i20
predClasses_THR56_T1_THR70i20model <- ifelse(THR56_T1_prob_THR70i20model >= thr_THR56_T1_THR70i20model$threshold, "longSurv", "shortSurv")
table(predClasses_THR56_T1_THR70i20model)
predClasses_THR56_T1_THR70i20model <- factor(predClasses_THR56_T1_THR70i20model, levels = c('longSurv', 'shortSurv'))

### Get predictions based on best threshold from ROC curve: ET60
predClasses_THR56_T1_ET60model <- ifelse(THR56_T1_prob_ET60model >= thr_THR56_T1_ET60model$threshold, "longSurv", "shortSurv")
table(predClasses_THR56_T1_ET60model)
predClasses_THR56_T1_ET60model <- factor(predClasses_THR56_T1_ET60model, levels = c('longSurv', 'shortSurv'))

### Get predictions based on best threshold from ROC curve: ET9
predClasses_THR56_T1_ET9model <- ifelse(THR56_T1_prob_ET9model >= thr_THR56_T1_ET9model$threshold, "longSurv", "shortSurv")
table(predClasses_THR56_T1_ET9model)
predClasses_THR56_T1_ET9model <- factor(predClasses_THR56_T1_ET9model, levels = c('longSurv', 'shortSurv'))

### Get predictions based on best threshold from ROC curve: Ki67
predClasses_THR56_T1_CCPmodel <- ifelse(THR56_T1_prob_CCPmodel >= thr_THR56_T1_CCPmodel$threshold, "longSurv", "shortSurv")
table(predClasses_THR56_T1_CCPmodel)
predClasses_THR56_T1_CCPmodel <- factor(predClasses_THR56_T1_CCPmodel, levels = c('longSurv', 'shortSurv'))

######################
## Keep only the relevant information (Metastasis Event and Time)
all(rownames(Data_T1_mod) == rownames(T1_pheno))
predClasses_i6_HER2 <- Data_T1_mod$predClasses_i6_HER2
predClasses_THR50i20_HER2 <- Data_T1_mod$predClasses_THR50i20_HER2
predClasses_THR50i20_HER2_4classes <- Data_T1_mod$predClasses_THR50i20_HER2_4classes

T1_pheno2 <- cbind(T1_pheno[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Relapse.Free.Status", "Relapse.Free.Status..Months.", "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC", "X3.Gene.classifier.subtype")], 
                        THR56_T1_prob_i20model, predClasses_THR56_T1_i20model,
                        THR56_T1_prob_THR50i20_HER2_logreg_model, predClasses_THR56_T1_THR50i20_HER2_logreg_model,
                        THR56_T1_prob_THR70i20model, predClasses_THR56_T1_THR70i20model,
                        THR56_T1_prob_ET60model, predClasses_THR56_T1_ET60model,
                        THR56_T1_prob_ET9model, predClasses_THR56_T1_ET9model,
                        THR56_T1_prob_CCPmodel, predClasses_THR56_T1_CCPmodel,
                        predClasses_THR50i20_HER2, predClasses_i6_HER2,
                        predClasses_THR50i20_HER2_4classes
)


CoxData_metabric_T1 <- data.frame(T1_pheno2)

##########################################################################################
##########################################################################################
##########################################################################################
## survival analysis for just T1

#########################################
# using the THR50 derived i20 model
#########################################

# OS
Fit_sig_metabric_os_THR56_T1_i20model <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR56_T1_i20model, data = CoxData_metabric_T1)

# RFS
Fit_sig_metabric_RFS_THR56_T1_i20model <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR56_T1_i20model, data = CoxData_metabric_T1)


# plot OS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_os_T1_i20model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR56_T1_i20model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR56: using the THR50-derived i20 model'
)
dev.off()

######################################
# plot RFS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_rfs_T1_i20model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR56_T1_i20model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR56: using the THR50-derived i20 model'
)
dev.off()

#########################################
# using the THR50-i20 + HER2 as log reg model
#########################################

# OS
Fit_sig_metabric_os_THR56_T1_THR50i20_HER2_logreg_model <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR56_T1_THR50i20_HER2_logreg_model, data = CoxData_metabric_T1)

# RFS
Fit_sig_metabric_RFS_THR56_T1_THR50i20_HER2_logreg_model <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR56_T1_THR50i20_HER2_logreg_model, data = CoxData_metabric_T1)


# plot OS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_os_T1_THR50i20_HER2_logreg_model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR56_T1_THR50i20_HER2_logreg_model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR56: using the THR50-i20 + HER2 logreg model'
)
dev.off()

##################
# plot RFS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_rfs_T1_THR50i20_HER2_logreg_model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR56_T1_THR50i20_HER2_logreg_model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR56: using the THR50-i20 + HER2 logreg model'
)
dev.off()

#########################################
# using i6 + HER2
#########################################

# OS
Fit_sig_metabric_os_THR56_T1_i6_HER2_model <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_i6_HER2, data = CoxData_metabric_T1)

# RFS
Fit_sig_metabric_RFS_THR56_T1_i6_HER2_model <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_i6_HER2, data = CoxData_metabric_T1)


# plot OS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_os_T1_i6_HER2_model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR56_T1_i6_HER2_model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR56: i6 + HER2'
)
dev.off()

###########################
# plot RFS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_rfs_T1_i6_HER2_model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR56_T1_i6_HER2_model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR56: i6 + HER2'
)
dev.off()

#########################################
# using the THR50 derived i20 + HER2
#########################################

# OS
Fit_sig_metabric_os_THR56_T1_THR50i20_HER2_model <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR50i20_HER2, data = CoxData_metabric_T1)

# RFS
Fit_sig_metabric_RFS_THR56_T1_THR50i20_HER2_model <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR50i20_HER2, data = CoxData_metabric_T1)


# plot OS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_os_T1_THR50i20_HER2_model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR56_T1_THR50i20_HER2_model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR56: THR50-i20 + HER2'
)
dev.off()

###########################
# plot RFS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_rfs_T1_THR50i20_HER2_model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR56_T1_THR50i20_HER2_model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR56: THR50-i20 + HER2'
)
dev.off()

#########################################
# using the THR50 derived i20 + HER2: 4 classes
#########################################

classes <- levels(factor(CoxData_metabric_T1$predClasses_THR50i20_HER2_4classes))

# OS
Fit_sig_metabric_os_THR56_T1_THR50i20_HER2_4classes_model <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR50i20_HER2_4classes, data = CoxData_metabric_T1)

# RFS
Fit_sig_metabric_RFS_THR56_T1_THR50i20_HER2_4classes_model <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR50i20_HER2_4classes, data = CoxData_metabric_T1)


# plot OS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_os_T1_THR50i20_HER2_4classes_model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR56_T1_THR50i20_HER2_4classes_model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = classes,
           ggtheme = theme_survminer(base_size = 20, font.x = c(20, 'bold.italic', 'black'), font.y = c(20, 'bold.italic', 'black'), font.tickslab = c(20, 'plain', 'black'), font.legend = c(20, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR56: THR50-i20 + HER2'
)
dev.off()

###########################
# plot RFS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_rfs_T1_THR50i20_HER2_4classes_model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR56_T1_THR50i20_HER2_4classes_model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = classes,
           ggtheme = theme_survminer(base_size = 20, font.x = c(20, 'bold.italic', 'black'), font.y = c(20, 'bold.italic', 'black'), font.tickslab = c(20, 'plain', 'black'), font.legend = c(20, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR56: THR50-i20 + HER2'
)
dev.off()

#########################################
# using the THR70 derived i20 model
#########################################

# OS
Fit_sig_metabric_os_THR56_T1_THR70i20model <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR56_T1_THR70i20model, data = CoxData_metabric_T1)

# RFS
Fit_sig_metabric_RFS_THR56_T1_THR70i20model <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR56_T1_THR70i20model, data = CoxData_metabric_T1)


# plot OS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_os_T1_THR70i20model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR56_T1_THR70i20model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR56: using the THR70-derived i20 model'
)
dev.off()

######################################
# plot RFS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_rfs_T1_THR70i20model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR56_T1_THR70i20model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR56: using the THR70-derived i20 model'
)
dev.off()

#########################################
# using ET60 model
#########################################

# OS
Fit_sig_metabric_os_THR56_T1_ET60model <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR56_T1_ET60model, data = CoxData_metabric_T1)

# RFS
Fit_sig_metabric_RFS_THR56_T1_ET60model <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR56_T1_ET60model, data = CoxData_metabric_T1)


# plot OS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_os_T1_ET60model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR56_T1_ET60model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR56: using ET60'
)
dev.off()

######################################
# plot RFS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_rfs_T1_ET60model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR56_T1_ET60model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR56: using ET60'
)
dev.off()

#########################################
# using ET9 model
#########################################

# OS
Fit_sig_metabric_os_THR56_T1_ET9model <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR56_T1_ET9model, data = CoxData_metabric_T1)

# RFS
Fit_sig_metabric_RFS_THR56_T1_ET9model <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR56_T1_ET9model, data = CoxData_metabric_T1)


# plot OS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_os_T1_ET9model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR56_T1_ET9model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR56: using ET9'
)
dev.off()

######################################
# plot RFS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_rfs_T1_ET9model.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR56_T1_ET9model,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR56: using ET9'
)
dev.off()

#########################################
# using CCP model
#########################################

# OS
Fit_sig_metabric_os_THR56_T1_CCPmodel <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR56_T1_CCPmodel, data = CoxData_metabric_T1)

# RFS
Fit_sig_metabric_RFS_THR56_T1_CCPmodel <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR56_T1_CCPmodel, data = CoxData_metabric_T1)


# plot OS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_os_T1_CCPmodel.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR56_T1_CCPmodel,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR56: using CCP'
)
dev.off()

######################################
# plot RFS
tiff("./figures/logreg/THR56_clusters/THR56_metabric_rfs_T1_CCPmodel.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR56_T1_CCPmodel,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR56: using CCP'
)
dev.off()

##########################################################################################
##########################################################################################
## recombine c3 with the rest 
T1_pheno2$THR_clusters_i20model <- as.factor(T1_pheno2$predClasses_THR56_T1_i20model)
levels(T1_pheno2$THR_clusters_i20model) <- c('T1a', 'T1b')

T1_pheno2$THR_clusters_i6_HER2 <- as.factor(T1_pheno2$predClasses_i6_HER2)
levels(T1_pheno2$THR_clusters_i6_HER2) <- c('T1a', 'T1b')

T1_pheno2$THR_clusters_THR50i20_HER2 <- as.factor(T1_pheno2$predClasses_THR50i20_HER2)
levels(T1_pheno2$THR_clusters_THR50i20_HER2) <- c('T1a', 'T1b')

T1_pheno2$THR_clusters_THR50i20_HER2_4classes <- as.factor(T1_pheno2$predClasses_THR50i20_HER2_4classes)
levels(T1_pheno2$THR_clusters_THR50i20_HER2_4classes) <- c("I20-HER2-", "I20-HER2+", "I20+HER2-", "I20+HER2+")

T1_pheno2$THR_clusters_THR70i20model <- as.factor(T1_pheno2$predClasses_THR56_T1_THR70i20model)
levels(T1_pheno2$THR_clusters_THR70i20model) <- c('T1a', 'T1b')

T1_pheno2$THR_clusters_ET60model <- as.factor(T1_pheno2$predClasses_THR56_T1_ET60model)
levels(T1_pheno2$THR_clusters_ET60model) <- c('T1a', 'T1b')

T1_pheno2$THR_clusters_ET9model <- as.factor(T1_pheno2$predClasses_THR56_T1_ET9model)
levels(T1_pheno2$THR_clusters_ET9model) <- c('T1a', 'T1b')

T1_pheno2$THR_clusters_CCPmodel <- as.factor(T1_pheno2$predClasses_THR56_T1_CCPmodel)
levels(T1_pheno2$THR_clusters_CCPmodel) <- c('T1a', 'T1b')


T1 <- data.frame(THR_clusters_i20model = T1_pheno2$THR_clusters_i20model,
                 THR_clusters_THR70i20model = T1_pheno2$THR_clusters_THR70i20model,
                 THR_clusters_ET60model = T1_pheno2$THR_clusters_ET60model,
                 THR_clusters_ET9model = T1_pheno2$THR_clusters_ET9model,
                 THR_clusters_CCPmodel = T1_pheno2$THR_clusters_CCPmodel,
                 THR_clusters_i6_HER2 = T1_pheno2$THR_clusters_i6_HER2,
                 THR_clusters_THR50i20_HER2 = T1_pheno2$THR_clusters_THR50i20_HER2,
                 THR_clusters_THR50i20_HER2_4classes = T1_pheno2$THR_clusters_THR50i20_HER2_4classes,
                 `Sample.ID` = rownames(T1_pheno))

rownames(T1) <- rownames(T1_pheno)


# merge
Pheno_metabric$`THR clusters`[Pheno_metabric$`THR clusters` == 'T1'] <- NA


Pheno_metabric2 <- merge(x = T1, y = Pheno_metabric, by="Sample.ID", all.y = TRUE)

Pheno_metabric2 <- Pheno_metabric2 %>% 
  mutate(`THR clusters` = as.factor(`THR clusters`), THR_clusters_i20model = as.factor(THR_clusters_i20model), THR_clusters_THR70i20model = as.factor(THR_clusters_THR70i20model), THR_clusters_i6_HER2 = as.factor(THR_clusters_i6_HER2), THR_clusters_THR50i20_HER2 = as.factor(THR_clusters_THR50i20_HER2), THR_clusters_THR50i20_HER2_4classes = as.factor(THR_clusters_THR50i20_HER2_4classes), THR_clusters_ET60model = as.factor(THR_clusters_ET60model), THR_clusters_ET9model = as.factor(THR_clusters_ET9model), THR_clusters_CCPmodel = as.factor(THR_clusters_CCPmodel)) %>%
  mutate(THR.clusters_i20model_Merged = coalesce(THR_clusters_i20model,`THR clusters`)) %>%
  mutate(THR.clusters_THR70i20model_Merged = coalesce(THR_clusters_THR70i20model,`THR clusters`)) %>%
  mutate(THR.clusters_ET60model_Merged = coalesce(THR_clusters_ET60model,`THR clusters`)) %>%
  mutate(THR.clusters_ET9model_Merged = coalesce(THR_clusters_ET9model,`THR clusters`)) %>%
  mutate(THR.clusters_CCPmodel_Merged = coalesce(THR_clusters_CCPmodel,`THR clusters`)) %>%
  mutate(THR.clusters_i6_HER2_Merged = coalesce(THR_clusters_i6_HER2,`THR clusters`)) %>%
  mutate(THR.clusters_THR50i20_HER2_Merged = coalesce(THR_clusters_THR50i20_HER2,`THR clusters`)) %>%
  mutate(THR.clusters_THR50i20_HER2_4classes_Merged = coalesce(THR_clusters_THR50i20_HER2_4classes,`THR clusters`))

table(Pheno_metabric2$THR.clusters_i20model_Merged, Pheno_metabric2$THR.clusters_THR70i20model_Merged)
table(Pheno_metabric2$THR.clusters_i20model_Merged, Pheno_metabric2$THR.clusters_i6_HER2_Merged)
table(Pheno_metabric2$THR.clusters_THR50i20_HER2_Merged, Pheno_metabric2$THR.clusters_i6_HER2_Merged)

###########################################################################################
##########################################################################################
## survival analysis

## Keep only the relevant information (Metastasis Event and Time)
survival_metabric <- Pheno_metabric2[, c("Overall.Survival.Status", "Overall.Survival..Months.", 
                                         "Relapse.Free.Status", "Relapse.Free.Status..Months.", 
                                         "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC",
                                         "X3.Gene.classifier.subtype", 
                                         "THR.clusters_i20model_Merged",
                                         "THR.clusters_THR70i20model_Merged",
                                         "THR.clusters_ET60model_Merged",
                                         "THR.clusters_ET9model_Merged",
                                         "THR.clusters_CCPmodel_Merged",
                                         "THR.clusters_i6_HER2_Merged",
                                         "THR.clusters_THR50i20_HER2_Merged",
                                         "THR.clusters_THR50i20_HER2_4classes_Merged"
)] 

survival_metabric$THR.clusters_i20model_Merged <- as.factor(survival_metabric$THR.clusters_i20model_Merged)
survival_metabric$THR.clusters_i20model_Merged <- droplevels(survival_metabric$THR.clusters_i20model_Merged)

survival_metabric$THR.clusters_i6_HER2_Merged <- as.factor(survival_metabric$THR.clusters_i6_HER2_Merged)
survival_metabric$THR.clusters_i6_HER2_Merged <- droplevels(survival_metabric$THR.clusters_i6_HER2_Merged)

survival_metabric$THR.clusters_THR50i20_HER2_Merged <- as.factor(survival_metabric$THR.clusters_THR50i20_HER2_Merged)
survival_metabric$THR.clusters_THR50i20_HER2_Merged <- droplevels(survival_metabric$THR.clusters_THR50i20_HER2_Merged)

survival_metabric$THR.clusters_THR50i20_HER2_4classes_Merged <- as.factor(survival_metabric$THR.clusters_THR50i20_HER2_4classes_Merged)
survival_metabric$THR.clusters_THR50i20_HER2_4classes_Merged <- droplevels(survival_metabric$THR.clusters_THR50i20_HER2_4classes_Merged)

survival_metabric$THR.clusters_THR70i20model_Merged <- as.factor(survival_metabric$THR.clusters_THR70i20model_Merged)
survival_metabric$THR.clusters_THR70i20model_Merged <- droplevels(survival_metabric$THR.clusters_THR70i20model_Merged)

survival_metabric$THR.clusters_ET60model_Merged <- as.factor(survival_metabric$THR.clusters_ET60model_Merged)
survival_metabric$THR.clusters_ET60model_Merged <- droplevels(survival_metabric$THR.clusters_ET60model_Merged)

survival_metabric$THR.clusters_ET9model_Merged <- as.factor(survival_metabric$THR.clusters_ET9model_Merged)
survival_metabric$THR.clusters_ET9model_Merged <- droplevels(survival_metabric$THR.clusters_ET9model_Merged)

survival_metabric$THR.clusters_CCPmodel_Merged <- as.factor(survival_metabric$THR.clusters_CCPmodel_Merged)
survival_metabric$THR.clusters_CCPmodel_Merged <- droplevels(survival_metabric$THR.clusters_CCPmodel_Merged)

cluster_colors <- as.vector(ann_colors$THR.clusters_i20model_Merged) # same order for the others

################################################################
## Survival curves: THR50 i20
################################################################

# OS
Fit_metabric_os_i20model <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_i20model_Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_i20model <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_i20model_Merged, data = survival_metabric)

pdf("./figures/logreg/THR56_clusters/metabric_os_5clusters_model20Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_i20model,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_i20model_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and OS: THR56 + THR50-derived i20')
dev.off()

## RFS: 
pdf("./figures/logreg/THR56_clusters/metabric_rfs_5clusters_model20Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_i20model,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_i20model_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and RFS: THR56 + THR50-derived i20')
dev.off()

################################################################
## Survival curves: i6 + HER2
################################################################

# OS
Fit_metabric_os_i6_HER2 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_i6_HER2_Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_i6_HER2 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_i6_HER2_Merged, data = survival_metabric)

pdf("./figures/logreg/THR56_clusters/metabric_os_5clusters_i6_HER2_Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_i6_HER2,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_i6_HER2_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and OS: THR56 + i6 + HER2')
dev.off()

## RFS: 
pdf("./figures/logreg/THR56_clusters/metabric_rfs_5clusters_i6_HER2_Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_i6_HER2,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_i6_HER2_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and RFS: THR56 + i6 + HER2')
dev.off()

################################################################
## Survival curves: THR50-i20 + HER2
################################################################

# OS
Fit_metabric_os_THR50i20_HER2 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_THR50i20_HER2_Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_THR50i20_HER2 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_THR50i20_HER2_Merged, data = survival_metabric)

pdf("./figures/logreg/THR56_clusters/metabric_os_5clusters_THR50i20_HER2_Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_THR50i20_HER2,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_THR50i20_HER2_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and OS: THR56 + THR50-i20 + HER2')
dev.off()

## RFS: 
pdf("./figures/logreg/THR56_clusters/metabric_RFS_5clusters_THR50i20_HER2_Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_THR50i20_HER2,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_THR50i20_HER2_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and RFS: THR56 + THR50-i20 + HER2')
dev.off()

################################################################
## Survival curves: THR50-i20 + HER2: 4 classes
################################################################

# OS
Fit_metabric_os_THR50i20_HER2_4classes <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_THR50i20_HER2_4classes_Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_THR50i20_HER2_4classes <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_THR50i20_HER2_4classes_Merged, data = survival_metabric)

pdf("./figures/logreg/THR56_clusters/metabric_os_5clusters_THR50i20_HER2_4classes_Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_THR50i20_HER2_4classes,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_THR50i20_HER2_4classes_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and OS: THR56 + THR50-i20 + HER2')
dev.off()

## RFS: 
pdf("./figures/logreg/THR56_clusters/metabric_RFS_5clusters_THR50i20_HER2_4classes_Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_THR50i20_HER2_4classes,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_THR50i20_HER2_4classes_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and RFS: THR56 + THR50-i20 + HER2')
dev.off()

################################################################
## Survival curves: THR70 i20
################################################################

# OS
Fit_metabric_os_THR70i20model <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_THR70i20model_Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_THR70i20model <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_THR70i20model_Merged, data = survival_metabric)

pdf("./figures/logreg/THR56_clusters/metabric_os_5clusters_THR70i20Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_THR70i20model,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_i20model_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and OS: THR56 + THR70-derived i20')
dev.off()

## RFS: 
pdf("./figures/logreg/THR56_clusters/metabric_rfs_5clusters_THR70i20Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_THR70i20model,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_i20model_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and RFS: THR56 + THR70-derived i20')
dev.off()

################################################################
## Survival curves: ET60
################################################################

# OS
Fit_metabric_os_ET60model <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_ET60model_Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_ET60model <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_ET60model_Merged, data = survival_metabric)

pdf("./figures/logreg/THR56_clusters/metabric_os_5clusters_ET60Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_ET60model,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_ET60model_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and OS: THR56 + ET60')
dev.off()

## RFS: 
pdf("./figures/logreg/THR56_clusters/metabric_rfs_5clusters_ET60Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_ET60model,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_ET60model_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and RFS: THR56 + ET60')
dev.off()


################################################################
## Survival curves: ET9
################################################################

# OS
Fit_metabric_os_ET9model <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_ET9model_Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_ET9model <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_ET9model_Merged, data = survival_metabric)

pdf("./figures/logreg/THR56_clusters/metabric_os_5clusters_ET9Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_ET9model,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_ET9model_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and OS: THR56 + ET9')
dev.off()

## RFS: 
pdf("./figures/logreg/THR56_clusters/metabric_rfs_5clusters_ET9Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_ET9model,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_ET9model_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and RFS: THR56 + ET9')
dev.off()

################################################################
## Survival curves: CCP
################################################################

# OS
Fit_metabric_os_CCPmodel <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_CCPmodel_Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_CCPmodel <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_CCPmodel_Merged, data = survival_metabric)

pdf("./figures/logreg/THR56_clusters/metabric_os_5clusters_CCPMerged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_CCPmodel,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_CCPmodel_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and OS: THR56 + CCP')
dev.off()

## RFS: 
pdf("./figures/logreg/THR56_clusters/metabric_rfs_5clusters_CCPMerged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_CCPmodel,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_CCPmodel_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and RFS: THR56 + CCP')
dev.off()


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
library(vtable)

#################
THR_signature <- readxl::read_xlsx("./data/THR_Signatures_Jan25_2023.xlsx")

# get the THR70 signature
THR_70 <- THR_signature$`THR-70`[!is.na(THR_signature$`THR-70`)]

THR_70 <- gsub('-', '', THR_70)

################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')

# fix gene names
rownames(Expr_metabric_refAll)[grep('^ZNF652', rownames(Expr_metabric_refAll))]

# filter the THR signatures to include only the genes present in the expr matrices
THR_70_fil <- THR_70[THR_70 %in% rownames(Expr_metabric_refAll)]

#############################################################################################
#############################################################################################
## heatmap (THR 50)
Expr_metabric_refAll_heatmap <- Expr_metabric_refAll[THR_70_fil, ] 

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
                          silent = F,
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

# add the cluster info to the Ann dataframe and re-plot the heatmap
all(rownames(clusters_metabric) == rownames(AnnAll_metabric))
AnnAll_metabric$`THR clusters` <- as.factor(paste0('c', clusters_metabric$`THR clusters`))
table(AnnAll_metabric$`THR clusters`)

# re-order the annotation dataframe then the expression matrix by cluster
#AnnAll_metabric <- AnnAll_metabric[order(AnnAll_metabric$cluster, decreasing = FALSE), ]
#Expr_metabric_refAll_heatmap <- Expr_metabric_refAll_heatmap[, rownames(AnnAll_metabric)]


ann_colors$`THR clusters` <- colorRampPalette(colors = rev(brewer.pal(5,"Dark2")))(5)
levels(AnnAll_metabric$`THR clusters`) <- c('E2b', 'E2a', 'E1', 'E3', 'PNBC')
names(ann_colors$`THR clusters`) <- levels(AnnAll_metabric$`THR clusters`)

# fix the cluster names in the pheno table
table(Pheno_metabric$`THR clusters`)
Pheno_metabric$`THR clusters` <- as.factor(Pheno_metabric$`THR clusters`)
levels(Pheno_metabric$`THR clusters`) <- c('E2b', 'E2a', 'E1', 'E3', 'PNBC')
table(Pheno_metabric$`THR clusters`)

#############################################################################################################
# get cluster 3 with crossing curves
cluster3_pheno <- Pheno_metabric[Pheno_metabric$`THR clusters` == 'PNBC', ]
cluster3_expr <- Expr_metabric_refAll[, rownames(cluster3_pheno)]

all(rownames(cluster3_pheno) == colnames(cluster3_expr))

######################################
## divide based on survival (RFS)
summary(cluster3_pheno$Relapse.Free.Status..Months.)

summary(cluster3_pheno$Relapse.Free.Status..Months. >= 50)

cluster3_pheno$c3_rfs_binary <- ifelse(cluster3_pheno$Relapse.Free.Status..Months. >= 50, 'longSurv', 'shortSurv')
table(cluster3_pheno$c3_rfs_binary)

######################################
# differential expression 
design <- model.matrix( ~ cluster3_pheno$c3_rfs_binary)
colnames(design)[2] <- "longVSshortSurvival"

fit <- lmFit(cluster3_expr, design)

fitted.ebayes <- eBayes(fit)

c3_top20 <- topTable(fitted.ebayes, number = 20)
c3_top20$gene <- rownames(c3_top20)

c3_top200 <- topTable(fitted.ebayes, number = 200)
c3_top200$gene <- rownames(c3_top200)

#write_csv(as.data.frame(c3_top20), file = './figures/c3_DE_THR70_RFS/c3_longVSshortSurv_DE.csv')
library("writexl")
write_xlsx(c3_top20,"./figures/T1_DE_THR70_RFS/THR70_T1_longVSshortSurv_DE.xlsx")

# save top200 DE genes
write_xlsx(c3_top200,"./figures/T1_DE_THR70_RFS/THR70_T1_longVSshortSurv_DE_top200.xlsx")

#########
# genes in common () THR50 I20 and THR70 I20
THR50_c3_top20 <- readxl::read_xlsx("./figures/c3_DE_THR50_RFS/THR50_c3_longVSshortSurv_DE.xlsx")

I20common <- intersect(THR50_c3_top20$gene, c3_top20$gene)

############

c3_gns <- rownames(topTable(fitted.ebayes, number = 20))

# genes in common with THR70
summary(c3_gns %in% THR_70)

#summary(decideTests(fitted.ebayes[,"longVSshortSurvival"],lfc=0))



####
# c3 summary
# sumtable(Pheno_metabric,
#          group = 'THR clusters',
#          file='metabric_clusters_summary',
#          out = 'browser',
#          title='METABRIC clusters Summary Statistics',
#          simple.kable=FALSE,
#          opts=list())

#################################################################################################
## training
#################################################################################################

### combine in 1 dataset: Training
RFS_c3 <- as.factor(cluster3_pheno$c3_rfs_binary)
Data_c3 <- as.data.frame(cbind(t(cluster3_expr), RFS_c3))
Data_c3$RFS_c3 <- as.factor(Data_c3$RFS_c3)
levels(Data_c3$RFS_c3) <- c('longSurv', 'shortSurv')
table(Data_c3$RFS_c3)

# fix HLA-DOB
c3_gns[c3_gns == 'HLA-DOB'] <- 'HLA_DOB'
colnames(Data_c3)[colnames(Data_c3) == 'HLA-DOB'] <- 'HLA_DOB'

####################################
# THR70 I20 model
#####################################
THR70_I20_c3_model <- glm(as.formula((paste("RFS_c3 ~", paste(c3_gns, collapse = "+")))), data = Data_c3, family = "binomial")
summary(THR70_I20_c3_model)
save(THR70_I20_c3_model, file = './objs/THR70_I20_model.rda')

####################################
# THR70 I20: Recursive Feature Elimination 
#####################################
library(caret)
ctrl <- rfeControl(functions=rfFuncs, method="cv", number=10)
#results <- rfe(Data_c3[,THR50_c3_top20$gene], Data_c3$RFS_c3, sizes=c(1:4), rfeControl=ctrl)
#print(results)

# Get the optimal number of variables
# Select the top 5 variables based on the ranking
top_vars <- c('KIR3DL3', 'TIMP1', 'CD511953', 'NCBP2', 'SPRED2', 'STX1A', 'HIST1H2BF', 'USP30')
# Build the model with the optimal variables only
THR70_I20_c3_model_rfe <- glm(as.formula(paste("RFS_c3 ~", paste(top_vars, collapse = "+"))), data = Data_c3, family = "binomial")
summary(THR70_I20_c3_model_rfe)

####################################
# THR50 I20 model
#####################################
THR50_I20_c3_model <- glm(as.formula((paste("RFS_c3 ~", paste(THR50_c3_top20$gene, collapse = "+")))), data = Data_c3, family = "binomial")
summary(THR50_I20_c3_model)

####################################
# CTLA4 model
#####################################
CTLA4_model_c3 <- glm(RFS_c3 ~ CTLA4, data = Data_c3, family = "binomial")
summary(CTLA4_model_c3)

####################################
# CTLA4_ICOS_CXCL13_model
#####################################
CTLA4_ICOS_CXCL13_model_c3 <- glm(RFS_c3 ~ CTLA4+ICOS+CXCL13, data = Data_c3, family = "binomial")
summary(CTLA4_ICOS_CXCL13_model_c3)

####################################
# KLF7 model
#####################################
KLF7_model_c3 <- glm(RFS_c3 ~ KLF7, data = Data_c3, family = "binomial")
summary(KLF7_model_c3)

####################################
# ZNF627 model
#####################################
ZNF627_model_c3 <- glm(RFS_c3 ~ ZNF627, data = Data_c3, family = "binomial")
summary(ZNF627_model_c3)


############################################################################
# Make predictions

Train_prob_THR70_c3_THR70_I20 <- THR70_I20_c3_model %>% predict(Data_c3 , type = "response")
Train_prob_THR70_c3_THR70_I20_rfe <- THR70_I20_c3_model_rfe %>% predict(Data_c3 , type = "response")
Train_prob_THR70_c3_THR50_I20 <- THR50_I20_c3_model %>% predict(Data_c3 , type = "response")
Train_prob_THR70_c3_CTLA4 <- CTLA4_model_c3 %>% predict(Data_c3 , type = "response")
Train_prob_THR70_c3_CTLA4_ICOS_CXCL13 <- CTLA4_ICOS_CXCL13_model_c3 %>% predict(Data_c3 , type = "response")
Train_prob_THR70_c3_KLF7 <- KLF7_model_c3 %>% predict(Data_c3 , type = "response")
Train_prob_THR70_c3_ZNF627 <- ZNF627_model_c3 %>% predict(Data_c3 , type = "response")

### Threshold
thr_THR70_c3_THR70_I20 <- coords(roc(RFS_c3, Train_prob_THR70_c3_THR70_I20, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR70_c3_THR70_I20_rfe <- coords(roc(RFS_c3, Train_prob_THR70_c3_THR70_I20_rfe, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR70_c3_THR50_I20 <- coords(roc(RFS_c3, Train_prob_THR70_c3_THR50_I20, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR70_c3_CTLA4 <- coords(roc(RFS_c3, Train_prob_THR70_c3_CTLA4, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR70_c3_CTLA4_ICOS_CXCL13 <- coords(roc(RFS_c3, Train_prob_THR70_c3_CTLA4_ICOS_CXCL13, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR70_c3_KLF7 <- coords(roc(RFS_c3, Train_prob_THR70_c3_KLF7, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR70_c3_ZNF627 <- coords(roc(RFS_c3, Train_prob_THR70_c3_ZNF627, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]


### ROC Curve
ROCTrain_THR70_c3_THR70_I20 <- roc(RFS_c3, Train_prob_THR70_c3_THR70_I20, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR70_c3_THR70_I20

ROCTrain_THR70_c3_THR70_I20_rfe <- roc(RFS_c3, Train_prob_THR70_c3_THR70_I20_rfe, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR70_c3_THR70_I20_rfe

ROCTrain_THR70_c3_THR50_I20 <- roc(RFS_c3, Train_prob_THR70_c3_THR50_I20, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR70_c3_THR50_I20

ROCTrain_THR70_c3_CTLA4 <- roc(RFS_c3, Train_prob_THR70_c3_CTLA4, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR70_c3_CTLA4

ROCTrain_THR70_c3_CTLA4_ICOS_CXCL13 <- roc(RFS_c3, Train_prob_THR70_c3_CTLA4_ICOS_CXCL13, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR70_c3_CTLA4_ICOS_CXCL13

ROCTrain_THR70_c3_KLF7 <- roc(RFS_c3, Train_prob_THR70_c3_KLF7, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR70_c3_KLF7

ROCTrain_THR70_c3_ZNF627 <- roc(RFS_c3, Train_prob_THR70_c3_ZNF627, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR70_c3_ZNF627

### Get predictions based on best threshold from ROC curve
predClasses_THR70_c3_THR70_I20 <- ifelse(Train_prob_THR70_c3_THR70_I20 >= thr_THR70_c3_THR70_I20$threshold, "longSurv", "shortSurv")
table(predClasses_THR70_c3_THR70_I20)
predClasses_THR70_c3_THR70_I20 <- factor(predClasses_THR70_c3_THR70_I20, levels = c('longSurv', 'shortSurv'))

predClasses_THR70_c3_THR70_I20_rfe <- ifelse(Train_prob_THR70_c3_THR70_I20_rfe >= thr_THR70_c3_THR70_I20_rfe$threshold, "longSurv", "shortSurv")
table(predClasses_THR70_c3_THR70_I20_rfe)
predClasses_THR70_c3_THR70_I20_rfe <- factor(predClasses_THR70_c3_THR70_I20_rfe, levels = c('longSurv', 'shortSurv'))

predClasses_THR70_c3_THR50_I20 <- ifelse(Train_prob_THR70_c3_THR50_I20 >= thr_THR70_c3_THR50_I20$threshold, "longSurv", "shortSurv")
table(predClasses_THR70_c3_THR50_I20)
predClasses_THR70_c3_THR50_I20 <- factor(predClasses_THR70_c3_THR50_I20, levels = c('longSurv', 'shortSurv'))

predClasses_THR70_c3_CTLA4 <- ifelse(Train_prob_THR70_c3_CTLA4 >= thr_THR70_c3_CTLA4$threshold, "longSurv", "shortSurv")
table(predClasses_THR70_c3_CTLA4)
predClasses_THR70_c3_CTLA4 <- factor(predClasses_THR70_c3_CTLA4, levels = c('longSurv', 'shortSurv'))

predClasses_THR70_c3_CTLA4_ICOS_CXCL13 <- ifelse(Train_prob_THR70_c3_CTLA4_ICOS_CXCL13 >= thr_THR70_c3_CTLA4_ICOS_CXCL13$threshold, "longSurv", "shortSurv")
table(predClasses_THR70_c3_CTLA4_ICOS_CXCL13)
predClasses_THR70_c3_CTLA4_ICOS_CXCL13 <- factor(predClasses_THR70_c3_CTLA4_ICOS_CXCL13, levels = c('longSurv', 'shortSurv'))

predClasses_THR70_c3_KLF7 <- ifelse(Train_prob_THR70_c3_KLF7 >= thr_THR70_c3_KLF7$threshold, "longSurv", "shortSurv")
table(predClasses_THR70_c3_KLF7)
predClasses_THR70_c3_KLF7 <- factor(predClasses_THR70_c3_KLF7, levels = c('longSurv', 'shortSurv'))

predClasses_THR70_c3_ZNF627 <- ifelse(Train_prob_THR70_c3_ZNF627 >= thr_THR70_c3_ZNF627$threshold, "longSurv", "shortSurv")
table(predClasses_THR70_c3_ZNF627)
predClasses_THR70_c3_ZNF627 <- factor(predClasses_THR70_c3_ZNF627, levels = c('longSurv', 'shortSurv'))

##########################
## Keep only the relevant information (Metastasis Event and Time)
cluster3_pheno <- cbind(cluster3_pheno[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Relapse.Free.Status", "Relapse.Free.Status..Months.", "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC", "X3.Gene.classifier.subtype", "HER2.Status")], 
                        Train_prob_THR70_c3_THR70_I20, predClasses_THR70_c3_THR70_I20,
                        Train_prob_THR70_c3_THR70_I20_rfe, predClasses_THR70_c3_THR70_I20_rfe,
                        Train_prob_THR70_c3_THR50_I20, predClasses_THR70_c3_THR50_I20,
                        Train_prob_THR70_c3_CTLA4, predClasses_THR70_c3_CTLA4,
                        Train_prob_THR70_c3_CTLA4_ICOS_CXCL13, predClasses_THR70_c3_CTLA4_ICOS_CXCL13,
                        Train_prob_THR70_c3_KLF7, predClasses_THR70_c3_KLF7,
                        Train_prob_THR70_c3_ZNF627, predClasses_THR70_c3_ZNF627
)


CoxData_metabric_c3 <- data.frame(cluster3_pheno)

#########################################
# THR70 I20
#########################################

# OS
Fit_sig_metabric_os_THR70_T1_THR70_I20 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR70_c3_THR70_I20, data = CoxData_metabric_c3)

##########
# RFS
Fit_sig_metabric_RFS_THR70_T1_THR70_I20 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR70_c3_THR70_I20, data = CoxData_metabric_c3)


# RFS COXPH
#CoxData_metabric_c3_forCox <- CoxData_metabric_c3[CoxData_metabric_c3$Train_prob_THR70_c3_THR70_I20]  
Fit_sig_metabric_RFS_THR70_T1_THR70_I20_coxph <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_prob_THR70_c3_THR70_I20, data = CoxData_metabric_c3)
summary(Fit_sig_metabric_RFS_THR70_T1_THR70_I20_coxph)


# plot OS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_os_T1_THR70_I20.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR70_T1_THR70_I20,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR70: THR70 I20 model'
)
dev.off()

#####################
# plot RFS
png("./figures/T1_DE_THR70_RFS/THR70_metabric_RFS_T1_THR70_I20_20yrs.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_sig_metabric_RFS_THR70_T1_THR70_I20,
           risk.table = FALSE,
           pval = FALSE,
           pval.size  = 12,
           xlim = c(0, 240),
           break.x.by = 40,
           legend.labs = c('PQNBC.i-', 'PQNBC.i+'),
           legend.title = c(''),
           ggtheme = theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_line(colour = "grey90"),
                           panel.grid.minor = element_line(colour = "grey90"),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.spacing.x = unit(0.5, "cm"),
                           legend.spacing.y = unit(0.5, "cm"),
                           legend.key.height = unit(1.3, "lines"),
                           axis.title = element_text(size = 14, face = 'bold.italic', color = 'black'),
                           axis.text = element_text(size = 12, face = 'bold.italic', color = 'black'), 
                           legend.text = element_text(size = 16, face = 'bold.italic', color = 'black'),
           ), 
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           #title = 'RFS in METABRIC in T1 class derived from THR70: THR70 I20 model'
)
dev.off()

#####################
# plot RFS by Her2
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_RFS_T1_THR70_I20_20yrs_byHer2.tiff", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_sig_metabric_RFS_THR70_T1_THR70_I20_byHer2,
           risk.table = FALSE,
           pval = TRUE,
           pval.size  = 12,
           xlim = c(0, 240),
           break.x.by = 80,
           legend.labs = c('PQNBC.i+', 'PQNBC.i-'),
           legend.title = c('ER-negative clusters'),
           facet.by = 'HER2.Status', 
           short.panel.labs	= TRUE,
           ggtheme = theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_line(colour = "grey90"),
                           panel.grid.minor = element_line(colour = "grey90"),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.spacing.x = unit(0.5, "cm"),
                           legend.spacing.y = unit(0.5, "cm"),
                           legend.key.height = unit(1.3, "lines"),
                           axis.title = element_text(size = 14, face = 'bold.italic', color = 'black'),
                           axis.text = element_text(size = 12, face = 'bold.italic', color = 'black'), 
                           legend.text = element_text(size = 16, face = 'bold.italic', color = 'black'),
           ), 
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = c('#6057cc', '#8a899a'),
           #title = 'RFS in METABRIC in T1 class derived from THR70: THR70 I20 model'
)
dev.off()


#########################################
# THR70 I20 rfe
#########################################

# OS
Fit_sig_metabric_os_THR70_T1_THR70_I0_rfe <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR70_c3_THR70_I20_rfe, data = CoxData_metabric_c3)

##########
# RFS
Fit_sig_metabric_RFS_THR70_T1_THR70_I20_rfe <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR70_c3_THR70_I20_rfe, data = CoxData_metabric_c3)


# RFS COXPH
#CoxData_metabric_c3_forCox <- CoxData_metabric_c3[CoxData_metabric_c3$Train_prob_THR70_c3_THR70_I20]  
Fit_sig_metabric_RFS_THR70_T1_THR70_I20_rfe_coxph <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_prob_THR70_c3_THR70_I20_rfe, data = CoxData_metabric_c3)
summary(Fit_sig_metabric_RFS_THR70_T1_THR70_I20_rfe_coxph)


# plot OS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_os_T1_THR70_I20_rfe.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR70_T1_THR70_I0_rfe,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR70: THR70 I20 rfe model'
)
dev.off()

#####################
# plot RFS
png("./figures/T1_DE_THR70_RFS/THR70_metabric_RFS_T1_THR70_I20_rfe_20yrs.png", width = 2200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR70_T1_THR70_I20_rfe,
           risk.table = FALSE,
           pval = FALSE,
           pval.size  = 12,
           xlim = c(0, 240),
           break.x.by = 40,
           legend.labs = c('PQNBCi+', 'PQNBCi-'),
           legend.title = c(''),
           ggtheme = theme_survminer(base_size = 25, font.x = c(25, 'bold.italic', 'black'), font.y = c(25, 'bold.italic', 'black'), font.tickslab = c(25, 'plain', 'black'), font.legend = c(25, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           #title = 'RFS in METABRIC in T1 class derived from THR70: THR70 I20 model'
)
dev.off()

#########################################
# THR50 I20
#########################################

# OS
Fit_sig_metabric_os_THR70_T1_THR50_I20 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR70_c3_THR50_I20, data = CoxData_metabric_c3)

# RFS
Fit_sig_metabric_RFS_THR70_T1_THR50_I20 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR70_c3_THR50_I20, data = CoxData_metabric_c3)


# plot OS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_os_T1_THR50_I20.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR70_T1_THR50_I20,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR70: THR50 I20 model'
)
dev.off()

###################
# plot RFS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_RFS_T1_THR50_I20.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR70_T1_THR50_I20,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('PQNBCi+', 'PQNBCi-'),
           legend.title = '',
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR70: THR50 I20 model'
)
dev.off()

#########################################
# CTLA4
#########################################

# OS
Fit_sig_metabric_os_THR70_T1_CTLA4 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR70_c3_CTLA4, data = CoxData_metabric_c3)

# RFS
Fit_sig_metabric_RFS_THR70_T1_CTLA4 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR70_c3_CTLA4, data = CoxData_metabric_c3)


# plot OS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_os_T1_CTLA4.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR70_T1_CTLA4,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR70: CTLA4'
)
dev.off()

################
# plot RFS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_RFS_T1_CTLA4.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR70_T1_CTLA4,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR70: CTLA4'
)
dev.off()


#########################################
# CTLA4_ICOS_CXCL13
#########################################

# OS
Fit_sig_metabric_os_THR70_T1_CTLA4_ICOS_CXCL13 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR70_c3_CTLA4_ICOS_CXCL13, data = CoxData_metabric_c3)

# RFS
Fit_sig_metabric_RFS_THR70_T1_CTLA4_ICOS_CXCL13 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR70_c3_CTLA4_ICOS_CXCL13, data = CoxData_metabric_c3)


# plot OS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_os_T1_CTLA4_ICOS_CXCL13.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR70_T1_CTLA4_ICOS_CXCL13,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR70: CTLA4 + ICOS + CXCL13'
)
dev.off()

######################################
# plot RFS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_RFS_T1_CTLA4_ICOS_CXCL13.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR70_T1_CTLA4_ICOS_CXCL13,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR70: CTLA4 + ICOS + CXCL13'
)
dev.off()

#########################################
# KLF7
#########################################

# OS
Fit_sig_metabric_os_THR70_T1_KLF7 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR70_c3_KLF7, data = CoxData_metabric_c3)

# RFS
Fit_sig_metabric_RFS_THR70_T1_KLF7 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR70_c3_KLF7, data = CoxData_metabric_c3)


# plot OS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_os_T1_KLF7.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR70_T1_KLF7,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR70: KLF7'
)
dev.off()

#####################
# plot RFS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_RFS_T1_KLF7.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR70_T1_KLF7,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR70: KLF7'
)
dev.off()

#########################################
# ZNF627
#########################################

# OS
Fit_sig_metabric_os_THR70_T1_ZNF627 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR70_c3_ZNF627, data = CoxData_metabric_c3)

# RFS
Fit_sig_metabric_RFS_THR70_T1_ZNF627 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR70_c3_ZNF627, data = CoxData_metabric_c3)


# plot OS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_os_T1_ZNF627.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR70_T1_ZNF627,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR70: ZNF627'
)
dev.off()

###################
# plot RFS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_RFS_T1_ZNF627.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR70_T1_ZNF627,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR70: ZNF627'
)
dev.off()

##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################
## recombine c3 with the rest 
cluster3_pheno$THR_clusters_THR70_I20 <- as.factor(cluster3_pheno$predClasses_THR70_c3_THR70_I20)
levels(cluster3_pheno$THR_clusters_THR70_I20) <- c('PNBC_A', 'PNBC_B')

cluster3_pheno$THR_clusters_THR50_I20 <- as.factor(cluster3_pheno$predClasses_THR70_c3_THR50_I20)
levels(cluster3_pheno$THR_clusters_THR50_I20) <-  c('PNBC_A', 'PNBC_B')

cluster3_pheno$THR_clusters_CTLA4 <- as.factor(cluster3_pheno$predClasses_THR70_c3_CTLA4)
levels(cluster3_pheno$THR_clusters_CTLA4) <-  c('PNBC_A', 'PNBC_B')

cluster3_pheno$THR_clusters_CTLA4_ICOS_CXCL13 <- as.factor(cluster3_pheno$predClasses_THR70_c3_CTLA4_ICOS_CXCL13)
levels(cluster3_pheno$THR_clusters_CTLA4_ICOS_CXCL13) <-  c('PNBC_A', 'PNBC_B')

cluster3_pheno$THR_clusters_KLF7 <- as.factor(cluster3_pheno$predClasses_THR70_c3_KLF7)
levels(cluster3_pheno$THR_clusters_KLF7) <-  c('PNBC_A', 'PNBC_B')

cluster3_pheno$THR_clusters_ZNF627 <- as.factor(cluster3_pheno$predClasses_THR70_c3_ZNF627)
levels(cluster3_pheno$THR_clusters_ZNF627) <-  c('PNBC_A', 'PNBC_B')


T1 <- data.frame(THR_clusters_THR70_I20 = cluster3_pheno$THR_clusters_THR70_I20,
                 THR_clusters_THR50_I20 = cluster3_pheno$THR_clusters_THR50_I20,
                 THR_clusters_CTLA4 = cluster3_pheno$THR_clusters_CTLA4,
                 THR_clusters_CTLA4_ICOS_CXCL13 = cluster3_pheno$THR_clusters_CTLA4_ICOS_CXCL13,
                 THR_clusters_KLF7 = cluster3_pheno$THR_clusters_KLF7,
                 THR_clusters_ZNF627 = cluster3_pheno$THR_clusters_ZNF627,
                 `Sample.ID` = rownames(cluster3_pheno))

rownames(T1) <- rownames(cluster3_pheno)


# merge
Pheno_metabric$`THR clusters`[Pheno_metabric$`THR clusters` == 'PNBC'] <- NA


Pheno_metabric2 <- merge(x = T1, y = Pheno_metabric, by="Sample.ID", all.y = TRUE)

Pheno_metabric2 <- Pheno_metabric2 %>% 
  mutate(`THR clusters` = as.factor(`THR clusters`), THR_clusters_THR70_I20 = as.factor(THR_clusters_THR70_I20), THR_clusters_THR50_I20 = as.factor(THR_clusters_THR50_I20)) %>%
  mutate(THR.clusters_THR70_I20_Merged = coalesce(THR_clusters_THR70_I20,`THR clusters`)) %>%
  mutate(THR.clusters_THR50_I20_Merged = coalesce(THR_clusters_THR50_I20,`THR clusters`)) %>%
  mutate(THR.clusters_CTLA4Merged = coalesce(THR_clusters_CTLA4,`THR clusters`)) %>%
  mutate(THR.clusters_CTLA4_ICOS_CXCL13Merged = coalesce(THR_clusters_CTLA4_ICOS_CXCL13,`THR clusters`)) %>%
  mutate(THR.clusters_KLF7Merged = coalesce(THR_clusters_KLF7,`THR clusters`)) %>%
  mutate(THR.clusters_ZNF627Merged = coalesce(THR_clusters_ZNF627,`THR clusters`)) 

###########################################################################################
# a table/venn diagram of clinical covariates in each THR70 cluster
Pheno_metabric3 <- Pheno_metabric2
levels(Pheno_metabric3$THR.clusters_THR70_I20_Merged) <- c("PNBC-A", "PNBC-B", "E2", "E2", "E1", "E3", "PNBC")

Pheno_metabric3$THR.clusters_THR70_I20_Merged <- factor(Pheno_metabric3$THR.clusters_THR70_I20_Merged, levels = c('E1', 'E2', 'E3', 'PNBC-A', 'PNBC-B'))
# sumtable(Pheno_metabric3,
#          vars = c('ER.status.measured.by.IHC', 'PR.Status', 'HER2.Status', 'Pam50...Claudin.low.subtype', 'X3.Gene.classifier.subtype'),
#          group = 'THR.clusters_THR70_I20_Merged',
#          file='metabric_clusters_summary',
#          out = 'csv',
#          title='METABRIC clusters Summary Statistics',
#          simple.kable=FALSE,
#          opts=list())

###########################################################################################
##########################################################################################
## survival analysis

## Keep only the relevant information (Metastasis Event and Time)
survival_metabric <- Pheno_metabric2[, c("Overall.Survival.Status", "Overall.Survival..Months.", 
                                         "Relapse.Free.Status", "Relapse.Free.Status..Months.", 
                                         "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC", 'HER2.Status',
                                         "X3.Gene.classifier.subtype", 
                                         "THR.clusters_THR70_I20_Merged",
                                         "THR.clusters_THR50_I20_Merged",
                                         'THR.clusters_CTLA4Merged',
                                         'THR.clusters_CTLA4_ICOS_CXCL13Merged',
                                         'THR.clusters_KLF7Merged',
                                         'THR.clusters_ZNF627Merged'
)] 

survival_metabric$THR.clusters_THR70_I20_Merged <- as.factor(survival_metabric$THR.clusters_THR70_I20_Merged)
survival_metabric$THR.clusters_THR50_I20_Merged <- as.factor(survival_metabric$THR.clusters_THR50_I20_Merged)
survival_metabric$THR.clusters_CTLA4Merged <- as.factor(survival_metabric$THR.clusters_CTLA4Merged)
survival_metabric$THR.clusters_CTLA4_ICOS_CXCL13Merged <- as.factor(survival_metabric$THR.clusters_CTLA4_ICOS_CXCL13Merged)
survival_metabric$THR.clusters_KLF7Merged <- as.factor(survival_metabric$THR.clusters_KLF7Merged)
survival_metabric$THR.clusters_ZNF627Merged <- as.factor(survival_metabric$THR.clusters_ZNF627Merged)

survival_metabric$THR.clusters_THR70_I20_Merged <- droplevels(survival_metabric$THR.clusters_THR70_I20_Merged)
survival_metabric$THR.clusters_THR50_I20_Merged <- droplevels(survival_metabric$THR.clusters_THR50_I20_Merged)
survival_metabric$THR.clusters_CTLA4Merged <- droplevels(survival_metabric$THR.clusters_CTLA4Merged)
survival_metabric$THR.clusters_CTLA4_ICOS_CXCL13Merged <- droplevels(survival_metabric$THR.clusters_CTLA4_ICOS_CXCL13Merged)
survival_metabric$THR.clusters_KLF7Merged <- droplevels(survival_metabric$THR.clusters_KLF7Merged)
survival_metabric$THR.clusters_ZNF627Merged <- droplevels(survival_metabric$THR.clusters_ZNF627Merged)

cluster_colors <- c('#6057cc', '#8a899a', '#66a61e', "#E7298A", "#1B9E77", "#D95F02")


# fix the levels
#levels(survival_metabric$THR.clusters_THR70_I20_Merged) <- c("PNBC_B", "PNBC_A", 'E2', 'E2', 'E1', 'E3') 

################################################################
## Survival curves: THR70 I20
################################################################
survival_metabric2 <- survival_metabric
levels(survival_metabric2$THR.clusters_THR70_I20_Merged) <- c("PNBC_A", "PNBC_B", 'E2', 'E2', 'E1', 'E3') 

# re-order the levels
survival_metabric$THR.clusters_THR70_I20_Merged <- factor(survival_metabric$THR.clusters_THR70_I20_Merged, levels = c('E1', 'E2a', 'E2b', 'E3', 'PNBC_A', 'PNBC_B'))
survival_metabric2$THR.clusters_THR70_I20_Merged <- factor(survival_metabric2$THR.clusters_THR70_I20_Merged, levels = c('E1', 'E2', 'E3', 'PNBC_A', 'PNBC_B'))

# colors
cluster_colors <- c('#1B9E77', "#E7298A", "#66A61E", "#D95F02", '#6057cc', '#8a899a')
cluster_colors2 <- c('#1B9E77', "#c1b026", "#D95F02", '#6057cc', '#8a899a')

# OS
Fit_metabric_os_THR70_I20 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_THR70_I20_Merged, data = survival_metabric)
Fit_metabric_os_THR70_I20_2 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_THR70_I20_Merged, data = survival_metabric2)

# RFS
Fit_metabric_RFS_THR70_I20 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_THR70_I20_Merged, data = survival_metabric)
Fit_metabric_RFS_THR70_I20_2 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_THR70_I20_Merged, data = survival_metabric2)

png("./figures/T1_DE_THR70_RFS/metabric_os_5clusters_THR70_I20_merged_20yrs.png", width = 2500, height = 2500, res = 300)
ggsurvplot(Fit_metabric_os_THR70_I20,
           risk.table = FALSE,
           pval = FALSE,
           palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = levels(survival_metabric$THR.clusters_THR70_I20_Merged),
           legend.title	= 'THR-70 clusters',
           pval.size = 12,
           break.x.by = 40,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE 
           #title = 'THR70 clusters and OS: THR70 + I20'
           )
dev.off()

## RFS: 
png("./figures/T1_DE_THR70_RFS/metabric_RFS_5clusters_THR70_I20_merged_20yrs.png", width = 2500, height = 2500, res = 300)
ggsurvplot(Fit_metabric_RFS_THR70_I20,
           risk.table = FALSE,
           pval = FALSE,
           palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = c("E1", "E2a", "E2b","E3", "PQNBCi+", "PQNBCi-"),
           legend.title	= '',
           pval.size = 12,
           break.x.by = 40,
           ggtheme = theme_survminer(base_size = 20, font.x = c(20, 'bold.italic', 'black'), font.y = c(20, 'bold.italic', 'black'), font.tickslab = c(20, 'plain', 'black'), font.legend = c(20, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR70 clusters and RFS: THR70 + I20'
           )
dev.off()

###############
# same with merged E2a and E2b
pdf("./figures/T1_DE_THR70_RFS/metabric_os_5clusters_THR70_I20_merged_20yrs_E2merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_THR70_I20_2,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors2,
           xlim = c(0,240),
           legend.labs = c("E1", "E2","E3", "Pi+", "Pi-"),
           legend.title	= '',
           pval.size = 12,
           break.x.by = 40,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR70 clusters and OS: THR70 + I20'
           )
dev.off()

## RFS: 

png("./figures/T1_DE_THR70_RFS/metabric_RFS_5clusters_THR70_I20_merged_20yrs_E2merged.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_metabric_RFS_THR70_I20_2,
           risk.table = FALSE,
           pval = FALSE,
           palette =  cluster_colors2,
           xlim = c(0,240),
           legend.labs = c("E1", "E2","E3", "Pi-", "Pi+"),
           legend.title	= '',
           pval.size = 11,
           break.x.by = 40,
           ggtheme = theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_line(colour = "grey90"),
                           panel.grid.minor = element_line(colour = "grey90"),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.spacing.x = unit(0.5, "cm"),
                           legend.spacing.y = unit(0.5, "cm"),
                           legend.key.height = unit(1.3, "lines"),
                           axis.title = element_text(size = 14, face = 'bold.italic', color = 'black'),
                           axis.text = element_text(size = 12, face = 'bold.italic', color = 'black'), 
                           legend.text = element_text(size = 16, face = 'bold.italic', color = 'black'),
           ), 
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR70 clusters and RFS: THR70 + I20'
           ) + guides(
             colour = guide_legend(ncol = 2))
dev.off()



##########################################################################################
# RFS HER2 by THR clusters
##########################################################################################

survival_metabric$HER2.Status <- as.factor(survival_metabric$HER2.Status)
table(survival_metabric$HER2.Status)
survival_metabric_Her2Pos <- survival_metabric[survival_metabric$HER2.Status == 'Positive', ]

#survival_metabric_Her2Pos <- survival_metabric[survival_metabric$X3.Gene.classifier.subtype == 'HER2+', ]
#survival_metabric_Her2Pos <- survival_metabric_Her2Pos[!is.na(survival_metabric_Her2Pos$X3.Gene.classifier.subtype), ]

# change the levels of THR clusters to group all the E together
levels(survival_metabric_Her2Pos$THR.clusters_THR70_I20_Merged) <- c('E1, E2, E3', 'E1, E2, E3', 'E1, E2, E3', 'E1, E2, E3', 'PNBC', 'PNBC')

# fit survival curves
Fit_sig_metabric_RFS_HER2pos <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_THR70_I20_Merged, data = survival_metabric_Her2Pos)

# plot
png("./figures/T1_DE_THR70_RFS/HER2_by_THR70clusters.png", width = 2200, height = 2000, res = 250)
ggsurvplot(Fit_sig_metabric_RFS_HER2pos,
           risk.table = FALSE,
           pval = TRUE,
           #palette =  cluster_colors,
           xlim = c(0,240),
           legend.labs = gsub('_', '-', levels(survival_metabric_Her2Pos$THR.clusters_THR70_I20_Merged)),
           legend.title	= 'THR-70 clusters',
           pval.size = 11,
           break.x.by = 40,
           ggtheme = theme_survminer(base_size = 18, font.x = c(20, 'bold.italic', 'black'), font.y = c(20, 'bold.italic', 'black'), font.tickslab = c(20, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR70 clusters and RFS: THR70 + I20'
)
dev.off()

##########################################################################################
# RFS PNBC HER2+ versus PNBC HER2- 
##########################################################################################

survival_metabric$THR.clusters_THR70_I20_Merged <- as.factor(survival_metabric$THR.clusters_THR70_I20_Merged)
table(survival_metabric$THR.clusters_THR70_I20_Merged)
survival_metabric_PNBC <- survival_metabric[survival_metabric$THR.clusters_THR70_I20_Merged %in% c('PNBC_B', 'PNBC_A'), ]

#survival_metabric_Her2Pos <- survival_metabric[survival_metabric$X3.Gene.classifier.subtype == 'HER2+', ]
#survival_metabric_Her2Pos <- survival_metabric_Her2Pos[!is.na(survival_metabric_Her2Pos$X3.Gene.classifier.subtype), ]

# change the levels of THR clusters to group all the E together
survival_metabric_PNBC$THR.clusters_THR70_I20_Merged <- droplevels(survival_metabric_PNBC$THR.clusters_THR70_I20_Merged)
levels(survival_metabric_PNBC$THR.clusters_THR70_I20_Merged) <- c('PNBC', 'PNBC')

# change the levels of HER2 status
survival_metabric_PNBC$HER2.Status <- factor(survival_metabric_PNBC$HER2.Status, levels = c('Positive', 'Negative'))

# fit survival curves
Fit_sig_metabric_RFS_PNBC_byHER2 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ HER2.Status, data = survival_metabric_PNBC)

# RFS COXPH

survival_metabric_PNBC_2 <- survival_metabric_PNBC
survival_metabric_PNBC_2$HER2.Status <- factor(survival_metabric_PNBC_2$HER2.Status, levels = c('Negative', 'Positive'))
table(survival_metabric_PNBC_2$HER2.Status)
Fit_sig_metabric_RFS_PNBC_byHER2_coxph <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ HER2.Status, data = survival_metabric_PNBC_2)
summary(Fit_sig_metabric_RFS_PNBC_byHER2_coxph)

# plot
png("./figures/T1_DE_THR70_RFS/THR70_PNBC_byHER2.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_sig_metabric_RFS_PNBC_byHER2,
           risk.table = FALSE,
           pval = FALSE,
           palette =  'jco',
           xlim = c(0,240),
           legend.labs = c('HER2+', 'HER2-'),
           legend.title	= '',
           pval.size = 11,
           break.x.by = 40,
           ggtheme = theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_line(colour = "grey90"),
                           panel.grid.minor = element_line(colour = "grey90"),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.spacing.x = unit(0.5, "cm"),
                           legend.spacing.y = unit(0.5, "cm"),
                           legend.key.height = unit(1.3, "lines"),
                           axis.title = element_text(size = 14, face = 'bold.italic', color = 'black'),
                           axis.text = element_text(size = 12, face = 'bold.italic', color = 'black'), 
                           legend.text = element_text(size = 16, face = 'bold.italic', color = 'black'),
           ), 
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR70 clusters and RFS: THR70 + I20'
)
dev.off()


##########################################################################################
# RFS E1 HER2+ versus E1 HER2- 
##########################################################################################

survival_metabric$THR.clusters_THR70_I20_Merged <- as.factor(survival_metabric$THR.clusters_THR70_I20_Merged)
table(survival_metabric$THR.clusters_THR70_I20_Merged)
survival_metabric_E1 <- survival_metabric[survival_metabric$THR.clusters_THR70_I20_Merged == 'E1', ]

#survival_metabric_Her2Pos <- survival_metabric[survival_metabric$X3.Gene.classifier.subtype == 'HER2+', ]
#survival_metabric_Her2Pos <- survival_metabric_Her2Pos[!is.na(survival_metabric_Her2Pos$X3.Gene.classifier.subtype), ]

# change the levels of THR clusters to group all the E together
survival_metabric_E1$THR.clusters_THR70_I20_Merged <- droplevels(survival_metabric_E1$THR.clusters_THR70_I20_Merged)

# change the levels of HER2 status
survival_metabric_E1$HER2.Status <- factor(survival_metabric_E1$HER2.Status, levels = c('Positive', 'Negative'))

# fit survival curves
Fit_sig_metabric_RFS_E1_byHER2 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ HER2.Status, data = survival_metabric_E1)

# RFS COXPH

survival_metabric_E1_coxph <- survival_metabric_E1
survival_metabric_E1_coxph$HER2.Status <- factor(survival_metabric_E1_coxph$HER2.Status, levels = c('Negative', 'Positive'))
table(survival_metabric_E1_coxph$HER2.Status)
Fit_sig_metabric_RFS_E1_byHER2_coxph <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ HER2.Status, data = survival_metabric_E1_coxph)
summary(Fit_sig_metabric_RFS_E1_byHER2_coxph)

# plot
png("./figures/T1_DE_THR70_RFS/THR70_E1_byHER2.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_sig_metabric_RFS_E1_byHER2,
           risk.table = FALSE,
           pval = TRUE,
           palette =  'jco',
           xlim = c(0,240),
           legend.labs = c('HER2+', 'HER2-'),
           legend.title	= '',
           pval.size = 11,
           break.x.by = 40,
           ggtheme = theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_line(colour = "grey90"),
                           panel.grid.minor = element_line(colour = "grey90"),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.spacing.x = unit(0.5, "cm"),
                           legend.spacing.y = unit(0.5, "cm"),
                           legend.key.height = unit(1.3, "lines"),
                           axis.title = element_text(size = 14, face = 'bold.italic', color = 'black'),
                           axis.text = element_text(size = 12, face = 'bold.italic', color = 'black'), 
                           legend.text = element_text(size = 16, face = 'bold.italic', color = 'black'),
           ), 
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR70 clusters and RFS: THR70 + I20'
)
dev.off()

##########################################################################################
# RFS E2 HER2+ versus E2 HER2- 
##########################################################################################

survival_metabric$THR.clusters_THR70_I20_Merged <- as.factor(survival_metabric$THR.clusters_THR70_I20_Merged)
table(survival_metabric$THR.clusters_THR70_I20_Merged)
survival_metabric_E2 <- survival_metabric[survival_metabric$THR.clusters_THR70_I20_Merged %in% c('E2a', 'E2b'), ]

#survival_metabric_Her2Pos <- survival_metabric[survival_metabric$X3.Gene.classifier.subtype == 'HER2+', ]
#survival_metabric_Her2Pos <- survival_metabric_Her2Pos[!is.na(survival_metabric_Her2Pos$X3.Gene.classifier.subtype), ]

# change the levels of THR clusters to group all the E together
survival_metabric_E2$THR.clusters_THR70_I20_Merged <- droplevels(survival_metabric_E2$THR.clusters_THR70_I20_Merged)

# change the levels of HER2 status
survival_metabric_E2$HER2.Status <- factor(survival_metabric_E2$HER2.Status, levels = c('Positive', 'Negative'))

# fit survival curves
Fit_sig_metabric_RFS_E2_byHER2 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ HER2.Status, data = survival_metabric_E2)

# RFS COXPH

survival_metabric_E2_coxph <- survival_metabric_E2
survival_metabric_E2_coxph$HER2.Status <- factor(survival_metabric_E2_coxph$HER2.Status, levels = c('Negative', 'Positive'))
table(survival_metabric_E2_coxph$HER2.Status)
Fit_sig_metabric_RFS_E2_byHER2_coxph <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ HER2.Status, data = survival_metabric_E2_coxph)
summary(Fit_sig_metabric_RFS_E2_byHER2_coxph)

# plot
png("./figures/T1_DE_THR70_RFS/THR70_E2_byHER2.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_sig_metabric_RFS_E2_byHER2,
           risk.table = FALSE,
           pval = TRUE,
           palette =  'jco',
           xlim = c(0,240),
           legend.labs = c('HER2+', 'HER2-'),
           legend.title	= '',
           pval.size = 11,
           break.x.by = 40,
           ggtheme = theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_line(colour = "grey90"),
                           panel.grid.minor = element_line(colour = "grey90"),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.spacing.x = unit(0.5, "cm"),
                           legend.spacing.y = unit(0.5, "cm"),
                           legend.key.height = unit(1.3, "lines"),
                           axis.title = element_text(size = 14, face = 'bold.italic', color = 'black'),
                           axis.text = element_text(size = 12, face = 'bold.italic', color = 'black'), 
                           legend.text = element_text(size = 16, face = 'bold.italic', color = 'black'),
           ), 
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR70 clusters and RFS: THR70 + I20'
)
dev.off()

##########################################################################################
# RFS E3 HER2+ versus E3 HER2- 
##########################################################################################

survival_metabric$THR.clusters_THR70_I20_Merged <- as.factor(survival_metabric$THR.clusters_THR70_I20_Merged)
table(survival_metabric$THR.clusters_THR70_I20_Merged)
survival_metabric_E3 <- survival_metabric[survival_metabric$THR.clusters_THR70_I20_Merged == 'E3', ]

#survival_metabric_Her2Pos <- survival_metabric[survival_metabric$X3.Gene.classifier.subtype == 'HER2+', ]
#survival_metabric_Her2Pos <- survival_metabric_Her2Pos[!is.na(survival_metabric_Her2Pos$X3.Gene.classifier.subtype), ]

# change the levels of THR clusters to group all the E together
survival_metabric_E3$THR.clusters_THR70_I20_Merged <- droplevels(survival_metabric_E3$THR.clusters_THR70_I20_Merged)

# change the levels of HER2 status
survival_metabric_E3$HER2.Status <- factor(survival_metabric_E3$HER2.Status, levels = c('Positive', 'Negative'))

# fit survival curves
Fit_sig_metabric_RFS_E3_byHER2 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ HER2.Status, data = survival_metabric_E3)

# RFS COXPH

survival_metabric_E3_coxph <- survival_metabric_E3
survival_metabric_E3_coxph$HER2.Status <- factor(survival_metabric_E3_coxph$HER2.Status, levels = c('Negative', 'Positive'))
table(survival_metabric_E3_coxph$HER2.Status)
Fit_sig_metabric_RFS_E3_byHER2_coxph <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ HER2.Status, data = survival_metabric_E3_coxph)
summary(Fit_sig_metabric_RFS_E3_byHER2_coxph)

# plot
png("./figures/T1_DE_THR70_RFS/THR70_E3_byHER2.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_sig_metabric_RFS_E3_byHER2,
           risk.table = FALSE,
           pval = TRUE,
           palette =  'jco',
           xlim = c(0,240),
           legend.labs = c('HER2+', 'HER2-'),
           legend.title	= '',
           pval.size = 11,
           break.x.by = 40,
           ggtheme = theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_line(colour = "grey90"),
                           panel.grid.minor = element_line(colour = "grey90"),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.spacing.x = unit(0.5, "cm"),
                           legend.spacing.y = unit(0.5, "cm"),
                           legend.key.height = unit(1.3, "lines"),
                           axis.title = element_text(size = 14, face = 'bold.italic', color = 'black'),
                           axis.text = element_text(size = 12, face = 'bold.italic', color = 'black'), 
                           legend.text = element_text(size = 16, face = 'bold.italic', color = 'black'),
           ), 
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR70 clusters and RFS: THR70 + I20'
)
dev.off()

##########################################################################################
# RFS all E HER2+ versus all E HER2- 
##########################################################################################

survival_metabric$THR.clusters_THR70_I20_Merged <- as.factor(survival_metabric$THR.clusters_THR70_I20_Merged)
table(survival_metabric$THR.clusters_THR70_I20_Merged)
survival_metabric_allE <- survival_metabric[survival_metabric$THR.clusters_THR70_I20_Merged %in% c('E1', 'E2a', 'E2b', 'E3'), ]

#survival_metabric_Her2Pos <- survival_metabric[survival_metabric$X3.Gene.classifier.subtype == 'HER2+', ]
#survival_metabric_Her2Pos <- survival_metabric_Her2Pos[!is.na(survival_metabric_Her2Pos$X3.Gene.classifier.subtype), ]

# change the levels of THR clusters to group all the E together
survival_metabric_allE$THR.clusters_THR70_I20_Merged <- droplevels(survival_metabric_allE$THR.clusters_THR70_I20_Merged)

# change the levels of HER2 status
survival_metabric_allE$HER2.Status <- factor(survival_metabric_allE$HER2.Status, levels = c('Positive', 'Negative'))

# fit survival curves
Fit_sig_metabric_RFS_allE_byHER2 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ HER2.Status, data = survival_metabric_allE)

# RFS COXPH

survival_metabric_allE_2 <- survival_metabric_allE
survival_metabric_allE_2$HER2.Status <- factor(survival_metabric_allE_2$HER2.Status, levels = c('Negative', 'Positive'))
table(survival_metabric_allE_2$HER2.Status)
Fit_sig_metabric_RFS_allE_byHER2_coxph <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ HER2.Status, data = survival_metabric_allE_2)
summary(Fit_sig_metabric_RFS_allE_byHER2_coxph)


# plot
png("./figures/T1_DE_THR70_RFS/THR70_allE_byHER2.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_sig_metabric_RFS_allE_byHER2,
           risk.table = FALSE,
           pval = FALSE,
           palette =  'jco',
           xlim = c(0,240),
           legend.labs = c('HER2+', 'HER2-'),
           legend.title	= '',
           pval.size = 11,
           break.x.by = 40,
           ggtheme = theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_line(colour = "grey90"),
                           panel.grid.minor = element_line(colour = "grey90"),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.spacing.x = unit(0.5, "cm"),
                           legend.spacing.y = unit(0.5, "cm"),
                           legend.key.height = unit(1.3, "lines"),
                           axis.title = element_text(size = 14, face = 'bold.italic', color = 'black'),
                           axis.text = element_text(size = 12, face = 'bold.italic', color = 'black'), 
                           legend.text = element_text(size = 16, face = 'bold.italic', color = 'black'),
           ), 
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR70 clusters and RFS: THR70 + I20'
)
dev.off()

##########################################################################################
# RFS E1 vs PNBC in HER2+ groups 
##########################################################################################

survival_metabric$THR.clusters_THR70_I20_Merged <- as.factor(survival_metabric$THR.clusters_THR70_I20_Merged)
table(survival_metabric$THR.clusters_THR70_I20_Merged)
survival_metabric_E1_PNBC_HER2pos <- survival_metabric[survival_metabric$THR.clusters_THR70_I20_Merged %in% c('E1', 'PNBC_A', 'PNBC_B') & survival_metabric$HER2.Status == 'Positive', ]

#survival_metabric_Her2Pos <- survival_metabric[survival_metabric$X3.Gene.classifier.subtype == 'HER2+', ]
#survival_metabric_Her2Pos <- survival_metabric_Her2Pos[!is.na(survival_metabric_Her2Pos$X3.Gene.classifier.subtype), ]

# fix the levels of THR clusters
survival_metabric_E1_PNBC_HER2pos$THR.clusters_THR70_I20_Merged <- droplevels(survival_metabric_E1_PNBC_HER2pos$THR.clusters_THR70_I20_Merged)
levels(survival_metabric_E1_PNBC_HER2pos$THR.clusters_THR70_I20_Merged) <- c('E1', 'PNBC', 'PNBC')

# fit survival curves
Fit_sig_metabric_RFS_HER2pos_E1versusPNBC <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_THR70_I20_Merged, data = survival_metabric_E1_PNBC_HER2pos)

# plot
png("./figures/T1_DE_THR70_RFS/THR70_HER2pos_E1versusPNBC.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_sig_metabric_RFS_HER2pos_E1versusPNBC,
           risk.table = FALSE,
           pval = TRUE,
           palette =  'jco',
           xlim = c(0,240),
           legend.labs = c('E1 HER2+', 'PQNBC HER2+'),
           #legend.title	= 'THR-70 clusters',
           pval.size = 11,
           break.x.by = 40,
           ggtheme = theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_line(colour = "grey90"),
                           panel.grid.minor = element_line(colour = "grey90"),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.spacing.x = unit(0.5, "cm"),
                           legend.spacing.y = unit(0.5, "cm"),
                           legend.key.height = unit(1.3, "lines"),
                           axis.title = element_text(size = 14, face = 'bold.italic', color = 'black'),
                           axis.text = element_text(size = 12, face = 'bold.italic', color = 'black'), 
                           legend.text = element_text(size = 16, face = 'bold.italic', color = 'black'),
           ), 
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR70 clusters and RFS: THR70 + I20'
)
dev.off()


##########################################################################################
## survival analysis for just ER positive clusters
##########################################################################################

# keep only ER positive clusters
survival_metabric_ERpos <- survival_metabric2[survival_metabric2$THR.clusters_THR70_I20_Merged %in% c("E2", "E1", "E3"), ] 
survival_metabric_ERpos$THR.clusters_THR70_I20_Merged <- droplevels(survival_metabric_ERpos$THR.clusters_THR70_I20_Merged)


# OS
Fit_sig_metabric_os_THR70_ERpos_THR70_I20 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_THR70_I20_Merged, data = survival_metabric_ERpos)

##########
# RFS
Fit_sig_metabric_RFS_THR70_ERpos_THR70_I20 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_THR70_I20_Merged, data = survival_metabric_ERpos)

# RFS COXPH

survival_metabric_ERpos2 <- survival_metabric_ERpos
survival_metabric_ERpos2$THR.clusters_THR70_I20_Merged <- factor(survival_metabric_ERpos2$THR.clusters_THR70_I20_Merged, levels = c('E3', 'E2', 'E1'))
Fit_sig_metabric_RFS_THR70_ERpos_THR70_I20_coxph <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_THR70_I20_Merged, data = survival_metabric_ERpos2)
summary(Fit_sig_metabric_RFS_THR70_ERpos_THR70_I20_coxph)


# RFS by Her2
#Fit_sig_metabric_RFS_THR70_ERpos_THR70_I20_byHer2 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_THR70_I20_Merged + HER2.Status, data = survival_metabric_ERpos)


# plot OS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_os_T1_THR70_I20_ERposClusters.tiff", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_sig_metabric_os_THR70_ERpos_THR70_I20,
           risk.table = FALSE,
           pval = FALSE,
           pval.size = 12,
           xlim = c(0,240),
           break.x.by = 40,
           legend.labs = levels(survival_metabric_ERpos$THR.clusters_THR70_I20_Merged),
           ggtheme = theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_line(colour = "grey90"),
                           panel.grid.minor = element_line(colour = "grey90"),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.spacing.x = unit(0.5, "cm"),
                           legend.spacing.y = unit(0.5, "cm"),
                           legend.key.height = unit(1.3, "lines"),
                           axis.title = element_text(size = 14, face = 'bold.italic', color = 'black'),
                           axis.text = element_text(size = 12, face = 'bold.italic', color = 'black'), 
                           legend.text = element_text(size = 16, face = 'bold.italic', color = 'black'),
           ), 
           risk.table.y.text.col = FALSE,
           palette = c('#1B9E77', "#c1b026", "#D95F02"),
           risk.table.y.text = FALSE, 
           #title = 'OS in METABRIC in T1 class derived from THR70: THR70 I20 model'
)
dev.off()

#####################
# plot RFS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_RFS_T1_THR70_I20_20yrs_ERposClusters.tiff", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_sig_metabric_RFS_THR70_ERpos_THR70_I20,
           risk.table = FALSE,
           pval = FALSE,
           pval.size  = 12,
           xlim = c(0, 240),
           break.x.by = 40,
           legend.labs = levels(survival_metabric_ERpos$THR.clusters_THR70_I20_Merged),
           legend.title = c(''),
           ggtheme = theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_line(colour = "grey90"),
                           panel.grid.minor = element_line(colour = "grey90"),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.spacing.x = unit(0.5, "cm"),
                           legend.spacing.y = unit(0.5, "cm"),
                           legend.key.height = unit(1.3, "lines"),
                           axis.title = element_text(size = 14, face = 'bold.italic', color = 'black'),
                           axis.text = element_text(size = 12, face = 'bold.italic', color = 'black'), 
                           legend.text = element_text(size = 16, face = 'bold.italic', color = 'black'),
           ), 
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = c('#1B9E77', "#c1b026", "#D95F02"),
           #title = 'RFS in METABRIC in T1 class derived from THR70: THR70 I20 model'
)
dev.off()



##########################################################################################
## survival analysis X3 TNBC vs HER2+
##########################################################################################

# keep only X3 TNBC and HER2+
table(survival_metabric2$X3.Gene.classifier.subtype)
survival_metabric_X3_TNBC_HER2 <- survival_metabric2[survival_metabric2$X3.Gene.classifier.subtype %in% c('ER-/HER2-', 'HER2+'), ] 
survival_metabric_X3_TNBC_HER2$X3.Gene.classifier.subtype <- factor(survival_metabric_X3_TNBC_HER2$X3.Gene.classifier.subtype, levels = c('HER2+', 'ER-/HER2-'))
table(survival_metabric_X3_TNBC_HER2$X3.Gene.classifier.subtype)

# OS
Fit_sig_metabric_os_THR70_X3tnbcHer2_THR70_I20 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ X3.Gene.classifier.subtype, data = survival_metabric_X3_TNBC_HER2)

##########
# RFS
Fit_sig_metabric_RFS_THR70_X3tnbcHer2_THR70_I20 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ X3.Gene.classifier.subtype, data = survival_metabric_X3_TNBC_HER2)

# RFS COXPH

survival_metabric_X3_TNBC_HER2_2 <- survival_metabric_X3_TNBC_HER2
survival_metabric_X3_TNBC_HER2_2$X3.Gene.classifier.subtype <- factor(survival_metabric_X3_TNBC_HER2_2$X3.Gene.classifier.subtype, levels = c('ER-/HER2-', 'HER2+'))
table(survival_metabric_X3_TNBC_HER2_2$X3.Gene.classifier.subtype)
Fit_sig_metabric_RFS_THR70_X3tnbcHer2_THR70_I20_coxph <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ X3.Gene.classifier.subtype, data = survival_metabric_X3_TNBC_HER2_2)
summary(Fit_sig_metabric_RFS_THR70_X3tnbcHer2_THR70_I20_coxph)


# RFS by Her2
#Fit_sig_metabric_RFS_THR70_ERpos_THR70_I20_byHer2 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_THR70_I20_Merged + HER2.Status, data = survival_metabric_ERpos)


# plot OS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_os_T1_THR70_I20_X3tnbcHER2.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR70_X3tnbcHer2_THR70_I20,
           risk.table = FALSE,
           pval = FALSE,
           pval.size = 12,
           xlim = c(0,240),
           break.x.by = 40,
           legend.labs = c('HER2+', 'TNBC'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           #palette = c('#1B9E77', "#c1b026", "#D95F02"),
           risk.table.y.text = FALSE, 
           #title = 'OS in METABRIC in T1 class derived from THR70: THR70 I20 model'
)
dev.off()

#####################
# plot RFS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_rfs_T1_THR70_I20_X3tnbcHER2.tiff", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_sig_metabric_RFS_THR70_X3tnbcHer2_THR70_I20,
           risk.table = FALSE,
           pval = FALSE,
           pval.size  = 12,
           xlim = c(0, 240),
           break.x.by = 40,
           legend.labs = c('HER2+', 'TNBC'),
           legend.title = c(''),
           ggtheme = theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_line(colour = "grey90"),
                           panel.grid.minor = element_line(colour = "grey90"),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.spacing.x = unit(0.5, "cm"),
                           legend.spacing.y = unit(0.5, "cm"),
                           legend.key.height = unit(1.3, "lines"),
                           axis.title = element_text(size = 14, face = 'bold.italic', color = 'black'),
                           axis.text = element_text(size = 12, face = 'bold.italic', color = 'black'), 
                           legend.text = element_text(size = 16, face = 'bold.italic', color = 'black'),
           ), 
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #palette = c('#1B9E77', "#c1b026", "#D95F02"),
           #title = 'RFS in METABRIC in T1 class derived from THR70: THR70 I20 model'
)
dev.off()

##########################################################################################
## survival analysis X3 ER+ high vs low prolif
##########################################################################################

# keep only ER+ high and low prolif
table(survival_metabric2$X3.Gene.classifier.subtype)
survival_metabric_X3_ERprolif <- survival_metabric2[survival_metabric2$X3.Gene.classifier.subtype %in% c('ER+/HER2- High Prolif', 'ER+/HER2- Low Prolif'), ] 
survival_metabric_X3_ERprolif$X3.Gene.classifier.subtype <- factor(survival_metabric_X3_ERprolif$X3.Gene.classifier.subtype, levels = c('ER+/HER2- High Prolif', 'ER+/HER2- Low Prolif'))
table(survival_metabric_X3_ERprolif$X3.Gene.classifier.subtype)

# OS
Fit_sig_metabric_os_THR70_X3_ERprolif_THR70_I20 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ X3.Gene.classifier.subtype, data = survival_metabric_X3_ERprolif)

##########
# RFS
Fit_sig_metabric_RFS_THR70_X3_ERprolif_THR70_I20 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ X3.Gene.classifier.subtype, data = survival_metabric_X3_ERprolif)

# RFS COXPH

survival_metabric_X3_ERprolif_2 <- survival_metabric_X3_ERprolif
survival_metabric_X3_ERprolif_2$X3.Gene.classifier.subtype <- factor(survival_metabric_X3_ERprolif_2$X3.Gene.classifier.subtype, levels = c('ER+/HER2- Low Prolif', 'ER+/HER2- High Prolif'))
table(survival_metabric_X3_ERprolif_2$X3.Gene.classifier.subtype)
Fit_sig_metabric_RFS_THR70_X3tnbcHer2_THR70_I20_coxph <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ X3.Gene.classifier.subtype, data = survival_metabric_X3_ERprolif_2)
summary(Fit_sig_metabric_RFS_THR70_X3tnbcHer2_THR70_I20_coxph)


# RFS by Her2
#Fit_sig_metabric_RFS_THR70_ERpos_THR70_I20_byHer2 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_THR70_I20_Merged + HER2.Status, data = survival_metabric_ERpos)


# plot OS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_os_T1_THR70_I20_X3ERprolif.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR70_X3_ERprolif_THR70_I20,
           risk.table = FALSE,
           pval = FALSE,
           pval.size = 12,
           xlim = c(0,240),
           break.x.by = 40,
           legend.labs = levels(survival_metabric_X3_ERprolif$X3.Gene.classifier.subtype),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           #palette = c('#1B9E77', "#c1b026", "#D95F02"),
           risk.table.y.text = FALSE, 
           #title = 'OS in METABRIC in T1 class derived from THR70: THR70 I20 model'
)
dev.off()

#####################
# plot RFS
tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_rfs_T1_THR70_I20_X3ERprolif.tiff", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_sig_metabric_RFS_THR70_X3_ERprolif_THR70_I20,
           risk.table = FALSE,
           pval = FALSE,
           pval.size  = 12,
           xlim = c(0, 240),
           break.x.by = 40,
           legend.labs = c('ER+ HP', 'ER+ LP'),
           legend.title = c(''),
           ggtheme = theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_line(colour = "grey90"),
                           panel.grid.minor = element_line(colour = "grey90"),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.spacing.x = unit(0.5, "cm"),
                           legend.spacing.y = unit(0.5, "cm"),
                           legend.key.height = unit(1.3, "lines"),
                           axis.title = element_text(size = 14, face = 'bold.italic', color = 'black'),
                           axis.text = element_text(size = 12, face = 'bold.italic', color = 'black'), 
                           legend.text = element_text(size = 16, face = 'bold.italic', color = 'black'),
           ), 
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           #title = 'RFS in METABRIC in T1 class derived from THR70: THR70 I20 model'
            ) + guides(
                colour = guide_legend(ncol =1))
dev.off()


#####################
# plot RFS by Her2
# tiff("./figures/T1_DE_THR70_RFS/THR70_metabric_RFS_T1_THR70_I20_20yrs_byHer2.tiff", width = 3200, height = 2200, res = 300)
# ggsurvplot(Fit_sig_metabric_RFS_THR70_T1_THR70_I20_byHer2,
#            risk.table = FALSE,
#            pval = TRUE,
#            pval.size  = 12,
#            xlim = c(0, 240),
#            break.x.by = 80,
#            legend.labs = c('PNBC-A', 'PNBC-B'),
#            legend.title = c('ER-negative clusters'),
#            facet.by = 'HER2.Status', 
#            short.panel.labs	= TRUE,
#            ggtheme = theme_survminer(base_size = 25, font.x = c(25, 'bold.italic', 'black'), font.y = c(25, 'bold.italic', 'black'), font.tickslab = c(25, 'plain', 'black'), font.legend = c(25, 'bold', 'black')),
#            risk.table.y.text.col = FALSE,
#            risk.table.y.text = FALSE, 
#            palette = c('#6057cc', '#8a899a'),
#            #title = 'RFS in METABRIC in T1 class derived from THR70: THR70 I20 model'
# )
# dev.off()

################################################################
## Survival curves: THR50 I20
################################################################

# OS
Fit_metabric_os_THR50_I20 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_THR50_I20_Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_THR50_I20 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_THR50_I20_Merged, data = survival_metabric)

pdf("./figures/T1_DE_THR70_RFS/metabric_os_5clusters_THR50_I20_merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_THR50_I20,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('T1_a', 'Ta_b', 'E3', 'E1', 'E2', 'E4'),
           legend.title	= 'THR70 clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR70 clusters and OS: THR50 + I20')
dev.off()

## RFS: 
pdf("./figures/T1_DE_THR70_RFS/metabric_RFS_5clusters_THR50_I20_merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_THR50_I20,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('T1_a', 'Ta_b', 'E3', 'E1', 'E2', 'E4'),
           legend.title	= 'THR70 clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR70 clusters and RFS: THR50 + I20')
dev.off()

################################################################
## Survival curves: CTLA4
################################################################

# OS
Fit_metabric_os_CTLA4 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_CTLA4Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_CTLA4 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_CTLA4Merged, data = survival_metabric)

pdf("./figures/T1_DE_THR70_RFS/metabric_os_5clusters_CTLA4Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_CTLA4,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('T1_a', 'Ta_b', 'E3', 'E1', 'E2', 'E4'),
           legend.title	= 'THR70 clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR70 clusters and OS: THR70 + CTLA4')
dev.off()

## RFS: 
pdf("./figures/T1_DE_THR70_RFS/metabric_rfs_5clusters_CTLA4Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_CTLA4,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('T1_a', 'Ta_b', 'E3', 'E1', 'E2', 'E4'),
           legend.title	= 'THR70 clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR70 clusters and RFS: THR70 + CTLA4')
dev.off()


################################################################
## Survival curves: CTLA4_ICOS_CXCL13
################################################################

# OS
Fit_metabric_os_CTLA4_ICOS_CXCL13 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_CTLA4_ICOS_CXCL13Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_CTLA4_ICOS_CXCL13 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_CTLA4_ICOS_CXCL13Merged, data = survival_metabric)

pdf("./figures/T1_DE_THR70_RFS/metabric_os_5clusters_CTLA4_ICOS_CXCL13Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_CTLA4_ICOS_CXCL13,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('T1_a', 'Ta_b', 'E3', 'E1', 'E2', 'E4'),
           legend.title	= 'THR70 clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR70 clusters and OS: THR70 + CTLA4 + ICOS + CXCL13')
dev.off()

## RFS: 
pdf("./figures/T1_DE_THR70_RFS/metabric_rfs_5clusters_CTLA4_ICOS_CXCL13Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_CTLA4_ICOS_CXCL13,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('T1_a', 'Ta_b', 'E3', 'E1', 'E2', 'E4'),
           legend.title	= 'THR70 clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR70 clusters and RFS: THR70 + CTLA4 + ICOS + CXCL13')
dev.off()



################################################################
## Survival curves: KLF7
################################################################

# OS
Fit_metabric_os_KLF7 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_KLF7Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_KLF7 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_KLF7Merged, data = survival_metabric)

pdf("./figures/T1_DE_THR70_RFS/metabric_os_5clusters_KLF7Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_KLF7,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('T1_a', 'Ta_b', 'E3', 'E1', 'E2', 'E4'),
           legend.title	= 'THR70 clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR70 clusters and OS: THR70 + KLF7')
dev.off()

## RFS: 
pdf("./figures/T1_DE_THR70_RFS/metabric_rfs_5clusters_KLF7Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_KLF7,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('T1_a', 'Ta_b', 'E3', 'E1', 'E2', 'E4'),
           legend.title	= 'THR70 clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR70 clusters and RFS: THR70 + KLF7')
dev.off()


################################################################
## Survival curves: ZNF627
################################################################

# OS
Fit_metabric_os_ZNF627 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_ZNF627Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_ZNF627 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_ZNF627Merged, data = survival_metabric)

pdf("./figures/T1_DE_THR70_RFS/metabric_os_5clusters_ZNF627Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_ZNF627,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('T1_a', 'Ta_b', 'E3', 'E1', 'E2', 'E4'),
           legend.title	= 'THR70 clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR70 clusters and OS: THR70 + ZNF627')
dev.off()

## RFS: 
pdf("./figures/T1_DE_THR70_RFS/metabric_rfs_5clusters_ZNF627Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_ZNF627,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('T1_a', 'Ta_b', 'E3', 'E1', 'E2', 'E4'),
           legend.title	= 'THR70 clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR70 clusters and RFS: THR70 + ZNF627')
dev.off()




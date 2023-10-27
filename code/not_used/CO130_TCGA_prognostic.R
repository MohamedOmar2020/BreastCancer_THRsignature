
############################################################################
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
library(vtable)
library(readxl)

################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')

rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'DULLARD'] <- 'HSA011916'
rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'GPR56'] <- 'ADGRG1'
rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'FAM116B'] <- 'DENND6B'
rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'FAM46B'] <- 'TENT5B'


#######################################
## read TCGA PAM50
Pam50_tcga <- read.delim('data/brca_tcga/data_clinical_patient.txt') 
Pam50_tcga <- Pam50_tcga[-c(1:4), ]
Pam50_tcga$Patient.ID <- gsub("\\-", "\\.", Pam50_tcga$X.Patient.Identifier)
Pam50_tcga <- Pam50_tcga[, c('Patient.ID', 'Subtype', 'Overall.Survival.Status', 'Overall.Survival..Months.', 'Disease.Free.Status', 'Disease.Free..Months.', 'Progression.Free.Status', 'Progress.Free.Survival..Months.')]

# merge with TCGA pheno table
summary(Pheno_tcga$Patient.ID %in% Pam50_tcga$Patient.ID)

# remove the old survival info 
Pheno_tcga <- Pheno_tcga[, !colnames(Pheno_tcga) %in% c('Overall.Survival..Months.', 'Overall.Survival.Status', 'Disease.Free.Status', 'Disease.Free..Months.')]

# merge
Pheno_tcga <- merge(Pheno_tcga, Pam50_tcga, by = "Patient.ID")
rownames(Pheno_tcga) <- Pheno_tcga$Patient.ID

table(Pheno_tcga$Subtype)
table(Pam50_tcga$Subtype)

##########################
# remove weird cancer types
##########################
table(Pheno_metabric$Cancer.Type.Detailed)
table(Pheno_tcga$Cancer.Type.Detailed)


Pheno_tcga <- Pheno_tcga[!(Pheno_tcga$Cancer.Type.Detailed %in% c('Adenoid Cystic Breast Cancer', 
                                                                  'Basal Cell Carcinoma', 
                                                                  'Malignant Phyllodes Tumor of the Breast', 
                                                                  'Paget Disease of the Nipple', 
                                                                  'Solid Papillary Carcinoma of the Breast')), 
]

Expr_tcga_refAll <- Expr_tcga_refAll[, colnames(Expr_tcga_refAll) %in% rownames(Pheno_tcga)]

##############################################
### combine in 1 dataset
##############################################
group_tcga <- as.factor(Pheno_tcga$Overall.Survival.Status)
table(group_tcga)
Data_tcga <- as.data.frame(cbind(t(Expr_tcga_refAll), group_tcga))
Data_tcga$group_tcga <- as.factor(Data_tcga$group_tcga)
levels(Data_tcga$group_tcga) <- c('0', '1')
table(Data_tcga$group_tcga)
colnames(Data_tcga)[colnames(Data_tcga) %in% c('group_tcga')] <- c('os')
table(Data_tcga$os)

group_tcga <- Data_tcga$os

#############################################################################################################
##############################################################################################################
# load the trained the model
load('objs/CO130_model_logreg.rda')


###########################################################################
############################################################################
### predict in the training dataset
# Make predictions

Train_prob_CO130 <- CO130_model %>% predict(Data_tcga , type = "response")

### Threshold
thr_CO130 <- coords(roc(group_tcga, Train_prob_CO130, levels = c("0", "1"), direction = "<"), "best")["threshold"]
thr_CO130

### ROC Curve
ROCTrain_CO130 <- roc(group_tcga, Train_prob_CO130, plot = F, print.thres=thr_CO130$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("0", "1"), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_CO130

### Get predictions based on best threshold from ROC curve
Train_predClasses_CO130 <- ifelse(Train_prob_CO130 >= thr_CO130$threshold, "1", "0")
table(Train_predClasses_CO130)
Train_predClasses_CO130 <- factor(Train_predClasses_CO130, levels = c('0', '1'))


### Resubstitution performance in the TRAINING set
ConfusionTrain_CO130 <- confusionMatrix(Train_predClasses_CO130, group_tcga, positive = "1", mode = "everything")
ConfusionTrain_CO130

## MCC
MCC_Train_CO130 <- mltools::mcc(pred = Train_predClasses_CO130, actuals = group_tcga)
MCC_Train_CO130

##########################
## Keep only the relevant information (Metastasis Event and Time)
Phenotype_tcga <- cbind(Pheno_tcga[, c("Overall.Survival.Status", "Overall.Survival..Months.", 
                                       "Disease.Free.Status", "Disease.Free..Months.", 
                                       "Subtype", 
                                       "ER.Status.By.IHC", 
                                       "PR.status.by.ihc", "HER2.fish.status", "IHC.HER2")], 
                        Train_prob_CO130, Train_predClasses_CO130)

# create a merged pdata and Z-scores object
CoxData_tcga <- data.frame(Phenotype_tcga)

# divide the probabilities into quartiles
CoxData_tcga <- CoxData_tcga %>%
  mutate(tcga_prob_CO130_quartiles = ntile(Train_prob_CO130, 4), 
         tcga_prob_CO130_quintiles = ntile(Train_prob_CO130, 5),
         tcga_prob_CO130_tertiles = ntile(Train_prob_CO130, 3)
  )


#############################################################################################
# convert time to numeric
CoxData_tcga$Overall.Survival..Months. <- as.numeric(CoxData_tcga$Overall.Survival..Months.)
CoxData_tcga$Disease.Free..Months. <- as.numeric(CoxData_tcga$Disease.Free..Months.)

# fix survival info
CoxData_tcga$Disease.Free.Status <- gsub("\\:.+", "", CoxData_tcga$Disease.Free.Status)
CoxData_tcga$Overall.Survival.Status <- gsub("\\:.+", "", CoxData_tcga$Overall.Survival.Status)

CoxData_tcga$Disease.Free.Status <- as.numeric(CoxData_tcga$Disease.Free.Status)
CoxData_tcga$Overall.Survival.Status <- as.numeric(CoxData_tcga$Overall.Survival.Status)

########################################################################  
## Fit survival curves
########################################################################  

# OS
## tcga all genes
Fit_sig_tcga_os <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_CO130, data = CoxData_tcga)


## by quartiles
Fit_sig_tcga_os_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_quartiles, data = CoxData_tcga)

## by quintiles
Fit_sig_tcga_os_quintiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_quintiles, data = CoxData_tcga)

## by tertiles
Fit_sig_tcga_os_tertiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_tertiles, data = CoxData_tcga)

#############
# DFS
## tcga all genes
Fit_sig_tcga_DFS <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ Train_predClasses_CO130, data = CoxData_tcga)

## by quartiles
Fit_sig_tcga_DFS_quartiles <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_quartiles, data = CoxData_tcga)

## by quintiles
Fit_sig_tcga_DFS_quintiles <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_quintiles, data = CoxData_tcga)

## by tertiles
Fit_sig_tcga_DFS_tertiles <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_tertiles, data = CoxData_tcga)

#################################
## by clinical groups

## pam50

# keep only the major pam50 subtypes
CoxData_tcga_PAM <- CoxData_tcga %>%
  filter(Subtype %in% c('BRCA_Basal', 'BRCA_Her2', 'BRCA_LumA', 'BRCA_LumB', 'BRCA_Normal'))

# keep only quartiles 1 and 4 (has to be in each of THR25, CO130_1, and CO130_2)

CoxData_tcga_PAM_Q1vsQ4 <- CoxData_tcga_PAM %>%
  filter(tcga_prob_CO130_quartiles %in% c('1', '4'))

CoxData_tcga_PAM_Q1vsQ5 <- CoxData_tcga_PAM %>%
  filter(tcga_prob_CO130_quintiles %in% c('1', '5'))


# os: all
Fit_sig_tcga_os_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_CO130 + Subtype, data = CoxData_tcga_PAM)

# os: quartiles: all
Fit_sig_tcga_os_quartiles_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_quartiles + Subtype, data = CoxData_tcga_PAM)

# os: quintiles: all
Fit_sig_tcga_os_quintiles_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_quintiles + Subtype, data = CoxData_tcga_PAM)

# os: quartiles: Q1 vs Q4
Fit_sig_tcga_os_Q1vsQ4_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_quartiles + Subtype, data = CoxData_tcga_PAM_Q1vsQ4)

# os: quintiles: Q1 vs Q5
Fit_sig_tcga_os_Q1vsQ5_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_quintiles + Subtype, data = CoxData_tcga_PAM_Q1vsQ5)

##########
# DFS: all
Fit_sig_tcga_dfs_PAM <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ Train_predClasses_CO130 + Subtype, data = CoxData_tcga_PAM)

# rfs: quartiles: all
Fit_sig_tcga_dfs_quartiles_PAM <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_quartiles + Subtype, data = CoxData_tcga_PAM)

# rfs: quintiles: all
Fit_sig_tcga_dfs_quintiles_PAM <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_quintiles + Subtype, data = CoxData_tcga_PAM)

# rfs: quartiles: Q1 vs Q4
Fit_sig_tcga_dfs_Q1vsQ4_PAM <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_quartiles + Subtype, data = CoxData_tcga_PAM_Q1vsQ4)

# rfs: quintiles: Q1 vs Q5
Fit_sig_tcga_dfs_Q1vsQ5_PAM <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_quintiles + Subtype, data = CoxData_tcga_PAM_Q1vsQ5)

#############
## ER

# keep only quartiles 1 and 4 (has to be in each of THR25, CO130_1, and CO130_2)

# Note: we can use this for both ER and X3
CoxData_tcga_Q1vsQ4 <- CoxData_tcga %>%
  filter(tcga_prob_CO130_quartiles %in% c('1', '4'))

CoxData_tcga_Q1vsQ5 <- CoxData_tcga %>%
  filter(tcga_prob_CO130_quintiles %in% c('1', '5'))

# os: all
Fit_sig_tcga_os_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_CO130 + ER.Status.By.IHC, data = CoxData_tcga)

# os: quartiles: all
Fit_sig_tcga_os_quartiles_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_quartiles + ER.Status.By.IHC, data = CoxData_tcga)

# os: quintiles: all
Fit_sig_tcga_os_quintiles_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_quintiles + ER.Status.By.IHC, data = CoxData_tcga)

# os: quartiles: Q1 vs Q4
Fit_sig_tcga_os_Q1vsQ4_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_quartiles + ER.Status.By.IHC, data = CoxData_tcga_Q1vsQ4)

# os: quintiles: Q1 vs Q5
Fit_sig_tcga_os_Q1vsQ5_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_quintiles + ER.Status.By.IHC, data = CoxData_tcga_Q1vsQ5)

##############
# dfs: all
Fit_sig_tcga_dfs_ER <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ Train_predClasses_CO130 + ER.Status.By.IHC, data = CoxData_tcga)

# rfs: quartiles: all
Fit_sig_tcga_dfs_quartiles_ER <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_quartiles + ER.Status.By.IHC, data = CoxData_tcga)

# rfs: quintiles: all
Fit_sig_tcga_dfs_quintiles_ER <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_quintiles + ER.Status.By.IHC, data = CoxData_tcga)

# rfs: quartiles: Q1 vs Q4
Fit_sig_tcga_dfs_Q1vsQ4_ER <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_quartiles + ER.Status.By.IHC, data = CoxData_tcga_Q1vsQ4)

# rfs: quintiles: Q1 vs Q5
Fit_sig_tcga_dfs_Q1vsQ5_ER <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_quintiles + ER.Status.By.IHC, data = CoxData_tcga_Q1vsQ5)


############################################################################
# plot OS
############################################################################

tiff("./figures/CO130_tcga/CO130_tcga_os_all.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_tcga_os,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           #title = 'THR 50 and tcga OS'
)
dev.off()

########
# by quartiles
tiff("./figures/CO130_tcga/CO130_tcga_os_quartiles_10yrs.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_tcga_os_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           xlim = c(0,120),
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and tcga OS: quartiles'
)
dev.off()

# Fit the Cox proportional hazards model
os_surv_obj <- Surv(CoxData_tcga$`Overall.Survival..Months.`, CoxData_tcga$`Overall.Survival.Status`)
cox_Fit_sig_tcga_os_CO130_quartiles <- coxph(os_surv_obj ~ as.factor(tcga_prob_CO130_quartiles), data = CoxData_tcga)
summary(cox_Fit_sig_tcga_os_CO130_quartiles)

#######
# merge Q2-Q3
CoxData_tcga$tcga_prob_CO130_quartiles_merged <- as.factor(CoxData_tcga$tcga_prob_CO130_quartiles)

levels(CoxData_tcga$tcga_prob_CO130_quartiles_merged) <- c('Q1', 'Q2-3', 'Q2-3', 'Q4')
table(CoxData_tcga$tcga_prob_CO130_quartiles_merged)

Fit_sig_tcga_os_CO130_quartiles_merged <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_quartiles_merged, data = CoxData_tcga)
tiff("./figures/CO130_tcga/CO130_tcga_os_quartiles_merged_10yrs.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_tcga_os_CO130_quartiles_merged,
           risk.table = TRUE,
           pval = TRUE,
           pval.size = 12,
           xlim = c(0,120),
           legend.labs = c('Q1', 'Q2-Q3', 'Q4'),
           ggtheme = theme_survminer(base_size = 20, font.x = c(20, 'bold.italic', 'black'), font.y = c(20, 'bold.italic', 'black'), font.tickslab = c(20, 'plain', 'black'), font.legend = c(20, 'bold', 'black')),
           palette = 'jco',
           #risk.table.y.text.col = FALSE,
           #risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles'
)
dev.off()

# Fit the Cox proportional hazards model
os_surv_obj <- Surv(CoxData_tcga$`Overall.Survival..Months.`, CoxData_tcga$`Overall.Survival.Status`)
cox_Fit_sig_tcga_os_CO130_quartiles_merged <- coxph(os_surv_obj ~ as.factor(tcga_prob_CO130_quartiles_merged), data = CoxData_tcga)
summary(cox_Fit_sig_tcga_os_CO130_quartiles_merged)

#############
# by quintiles
tiff("./figures/CO130_tcga/CO130_tcga_os_quintiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_tcga_os_quintiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and tcga OS: quartiles'
)
dev.off()

#############
# by tertiles
tiff("./figures/CO130_tcga/CO130_tcga_os_tertiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_tcga_os_tertiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('Q1', 'Q2', 'Q3'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and tcga OS: quartiles'
)
dev.off()

######################################
# plot RFS
######################################
tiff("./figures/CO130_tcga/CO130_tcga_DFS_allpairs.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_tcga_DFS,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           #title = 'THR 50_1 (logistic regression) and tcga RFS'
)
dev.off()


########
# by quartiles

tiff("./figures/CO130_tcga/CO130_tcga_DFS_quartiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_tcga_DFS_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           #risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and tcga RFS: quartiles'
)
dev.off()

########
# by quintiles

tiff("./figures/CO130_tcga/CO130_tcga_DFS_quintiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_tcga_DFS_quintiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and tcga OS: quartiles'
)
dev.off()

########
# by tertiles

tiff("./figures/CO130_tcga/CO130_tcga_DFS_tertiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_tcga_DFS_tertiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('Q1', 'Q2', 'Q3'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and tcga OS: quartiles'
)
dev.off()
############################################################################
############################################################################
### by clinical group

## PAM50

# OS

pdf("./figures/CO130_prognostic/tcga/CO130_tcga_os_PAM50.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_os_CO130_PAM,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           short.panel.labs = T,
           facet.by = "Subtype",
           ggtheme = theme_minimal(),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'CO130 and tcga OS by PAM50 subtypes')
dev.off()


#########
# RFS

pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_RFS_PAM50.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_dfs_CO130_PAM,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           short.panel.labs = T,
           facet.by = "Subtype",
           ggtheme = theme_minimal(),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'CO130 and tcga DFS by PAM50 subtypes')
dev.off()


######################################################
# OS: quartiles: all


pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_os_PAM50_quartiles.pdf", width = 12, height = 10, onefile = F)
ggsurvplot(Fit_sig_tcga_os_CO130_quartiles_PAM,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           short.panel.labs = T,
           facet.by = "Subtype",
           ggtheme = theme_survminer(base_size = 25, font.x = c(25, 'bold.italic', 'black'), font.y = c(25, 'bold.italic', 'black'), font.tickslab = c(25, 'plain', 'black'), font.legend = c(25, 'bold', 'black')),
           palette = 'jco',
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'CO130 and tcga OS by PAM50 subtypes: quartiles')
dev.off()

####################################################
## OS: quintiles: all

tiff("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_os_PAM50_quintiles.tiff", width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_tcga_os_CO130_quintiles_PAM,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "Subtype",
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
           legend.title	= 'Quintiles',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and tcga OS: quartiles'
)
dev.off()

######################################################
# OS: quartiles: Q1 vs Q4

tiff("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_os_PAM50_Q1vsQ4.tiff",  width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_tcga_os_CO130_Q1vsQ4_PAM,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "Subtype",
           legend.title	= 'Quartiles',
           pval.size = 15,
           #break.x.by = 20,
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           legend.labs = c('Q1', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and tcga OS by PAM50 subtypes: Q1 vs Q4'
)
dev.off()

######################################################
# OS: quintiles: Q1 vs Q5
tiff("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_os_PAM50_Q1vsQ5.tiff",  width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_tcga_os_CO130_Q1vsQ5_PAM,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "Subtype",
           legend.labs = c('Q1', 'Q5'),
           legend.title	= 'Quintiles',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and tcga OS: quartiles'
)
dev.off()

#############################
# RFS: quartiles: all

tiff("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_dfs_PAM50_quartiles.tiff",  width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_tcga_dfs_CO130_quartiles_PAM,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "Subtype",
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           legend.title	= 'Quartiles',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text = FALSE, title = 'CO130 and tcga DFS by PAM50 subtypes: quartiles')
dev.off()


####################################################
## RFS: quintiles: all

tiff("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_dfs_PAM50_quintiles.tiff",  width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_tcga_dfs_CO130_quintiles_PAM,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "Subtype",
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
           legend.title	= 'Quintiles',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and tcga OS: quartiles'
)
dev.off()

#############################
# RFS: quartiles: Q1 vs Q4
tiff("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_dfs_PAM50_Q1vsQ4.tiff",  width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_tcga_dfs_CO130_Q1vsQ4_PAM,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           short.panel.labs = T,
           facet.by = "Subtype",
           #legend.labs = c('Basal', 'Claudin-low', 'Her2+', 'Luminal A', 'Luminal B'),
           legend.title	= 'Quartiles',
           #break.x.by = 20,
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q4'),
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and tcga RFS by PAM50 subtypes: Q1 vs Q4'
)
dev.off()

#############################
# RFS: quintiles: Q1 vs Q5
tiff("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_dfs_PAM50_Q1vsQ5.tiff",  width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_tcga_dfs_CO130_Q1vsQ5_PAM,
           risk.table = FALSE,
           pval = T,
           short.panel.labs = T,
           facet.by = "Subtype",
           #legend.labs = c('Basal', 'Claudin-low', 'Her2+', 'Luminal A', 'Luminal B'),
           legend.title	= 'Quintiles',
           pval.size = 12,
           #break.x.by = 20,
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q5'),
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and tcga RFS by PAM50 subtypes: Q1 vs Q4'
)
dev.off()
##############################################################################################
##############################################################################################
## ER

# OS

pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_os_ER.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_os_CO130_ER,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "ER.Status.By.IHC",
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and tcga OS by ER status'
)
dev.off()


######################
# RFS

pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_DFS_ER.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_dfs_CO130_ER,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "ER.Status.By.IHC",
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and tcga OS by ER status'
)
dev.off()



####################################################
# OS: quartiles: all

pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_os_ER_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_os_CO130_quartiles_ER,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "ER.Status.By.IHC",
           ggtheme = theme_minimal(),
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and tcga OS by ER status: quartiles'
)
dev.off()



####################################################
# OS: quartiles: Q1 vs Q4
pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_os_ER_Q1vsQ4.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_os_CO130_Q1vsQ4_ER,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "ER.Status.By.IHC",
           ggtheme = theme_minimal(),
           legend.labs = c('Q1', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = THR 50_1 and tcga OS by ER status: Q1 vs Q4'
)
dev.off()

####################################################
# RFS: quartiles: all
pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_rfs_ER_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_rfs_CO130_quartiles_ER,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "ER.Status.By.IHC",
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and tcga RFS by ER status: quartiles'
)
dev.off()

####################################################
# RFS: quartiles: Q1 vs Q4
pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_rfs_ER_Q1vsQ4.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_rfs_CO130_Q1vsQ4_ER,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "ER.Status.By.IHC",
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q4'),
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and tcga RFS by ER status: Q1 vs Q4'
)
dev.off()

##############################################################################################
##############################################################################################
## X3

# OS
pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_os_X3.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_os_CO130_X3,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and tcga OS by X3 classifier subtypes'
)
dev.off()


######################
# RFS
pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_rfs_X3.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_rfs_CO130_X3,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and tcga RFS by X3 classifier subtypes'
)
dev.off()


####################################################
# OS: quartiles: all

pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_os_X3_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_os_CO130_quartiles_X3,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           ggtheme = theme_minimal(),
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and tcga OS by X3 classifier subtypes: quartiles'
)
dev.off()

####################################################
## OS: quintiles: all
pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_os_X3_quintiles.tiff", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_os_CO130_quintiles_X3,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
           legend.title	= 'Quintiles',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and tcga OS: quartiles'
)
dev.off()

####################################################
# OS: quartiles: Q1 vs Q4
pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_os_X3_Q1vsQ4.tiff", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_os_CO130_Q1vsQ4_X3,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           legend.title	= 'Quartiles',
           pval.size = 15,
           #break.x.by = 20,
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           legend.labs = c('Q1', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and tcga OS by PAM50 subtypes: Q1 vs Q4'
)
dev.off()

####################################################
## OS: quintiles: Q1 vs Q5
pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_os_X3_Q1vsQ5.tiff", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_os_CO130_Q1vsQ5_X3,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           legend.labs = c('Q1', 'Q5'),
           legend.title	= 'Quintiles',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and tcga OS: quartiles'
)
dev.off()


###############################################################
# RFS: quartiles: all
pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_rfs_X3_quartiles.tiff", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_rfs_CO130_quartiles_X3,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           legend.title	= 'Quartiles',
           pval.size = 15,
           #break.x.by = 20,
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and tcga RFS by X3 classifier subtypes: quartiles'
)
dev.off()

####################################################
## RFS: quintiles: all
pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_rfs_X3_quintiles.tiff", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_rfs_CO130_quintiles_X3,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
           legend.title	= 'Quintiles',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and tcga OS: quartiles'
)
dev.off()


###############################################################
# RFS: quartiles: Q1 vs Q4

pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_rfs_X3_Q1vsQ4.tiff", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_rfs_CO130_Q1vsQ4_X3,
           risk.table = FALSE,
           pval = F,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           legend.title	= 'Quartiles',
           pval.size = 15,
           #break.x.by = 20,
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q4'),
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and tcga RFS by X3 classifier subtypes: Q1 vs Q4'
)
dev.off()


####################################################
## RFS: quintiles: Q1 vs Q5
pdf("./figures/CO130_tcga/logreg/byClinicalGroup/CO130_tcga_rfs_X3_Q1vsQ5.tiff", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_rfs_CO130_Q1vsQ5_X3,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           legend.labs = c('Q1', 'Q5'),
           legend.title	= 'Quintiles',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and tcga OS: quartiles'
)
dev.off()


##############################################################################
##############################################################################
##############################################################################
## fit coxph model:

########
## OS

# by probaility

#Fit_sig_tcga_os_coxph_THR25 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_prob_THR25, data = CoxData_tcga)
#summary(Fit_sig_tcga_os_coxph_THR25)

Fit_sig_tcga_os_coxph_CO130_1 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_prob_CO130_1, data = CoxData_tcga)
summary(Fit_sig_tcga_os_coxph_CO130_1)

#Fit_sig_tcga_os_coxph_CO130_2 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_prob_CO130_2, data = CoxData_tcga)
#summary(Fit_sig_tcga_os_coxph_CO130_2)

#png('./figures/logreg/THR25_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
#ggforest(Fit_sig_tcga_coxph_THR25, fontsize = 0.5)
#dev.off()

#png('./figures/logreg/CO130_1_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
#ggforest(Fit_sig_tcga_coxph_CO130_1, fontsize = 0.5)
#dev.off()

#png('./figures/logreg/CO130_2_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
#ggforest(Fit_sig_tcga_coxph_CO130_2, fontsize = 0.5)
#dev.off()

########
## by quartiles

# make a factor with Q1 (lowest risk) being the reference
#CoxData_tcga$tcga_prob_THR25_quartiles <- factor(CoxData_tcga$tcga_prob_THR25_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_tcga$tcga_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_tcga$tcga_prob_THR25_quartiles))

CoxData_tcga$tcga_prob_CO130_1_quartiles <- factor(CoxData_tcga$tcga_prob_CO130_1_quartiles, levels = c('1', '2', '3', '4'))
levels(CoxData_tcga$tcga_prob_CO130_1_quartiles) <- paste0('Q', levels(CoxData_tcga$tcga_prob_CO130_1_quartiles))

#CoxData_tcga$tcga_prob_CO130_2_quartiles <- factor(CoxData_tcga$tcga_prob_CO130_2_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_tcga$tcga_prob_CO130_2_quartiles) <- paste0('Q', levels(CoxData_tcga$tcga_prob_CO130_2_quartiles))

# fit
#Fit_sig_tcga_os_coxph_THR25_quartiles <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_THR25_quartiles, data = CoxData_tcga)
#summary_tcga_os_coxph_THR25_quartiles <- summary(Fit_sig_tcga_os_coxph_THR25_quartiles)

Fit_sig_tcga_os_coxph_CO130_1_quartiles <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_1_quartiles, data = CoxData_tcga)
summary_tcga_os_coxph_CO130_1_quartiles <- summary(Fit_sig_tcga_os_coxph_CO130_1_quartiles)

#Fit_sig_tcga_os_coxph_CO130_2_quartiles <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_2_quartiles, data = CoxData_tcga)
#summary_tcga_os_coxph_CO130_2_quartiles <- summary(Fit_sig_tcga_os_coxph_CO130_2_quartiles)

#summary_list_tcga_os_quartiles <- list(THR25 = summary_tcga_os_coxph_THR25_quartiles, CO130_1 = summary_tcga_os_coxph_CO130_1_quartiles, CO130_2 = summary_tcga_os_coxph_CO130_2_quartiles) 

# get the HR
# HR_list_tcga_os_coxph_quartiles <- lapply(summary_list_tcga_os_quartiles, function(x){
#   HR <- x$conf.int[, 'exp(coef)']
#   lower_95CI <- x$conf.int[, 'lower .95']
#   upper_95CI <- x$conf.int[, 'upper .95']
#   Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
#   Pvalue_logrank_test <- x$sctest['pvalue']
#   Pvalue_wald_test <- x$waldtest['pvalue']
#   data.frame(HR = HR, lower_95CI = lower_95CI, upper_95CI = upper_95CI, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test,
#              Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
# })

# HR_df_tcga_os_coxph_quartiles_THR25 <- as.data.frame(HR_list_tcga_os_coxph_quartiles$THR25)
# HR_df_tcga_os_coxph_quartiles_THR25$quartile <- gsub('.+quartiles', '', rownames(HR_df_tcga_os_coxph_quartiles_THR25))
# 
# HR_df_tcga_os_coxph_quartiles_CO130_1 <- as.data.frame(HR_list_tcga_os_coxph_quartiles$CO130_1)
# HR_df_tcga_os_coxph_quartiles_CO130_1$quartile <- gsub('.+quartiles', '', rownames(HR_df_tcga_os_coxph_quartiles_CO130_1))
# 
# HR_df_tcga_os_coxph_quartiles_CO130_2 <- as.data.frame(HR_list_tcga_os_coxph_quartiles$CO130_2)
# HR_df_tcga_os_coxph_quartiles_CO130_2$quartile <- gsub('.+quartiles', '', rownames(HR_df_tcga_os_coxph_quartiles_CO130_2))
# 
# # save the results
# write.csv(HR_df_tcga_os_coxph_quartiles_THR25, 'objs/HR/tcga/OS/THR25_quartiles_HR.csv')
# write.csv(HR_df_tcga_os_coxph_quartiles_CO130_1, 'objs/HR/tcga/OS/CO130_1_quartiles_HR.csv')
# write.csv(HR_df_tcga_os_coxph_quartiles_CO130_2, 'objs/HR/tcga/OS/CO130_2_quartiles_HR.csv')

################
# RFS

# by probaility

#Fit_sig_tcga_RFS_coxph_THR25 <- coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ Train_prob_THR25, data = CoxData_tcga)
#summary(Fit_sig_tcga_RFS_coxph_THR25)

Fit_sig_tcga_RFS_coxph_CO130_1 <- coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ Train_prob_CO130_1, data = CoxData_tcga)
summary(Fit_sig_tcga_RFS_coxph_CO130_1)

#Fit_sig_tcga_RFS_coxph_CO130_2 <- coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ Train_prob_CO130_2, data = CoxData_tcga)
#summary(Fit_sig_tcga_RFS_coxph_CO130_2)


########
## by quartiles

# fit
#Fit_sig_tcga_RFS_coxph_THR25_quartiles <- coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_THR25_quartiles, data = CoxData_tcga)
#summary_tcga_RFS_coxph_THR25_quartiles <- summary(Fit_sig_tcga_RFS_coxph_THR25_quartiles)

Fit_sig_tcga_RFS_coxph_CO130_1_quartiles <- coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_1_quartiles, data = CoxData_tcga)
summary_tcga_RFS_coxph_CO130_1_quartiles <- summary(Fit_sig_tcga_RFS_coxph_CO130_1_quartiles)

#Fit_sig_tcga_RFS_coxph_CO130_2_quartiles <- coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_2_quartiles, data = CoxData_tcga)
#summary_tcga_RFS_coxph_CO130_2_quartiles <- summary(Fit_sig_tcga_RFS_coxph_CO130_2_quartiles)

#summary_list_tcga_RFS_quartiles <- list(THR25 = summary_tcga_RFS_coxph_THR25_quartiles, CO130_1 = summary_tcga_RFS_coxph_CO130_1_quartiles, CO130_2 = summary_tcga_RFS_coxph_CO130_2_quartiles) 

# get the HR
# HR_list_tcga_RFS_coxph_quartiles <- lapply(summary_list_tcga_RFS_quartiles, function(x){
#   HR <- x$conf.int[, 'exp(coef)']
#   lower_95CI <- x$conf.int[, 'lower .95']
#   upper_95CI <- x$conf.int[, 'upper .95']
#   Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
#   Pvalue_logrank_test <- x$sctest['pvalue']
#   Pvalue_wald_test <- x$waldtest['pvalue']
#   data.frame(HR = HR, lower_95CI = lower_95CI, upper_95CI = upper_95CI, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test,
#              Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
# })
# 
# HR_df_tcga_RFS_coxph_quartiles_THR25 <- as.data.frame(HR_list_tcga_RFS_coxph_quartiles$THR25)
# HR_df_tcga_RFS_coxph_quartiles_THR25$quartile <- gsub('.+quartiles', '', rownames(HR_df_tcga_RFS_coxph_quartiles_THR25))
# 
# HR_df_tcga_RFS_coxph_quartiles_CO130_1 <- as.data.frame(HR_list_tcga_RFS_coxph_quartiles$CO130_1)
# HR_df_tcga_RFS_coxph_quartiles_CO130_1$quartile <- gsub('.+quartiles', '', rownames(HR_df_tcga_RFS_coxph_quartiles_CO130_1))
# 
# HR_df_tcga_RFS_coxph_quartiles_CO130_2 <- as.data.frame(HR_list_tcga_RFS_coxph_quartiles$CO130_2)
# HR_df_tcga_RFS_coxph_quartiles_CO130_2$quartile <- gsub('.+quartiles', '', rownames(HR_df_tcga_RFS_coxph_quartiles_CO130_2))
# 
# # save the results
# write.csv(HR_df_tcga_RFS_coxph_quartiles_THR25, 'objs/HR/tcga/RFS/THR25_quartiles_HR.csv')
# write.csv(HR_df_tcga_RFS_coxph_quartiles_CO130_1, 'objs/HR/tcga/RFS/CO130_1_quartiles_HR.csv')
# write.csv(HR_df_tcga_RFS_coxph_quartiles_CO130_2, 'objs/HR/tcga/RFS/CO130_2_quartiles_HR.csv')


########################################################################  
########################################################################  
## by clinical groups

########
## by quartiles

## make a factor with Q1 (lowest risk) being the reference

# CoxData_tcga
#CoxData_tcga$tcga_prob_THR25_quartiles <- factor(CoxData_tcga$tcga_prob_THR25_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_tcga$tcga_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_tcga$tcga_prob_THR25_quartiles))

#CoxData_tcga$tcga_prob_CO130_1_quartiles <- factor(CoxData_tcga$tcga_prob_CO130_1_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_tcga$tcga_prob_CO130_1_quartiles) <- paste0('Q', levels(CoxData_tcga$tcga_prob_CO130_1_quartiles))

#CoxData_tcga$tcga_prob_CO130_2_quartiles <- factor(CoxData_tcga$tcga_prob_CO130_2_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_tcga$tcga_prob_CO130_2_quartiles) <- paste0('Q', levels(CoxData_tcga$tcga_prob_CO130_2_quartiles))

######
# CoxData_tcga_PAM
#CoxData_tcga_PAM$tcga_prob_THR25_quartiles <- factor(CoxData_tcga_PAM$tcga_prob_THR25_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_tcga_PAM$tcga_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_tcga_PAM$tcga_prob_THR25_quartiles))

CoxData_tcga_PAM$tcga_prob_CO130_1_quartiles <- factor(CoxData_tcga_PAM$tcga_prob_CO130_1_quartiles, levels = c('1', '2', '3', '4'))
levels(CoxData_tcga_PAM$tcga_prob_CO130_1_quartiles) <- paste0('Q', levels(CoxData_tcga_PAM$tcga_prob_CO130_1_quartiles))

#CoxData_tcga_PAM$tcga_prob_CO130_2_quartiles <- factor(CoxData_tcga_PAM$tcga_prob_CO130_2_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_tcga_PAM$tcga_prob_CO130_2_quartiles) <- paste0('Q', levels(CoxData_tcga_PAM$tcga_prob_CO130_2_quartiles))

######
# CoxData_tcga_PAM_Q1vsQ4: 
#CoxData_tcga_PAM_Q1vsQ4_THR25$tcga_prob_THR25_quartiles <- factor(CoxData_tcga_PAM_Q1vsQ4_THR25$tcga_prob_THR25_quartiles, levels = c('1', '4'))
#levels(CoxData_tcga_PAM_Q1vsQ4_THR25$tcga_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_tcga_PAM_Q1vsQ4_THR25$tcga_prob_THR25_quartiles))

CoxData_tcga_PAM_Q1vsQ4_CO130_1$tcga_prob_CO130_1_quartiles <- factor(CoxData_tcga_PAM_Q1vsQ4_CO130_1$tcga_prob_CO130_1_quartiles, levels = c('1', '4'))
levels(CoxData_tcga_PAM_Q1vsQ4_CO130_1$tcga_prob_CO130_1_quartiles) <- paste0('Q', levels(CoxData_tcga_PAM_Q1vsQ4_CO130_1$tcga_prob_CO130_1_quartiles))

#CoxData_tcga_PAM_Q1vsQ4_CO130_2$tcga_prob_CO130_2_quartiles <- factor(CoxData_tcga_PAM_Q1vsQ4_CO130_2$tcga_prob_CO130_2_quartiles, levels = c('1', '4'))
#levels(CoxData_tcga_PAM_Q1vsQ4_CO130_2$tcga_prob_CO130_2_quartiles) <- paste0('Q', levels(CoxData_tcga_PAM_Q1vsQ4_CO130_2$tcga_prob_CO130_2_quartiles))


######
# CoxData_tcga_Q1vsQ4
#CoxData_tcga_Q1vsQ4_THR25$tcga_prob_THR25_quartiles <- factor(CoxData_tcga_Q1vsQ4_THR25$tcga_prob_THR25_quartiles, levels = c('1', '4'))
#levels(CoxData_tcga_Q1vsQ4_THR25$tcga_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_tcga_Q1vsQ4_THR25$tcga_prob_THR25_quartiles))

CoxData_tcga_Q1vsQ4_CO130_1$tcga_prob_CO130_1_quartiles <- factor(CoxData_tcga_Q1vsQ4_CO130_1$tcga_prob_CO130_1_quartiles, levels = c('1', '4'))
levels(CoxData_tcga_Q1vsQ4_CO130_1$tcga_prob_CO130_1_quartiles) <- paste0('Q', levels(CoxData_tcga_Q1vsQ4_CO130_1$tcga_prob_CO130_1_quartiles))

#CoxData_tcga_Q1vsQ4_CO130_2$tcga_prob_CO130_2_quartiles <- factor(CoxData_tcga_Q1vsQ4_CO130_2$tcga_prob_CO130_2_quartiles, levels = c('1', '4'))
#levels(CoxData_tcga_Q1vsQ4_CO130_2$tcga_prob_CO130_2_quartiles) <- paste0('Q', levels(CoxData_tcga_Q1vsQ4_CO130_2$tcga_prob_CO130_2_quartiles))


#CoxData_tcga_PAM$Pam50...Claudin.low.subtype <- factor(CoxData_tcga_PAM$Pam50...Claudin.low.subtype)

##############################################
## fit

## PAM50

# os: quartiles: all
# fit on each PAM50 subtype
# lapply(split(CoxData_tcga_PAM, CoxData_tcga_PAM$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_THR25_quartiles, data = x)))

lapply(split(CoxData_tcga_PAM, CoxData_tcga_PAM$Pam50...Claudin.low.subtype),
       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_1_quartiles, data = x)))

#lapply(split(CoxData_tcga_PAM, CoxData_tcga_PAM$Pam50...Claudin.low.subtype),
#       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_2_quartiles, data = x)))

################################
# os: quartiles: Q1 vs Q4

# fit on each PAM50 subtype

# lapply(split(CoxData_tcga_PAM_Q1vsQ4_THR25, CoxData_tcga_PAM_Q1vsQ4_THR25$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_THR25_quartiles, data = x)))
# 
# lapply(split(CoxData_tcga_PAM_Q1vsQ4_CO130_1, CoxData_tcga_PAM_Q1vsQ4_CO130_1$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_1_quartiles, data = x)))
# 
# lapply(split(CoxData_tcga_PAM_Q1vsQ4_CO130_2, CoxData_tcga_PAM_Q1vsQ4_CO130_2$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_2_quartiles, data = x)))

################################
# rfs: quartiles: all
# fit on each PAM50 subtype
#lapply(split(CoxData_tcga_PAM, CoxData_tcga_PAM$Pam50...Claudin.low.subtype),
#       function(x) summary(coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_THR25_quartiles, data = x)))

lapply(split(CoxData_tcga_PAM, CoxData_tcga_PAM$Pam50...Claudin.low.subtype),
       function(x) summary(coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_1_quartiles, data = x)))

#lapply(split(CoxData_tcga_PAM, CoxData_tcga_PAM$Pam50...Claudin.low.subtype),
#       function(x) summary(coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_2_quartiles, data = x)))

################################
# rfs: quartiles: Q1 vs Q4

# fit on each PAM50 subtype

# lapply(split(CoxData_tcga_PAM_Q1vsQ4_THR25, CoxData_tcga_PAM_Q1vsQ4_THR25$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_THR25_quartiles, data = x)))
# 
# lapply(split(CoxData_tcga_PAM_Q1vsQ4_CO130_1, CoxData_tcga_PAM_Q1vsQ4_CO130_1$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_1_quartiles, data = x)))
# 
# lapply(split(CoxData_tcga_PAM_Q1vsQ4_CO130_2, CoxData_tcga_PAM_Q1vsQ4_CO130_2$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_2_quartiles, data = x)))



##############################################
## fit

## X3 classifier

# os: quartiles: all
# fit on each X3 subtypes
#lapply(split(CoxData_tcga, CoxData_tcga$X3.Gene.classifier.subtype),
#       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_THR25_quartiles, data = x)))

lapply(split(CoxData_tcga, CoxData_tcga$X3.Gene.classifier.subtype),
       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_1_quartiles, data = x)))

#lapply(split(CoxData_tcga, CoxData_tcga$X3.Gene.classifier.subtype),
#       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_2_quartiles, data = x)))

################################
# os: quartiles: Q1 vs Q4

# fit on each X3 subtypes

# lapply(split(CoxData_tcga_PAM_Q1vsQ4_THR25, CoxData_tcga_PAM_Q1vsQ4_THR25$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_THR25_quartiles, data = x)))
# 
# lapply(split(CoxData_tcga_PAM_Q1vsQ4_CO130_1, CoxData_tcga_PAM_Q1vsQ4_CO130_1$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_1_quartiles, data = x)))
# 
# lapply(split(CoxData_tcga_PAM_Q1vsQ4_CO130_2, CoxData_tcga_PAM_Q1vsQ4_CO130_2$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_2_quartiles, data = x)))

################################
# rfs: quartiles: all
# fit on each X3 subtypes
#lapply(split(CoxData_tcga, CoxData_tcga$X3.Gene.classifier.subtype),
#       function(x) summary(coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_THR25_quartiles, data = x)))

lapply(split(CoxData_tcga, CoxData_tcga$X3.Gene.classifier.subtype),
       function(x) summary(coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_1_quartiles, data = x)))

#lapply(split(CoxData_tcga, CoxData_tcga$X3.Gene.classifier.subtype),
#       function(x) summary(coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_2_quartiles, data = x)))

################################
# rfs: quartiles: Q1 vs Q4

# fit on each X3 subtypes

# lapply(split(CoxData_tcga_PAM_Q1vsQ4_THR25, CoxData_tcga_PAM_Q1vsQ4_THR25$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_THR25_quartiles, data = x)))
# 
# lapply(split(CoxData_tcga_PAM_Q1vsQ4_CO130_1, CoxData_tcga_PAM_Q1vsQ4_CO130_1$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_1_quartiles, data = x)))
# 
# lapply(split(CoxData_tcga_PAM_Q1vsQ4_CO130_2, CoxData_tcga_PAM_Q1vsQ4_CO130_2$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Disease.Free..Months., Disease.Free.Status) ~ tcga_prob_CO130_2_quartiles, data = x)))





















########################################################################  
########################################################################  
## Fit survival curves: TCGA: 

# fit surv curves: os

# surv_func_TCGA_os_THR25 <- function(x){
#   f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
#   return(surv_fit(f, data = CoxData_tcga_THR25))
# }
# 
# fit_list_tcga_os_THR25 <- lapply(pairs_list_THR25, surv_func_TCGA_os_THR25)
# 
# # calculate the pvalue
# Pval_list_tcga_os_THR25 <- surv_pvalue(fit_list_tcga_os_THR25)
# 
# Pval_df_tcga_os_THR25 <- do.call(rbind.data.frame, Pval_list_tcga_os_THR25)

#Pval_df_fil <- Pval_df[Pval_df$pval < 0.05, ] 

############
## Plot survival curves

# plot_list_tcga_os_THR25 <- ggsurvplot_list(fit_list_tcga_os_THR25, CoxData_tcga_THR25, legend.title = names(fit_list_tcga_os_THR25), pval = TRUE)
# 
# 
# Splot_tcga_os_THR25 <- arrange_ggsurvplots(plot_list_tcga_os_THR25, title = "Overall survival in the TCGA using the 10 pairs individually (THR25)", ncol = 2, nrow = 5)
# ggsave("./figures/sep28/THR25_10TSPs_tcga_os.pdf", Splot_tcga_os_THR25, width = 40, height = 40, units = "cm")


## TCGA all pairs
Fit_sig_TCGA_os_THR25 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_predClasses_THR25, data = CoxData_tcga)
Fit_sig_TCGA_os_CO130_1 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_predClasses_CO130_1, data = CoxData_tcga)
Fit_sig_TCGA_os_CO130_2 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_predClasses_CO130_2, data = CoxData_tcga)

## by quartiles
Fit_sig_TCGA_os_THR25_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_THR25_quartiles, data = CoxData_tcga)
Fit_sig_TCGA_os_CO130_1_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_1_quartiles, data = CoxData_tcga)
Fit_sig_TCGA_os_CO130_2_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_2_quartiles, data = CoxData_tcga)

pdf("./figures/logreg/THR25_Allpairs_TCGA_os.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_THR25,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 25 (logistic regression) and TCGA OS')
dev.off()


pdf("./figures/logreg/CO130_1_Allpairs_TCGA_os.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_CO130_1,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and TCGA OS')
dev.off()


pdf("./figures/logreg/CO130_2_Allpairs_TCGA_os.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_CO130_1,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_2 (logistic regression) and TCGA OS')
dev.off()

########
# by quartiles
pdf("./figures/logreg/THR25_tcga_os_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_THR25_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 25 (logistic regression) and TCGA OS: quartiles')
dev.off()

pdf("./figures/logreg/CO130_1_tcga_os_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_CO130_1_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and TCGA OS: quartiles')
dev.off()

pdf("./figures/logreg/CO130_2_tcga_os_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_CO130_2_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_2 (logistic regression) and TCGA OS: quartiles')
dev.off()


#############
## fit coxph model:

# by pair
# surv_func_TCGA_os_coxph_THR25 <- function(x){
#   f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
#   return(coxph(f, data = CoxData_tcga_THR25))
# }
# 
# fit_list_TCGA_os_coxph_THR25 <- lapply(pairs_list_THR25, surv_func_TCGA_os_coxph_THR25)
# names(fit_list_TCGA_os_coxph_THR25) <- pairs_list_THR25
# 
# summary_list_TCGA_os_coxph_THR25 <- lapply(fit_list_TCGA_os_coxph_THR25, summary)
# 
# # get the HR
# HR_list_TCGA_os_coxph_THR25 <- lapply(summary_list_TCGA_os_coxph_THR25, function(x){
#   HR <- x$conf.int[, 'exp(coef)']
#   Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
#   Pvalue_logrank_test <- x$sctest['pvalue']
#   Pvalue_wald_test <- x$waldtest['pvalue']
#   data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
#              Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
# })
# 
# 
# HR_df_TCGA_os_coxph_THR25 <- as.data.frame(do.call(rbind, HR_list_TCGA_os_coxph_THR25))
# HR_df_TCGA_os_coxph_THR25$variable <- rownames(HR_df_TCGA_os_coxph_THR25)
# HR_df_TCGA_os_coxph_THR25 <- HR_df_TCGA_os_coxph_THR25[order(HR_df_TCGA_os_coxph_THR25$HR, decreasing = T), ]
# 
# # save the results
# write.csv(HR_df_TCGA_os_coxph_THR25, 'objs/sep28/THR25_HR_df_TCGA_os_coxph.csv')


########
# by predictions
Fit_sig_TCGA_coxph_THR25 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_THR25, data = CoxData_tcga)
summary(Fit_sig_TCGA_coxph_THR25)

Fit_sig_TCGA_coxph_CO130_1 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_1, data = CoxData_tcga)
summary(Fit_sig_TCGA_coxph_CO130_1)

Fit_sig_TCGA_coxph_CO130_2 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO130_2, data = CoxData_tcga)
summary(Fit_sig_TCGA_coxph_CO130_2)


png('./figures/logreg/THR25_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_coxph_THR25, fontsize = 0.5)
dev.off()

png('./figures/logreg/CO130_1_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_coxph_CO130_1, fontsize = 0.5)
dev.off()

png('./figures/logreg/CO130_2_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_coxph_CO130_2, fontsize = 0.5)
dev.off()
########################################################################  
########################################################################  
## Fit survival curves: TCGA: 

# fit surv curves: PFS

# surv_func_TCGA_pfs_THR25 <- function(x){
#   f <- as.formula(paste("Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~", x))
#   return(surv_fit(f, data = CoxData_tcga_THR25))
# }
# 
# fit_list_tcga_pfs_THR25 <- lapply(pairs_list_THR25, surv_func_TCGA_pfs_THR25)
# 
# # calculate the pvalue
# Pval_list_tcga_pfs_THR25 <- surv_pvalue(fit_list_tcga_pfs_THR25)
# 
# Pval_df_tcga_pfs_THR25 <- do.call(rbind.data.frame, Pval_list_tcga_pfs_THR25)

#Pval_df_fil <- Pval_df[Pval_df$pval < 0.05, ] 

############
## Plot survival curves

# plot_list_tcga_pfs_THR25 <- ggsurvplot_list(fit_list_tcga_pfs_THR25, CoxData_tcga_THR25, legend.title = names(fit_list_tcga_pfs_THR25), pval = TRUE)
# 
# 
# Splot_tcga_pfs_THR25 <- arrange_ggsurvplots(plot_list_tcga_pfs_THR25, title = "Progression free survival in the TCGA using the 10 pairs individually (THR25)", ncol = 2, nrow = 5)
# ggsave("./figures/sep28/10TSPs_tcga_pfs_THR25.pdf", Splot_tcga_pfs_THR25, width = 40, height = 40, units = "cm")


## TCGA all genes
Fit_sig_TCGA_pfs_THR25 <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_predClasses_THR25, data = CoxData_tcga)
Fit_sig_TCGA_pfs_CO130_1 <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_predClasses_CO130_1, data = CoxData_tcga)
Fit_sig_TCGA_pfs_CO130_2 <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_predClasses_CO130_2, data = CoxData_tcga)

## by quartiles
Fit_sig_TCGA_pfs_THR25_quartiles <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_THR25_quartiles, data = CoxData_tcga)
Fit_sig_TCGA_pfs_CO130_1_quartiles <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_CO130_1_quartiles, data = CoxData_tcga)
Fit_sig_TCGA_pfs_CO130_2_quartiles <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_CO130_2_quartiles, data = CoxData_tcga)

pdf("./figures/logreg/THR25_Allpairs_TCGA_PFS.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_THR25,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 25 (logistic regression) and TCGA PFS')
dev.off()

pdf("./figures/logreg/CO130_1_Allpairs_TCGA_PFS.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_CO130_1,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and TCGA PFS')
dev.off()

pdf("./figures/logreg/CO130_2_Allpairs_TCGA_PFS.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_CO130_2,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_2 (logistic regression) and TCGA PFS')
dev.off()

########
# by quartiles
pdf("./figures/logreg/THR25_tcga_pfs_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_THR25_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 25 (logistic regression) and TCGA PFS: quartiles')
dev.off()

pdf("./figures/logreg/CO130_1_tcga_pfs_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_CO130_1_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and TCGA PFS: quartiles')
dev.off()

pdf("./figures/logreg/CO130_2_tcga_pfs_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_CO130_2_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_2 (logistic regression) and TCGA PFS: quartiles')
dev.off()
#############
## fit coxph model:

# by gene
# surv_func_TCGA_pfs_coxph_THR25 <- function(x){
#   f <- as.formula(paste("Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~", x))
#   return(coxph(f, data = CoxData_tcga_THR25))
# }
# 
# fit_list_TCGA_pfs_coxph_THR25 <- lapply(pairs_list_THR25, surv_func_TCGA_pfs_coxph_THR25)
# names(fit_list_TCGA_pfs_coxph_THR25) <- pairs_list_THR25
# 
# summary_list_TCGA_pfs_coxph_THR25 <- lapply(fit_list_TCGA_pfs_coxph_THR25, summary)
# 
# # get the HR
# HR_list_TCGA_pfs_coxph_THR25 <- lapply(summary_list_TCGA_pfs_coxph_THR25, function(x){
#   HR <- x$conf.int[, 'exp(coef)']
#   Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
#   Pvalue_logrank_test <- x$sctest['pvalue']
#   Pvalue_wald_test <- x$waldtest['pvalue']
#   data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
#              Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
# })
# 
# 
# HR_df_TCGA_pfs_coxph_THR25 <- as.data.frame(do.call(rbind, HR_list_TCGA_pfs_coxph_THR25))
# HR_df_TCGA_pfs_coxph_THR25$variable <- rownames(HR_df_TCGA_pfs_coxph_THR25)
# HR_df_TCGA_pfs_coxph_THR25 <- HR_df_TCGA_pfs_coxph_THR25[order(HR_df_TCGA_pfs_coxph_THR25$HR, decreasing = T), ]
# 
# # save the results
# write.csv(HR_df_TCGA_pfs_coxph_THR25, 'objs/sep28/HR_df_TCGA_pfs_coxph_THR25.csv')


###########
# by predictions
Fit_sig_TCGA_pfs_coxph_THR25 <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_THR25, data = CoxData_tcga)
summary(Fit_sig_TCGA_pfs_coxph_THR25)

Fit_sig_TCGA_pfs_coxph_CO130_1 <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_CO130_1, data = CoxData_tcga)
summary(Fit_sig_TCGA_pfs_coxph_CO130_1)

Fit_sig_TCGA_pfs_coxph_CO130_2 <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_CO130_2, data = CoxData_tcga)
summary(Fit_sig_TCGA_pfs_coxph_CO130_2)

png('./figures/logreg/THR25_HR_tcga_pfs.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_pfs_coxph_THR25, fontsize = 0.5)
dev.off()

png('./figures/logreg/CO130_1_HR_tcga_pfs.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_pfs_coxph_CO130_1, fontsize = 0.5)
dev.off()

png('./figures/logreg/CO130_2_HR_tcga_pfs.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_pfs_coxph_CO130_2, fontsize = 0.5)
dev.off()



###########################
## heatmap

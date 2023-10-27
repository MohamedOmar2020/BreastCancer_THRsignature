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

#############################
## load the ET signatures
All <- read_xlsx('./data/ET-9 Selection Steps.xlsx')
ET60 <- All$`ET-60`[!is.na(All$`ET-60`)]

#################
# THR signature
THR_signature <- readxl::read_xlsx("./data/THR Signatures_sep23.xlsx")

# get THR 50.1
THR_50 <- THR_signature$`THR-50.1`[!is.na(THR_signature$`THR-50.1`)]

THR_50 <- gsub('-', '', THR_50)


# combine
CO110 <- c(ET60, THR_50)
CO110 <- CO110[!CO110 %in% 'Gene']

#############################################
# Load the  expression and pheno data
load('./data/prostate/bulk/MetastasisDataGood.rda')
load("./data/prostate//for_survival.rda")

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


##############################################
# fix gene names
##############################################
setdiff(CO110, rownames(usedTrainMat))

# filter the CO110 signatures to include only the genes present in the expr matrices
CO110_fil <- CO110[CO110 %in% rownames(usedTrainMat)]
CO110_fil <- CO110[CO110 %in% rownames(usedTrainMat) & CO110 %in% rownames(Expr_tcga)]

##############################################
### combine in 1 dataset
##############################################
Data_train <- as.data.frame(cbind(t(usedTrainMat), usedTrainGroup))
Data_train$usedTrainGroup <- as.factor(Data_train$usedTrainGroup)
levels(Data_train$usedTrainGroup) <- c(0, 1)
colnames(Data_train)[colnames(Data_train) %in% c('usedTrainGroup')] <- c('label')

Data_test <- as.data.frame(cbind(t(usedTestMat), usedTestGroup))
Data_test$usedTestGroup <- as.factor(Data_test$usedTestGroup)
levels(Data_test$usedTestGroup) <- c(0, 1)
colnames(Data_test)[colnames(Data_test) %in% c('usedTestGroup')] <- c('label')

Data_tcga <- as.data.frame(cbind(t(Expr_tcga), group_tcga))
Data_tcga$group_tcga <- as.factor(Data_tcga$group_tcga)
levels(Data_tcga$group_tcga) <- c(0, 1)
colnames(Data_tcga)[colnames(Data_tcga) %in% c('group_tcga')] <- c('label')

#############################################################################################################
##############################################################################################################
# the model: trained on mets
CO110_model_prostate <- glm(as.formula((paste("label ~", paste(CO110_fil, collapse = "+")))), data = Data_train, family = "binomial")
summary(CO110_model_prostate)

save(CO110_model_prostate, file = 'objs/CO110_model_prostate.rda')

# the model: trained on survival
CO110_model_prostate_TCGA_PFS <- glm(as.formula((paste("label ~", paste(CO110_fil, collapse = "+")))), data = Data_tcga, family = "binomial")
summary(CO110_model_prostate_TCGA_PFS)

save(CO110_model_prostate_TCGA_PFS, file = 'objs/CO110_model_prostate_TCGA_PFS.rda')

###########################################################################
############################################################################
### predict in the training dataset
# Make predictions

Train_prob <- CO110_model_prostate %>% predict(Data_train , type = "response")

### Threshold
thr <- coords(roc(Data_train$label, Train_prob, levels = c(0, 1), direction = "<"), "best")["threshold"]
thr

### ROC Curve
ROCTrain <- roc(Data_train$label, Train_prob, plot = F, print.thres=thr$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c(0,1), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain

### Get predictions based on best threshold from ROC curve
Train_predClasses <- ifelse(Train_prob >= thr$threshold, '1', '0')
table(Train_predClasses)
Train_predClasses <- factor(Train_predClasses, levels = c('0', '1'))


### Resubstitution performance in the TRAINING set
ConfusionTrain <- confusionMatrix(Train_predClasses, Data_train$label, positive = '1', mode = "everything")
ConfusionTrain

## MCC
MCC_Train <- mltools::mcc(pred = Train_predClasses, actuals = Data_train$label)
MCC_Train


############################################################################
### predict in the testing dataset
# Make predictions

Test_prob <- CO110_model_prostate %>% predict(Data_test , type = "response")


### ROC Curve
ROCTest <- roc(Data_test$label, Test_prob, plot = F, print.thres=thr$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTest

### Get predictions based on best threshold from ROC curve
Test_predClasses <- ifelse(Test_prob >= thr$threshold, "1", "0")
table(Test_predClasses)
Test_predClasses <- factor(Test_predClasses, levels = c('0', '1'))


### Resubstitution performance in the TestING set
ConfusionTest <- confusionMatrix(Test_predClasses, Data_test$label, positive = "1", mode = "everything")
ConfusionTest

## MCC
MCC_Test <- mltools::mcc(pred = Test_predClasses, actuals = Data_test$label)
MCC_Test


#########################################################################
#########################################################################
### TCGA from CO110 trained on mets

tcga_prob <- CO110_model_prostate %>% predict(Data_tcga , type = "response")

### ROC
ROC_tcga <- roc(Data_tcga$label, tcga_prob, plot = F, print.thres=thr$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_tcga


### Get predictions based on best threshold from ROC curve
tcga_predClasses <- ifelse(tcga_prob >= thr$threshold, "1", "0")
table(tcga_predClasses)
tcga_predClasses <- factor(tcga_predClasses, levels = c('0', '1'))


### CI in tcga
Confusion_tcga <- confusionMatrix(tcga_predClasses, Data_tcga$label, positive = "1", mode = "everything")
Confusion_tcga


## MCC
MCC_tcga <- mltools::mcc(pred = tcga_predClasses, actuals = Data_tcga$label)
MCC_tcga

##############
### TCGA from CO110 trained on PFS

tcga_prob_PFS <- CO110_model_prostate_TCGA_PFS %>% predict(Data_tcga , type = "response")

### ROC
ROC_tcga_PFS <- roc(Data_tcga$label, tcga_prob_PFS, plot = F, print.thres=thr$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_tcga_PFS


### Get predictions based on best threshold from ROC curve
tcga_predClasses_PFS <- ifelse(tcga_prob_PFS >= thr$threshold, "1", "0")
table(tcga_predClasses_PFS)
tcga_predClasses_PFS <- factor(tcga_predClasses_PFS, levels = c('0', '1'))


### CI in tcga
Confusion_tcga_PFS <- confusionMatrix(tcga_predClasses_PFS, Data_tcga$label, positive = "1", mode = "everything")
Confusion_tcga_PFS


## MCC
MCC_tcga_PFS <- mltools::mcc(pred = tcga_predClasses_PFS, actuals = Data_tcga$label)
MCC_tcga_PFS

##########################
# load another TCGA metadata file to get the tumor purity score
Pheno_tcga2 <- read.delim2('./data/prostate/bulk/TCGA/tcga_meta_data_all.csv', sep = ',')
Pheno_tcga2 <- Pheno_tcga2[!duplicated(Pheno_tcga2$patient_id), ]
rownames(Pheno_tcga2) <- Pheno_tcga2$patient_id
CommonSamples <- intersect(rownames(Pheno_tcga), rownames(Pheno_tcga2))
Pheno_tcga2 <- Pheno_tcga2[CommonSamples, ]
all(rownames(Pheno_tcga2) == rownames(Pheno_tcga))

# add the gleason score to CoxData_tcga
Pheno_tcga$purity <- Pheno_tcga2$purity
Pheno_tcga$purity <- as.numeric(Pheno_tcga$purity)

## Keep only the relevant information (Metastasis Event and Time)
Phenotype_tcga <- cbind(Pheno_tcga[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Progression.Free.Status", "Progress.Free.Survival..Months.", "purity")], 
                        tcga_prob, tcga_predClasses, 
                        tcga_prob_PFS, tcga_predClasses_PFS)

#Expr_metabric <- Expr_metabric[ClassifierGenes, ]
#Expr_tcga <- Expr_tcga[ClassifierGenes, ]

# create a merged pdata and Z-scores object
CoxData_tcga <- data.frame(Phenotype_tcga)

# divide the CO110 probability into quartiles
CoxData_tcga <- CoxData_tcga %>%
  mutate(CO110_quartiles = ntile(tcga_prob, 4)) %>%
  mutate(CO110_quartiles_PFS = ntile(tcga_prob_PFS, 4))


########################################################################  
## Fit survival curves

# OS
## TCGA all genes
Fit_sig_tcga_os <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_predClasses, data = CoxData_tcga)
Fit_sig_tcga_CO110PFS_os <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_predClasses_PFS, data = CoxData_tcga)


## by PRN quartiles
Fit_sig_tcga_os_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ CO110_quartiles, data = CoxData_tcga)
Fit_sig_tcga_CO110PFS_os_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ CO110_quartiles_PFS, data = CoxData_tcga)

################
# PFS
## metabric all genes
Fit_sig_tcga_pfs <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_predClasses, data = CoxData_tcga)
Fit_sig_tcga_CO110PFS_pfs <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_predClasses_PFS, data = CoxData_tcga)

## by PRN quartiles
Fit_sig_tcga_pfs_quartiles <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ CO110_quartiles, data = CoxData_tcga)
Fit_sig_tcga_CO110PFS_pfs_quartiles <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ CO110_quartiles_PFS, data = CoxData_tcga)


############################################################################
############################################################################
# plot OS

pdf("./figures/CO110_prostate/CO110_tcga_os_all.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_os,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 10,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'CO110 signature and TCGA OS')
dev.off()

pdf("./figures/CO110_prostate/CO110_tcga_CO110PFS_os_all.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_CO110PFS_os,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 10,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'CO110 signature and TCGA OS')
dev.off()


#############
# by quartiles
pdf("./figures/CO110_prostate/CO110_tcga_os_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_os_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 10,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'CO110 signature and TCGA OS: quartiles')
dev.off()

pdf("./figures/CO110_prostate/CO110_tcga_CO110PFS_os_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga_CO110PFS_os_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 10,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'CO110 signature and TCGA OS: quartiles')
dev.off()


######################################
# plot PFS

tiff("./figures/CO110_prostate/CO110_tcga_pfs_all.tiff", width = 2000, height = 2000, res = 400)
ggsurvplot(Fit_sig_tcga_pfs,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 8,
           legend.labs = c('predicted PFS: 0', 'predicted PFS: 1'),
           legend.title	= '',
           ggtheme = theme_survminer(base_size = 12, font.x = c(12, 'bold.italic', 'black'), font.y = c(12, 'bold.italic', 'black'), font.tickslab = c(12, 'plain', 'black'), font.legend = c(12, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'PRN signature and TCGA PFS'
)
dev.off()

tiff("./figures/CO110_prostate/CO110_tcga_CO110PFS_pfs_all.tiff", width = 2000, height = 2000, res = 400)
ggsurvplot(Fit_sig_tcga_CO110PFS_pfs,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 8,
           legend.labs = c('predicted PFS: 0', 'predicted PFS: 1'),
           legend.title	= '',
           ggtheme = theme_survminer(base_size = 12, font.x = c(12, 'bold.italic', 'black'), font.y = c(12, 'bold.italic', 'black'), font.tickslab = c(12, 'plain', 'black'), font.legend = c(12, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'PRN signature and TCGA PFS'
)
dev.off()


########
# by quartiles
tiff("./figures/CO110_prostate/tcga_PFS_quartiles.tiff", width = 2000, height = 2000, res = 400)
ggsurvplot(Fit_sig_tcga_pfs_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 10,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           legend.title	= '',
           ggtheme = theme_survminer(base_size = 12, font.x = c(12), font.y = c(12, 'bold.italic', 'black'), font.tickslab = c(12, 'plain', 'black'), font.legend = c(12, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'PRN signature and TCGA PFS: quartiles'
)
dev.off()

tiff("./figures/CO110_prostate/tcga_CO110PFS_PFS_quartiles.tiff", width = 2000, height = 2000, res = 400)
ggsurvplot(Fit_sig_tcga_CO110PFS_pfs_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 10,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           legend.title	= '',
           ggtheme = theme_survminer(base_size = 12, font.x = c(12), font.y = c(12, 'bold.italic', 'black'), font.tickslab = c(12, 'plain', 'black'), font.legend = c(12, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'PRN signature and TCGA PFS: quartiles'
)
dev.off()

##############################################################################
## fit coxph model:

########
## PFS

# by probability

Fit_sig_tcga_PFS_coxph_logReg <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_logReg, data = CoxData_tcga)
summary(Fit_sig_tcga_PFS_coxph_logReg)

########
## by quartiles

# make a factor with Q1 (lowest risk) being the reference
CoxData_tcga$quartiles <- factor(CoxData_tcga$quartiles, levels = c('1', '2', '3', '4'))
levels(CoxData_tcga$quartiles) <- paste0('Q', levels(CoxData_tcga$quartiles))

# fit
Fit_sig_tcga_pfs_coxph_logReg_quartiles <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ quartiles, data = CoxData_tcga)
summary_tcga_pfs_coxph_logReg_quartiles <- summary(Fit_sig_tcga_pfs_coxph_logReg_quartiles)

##########################
## multivariate cox with gleason
# read another version of the clinical data (Firehose legacy) which contains gleason score
pheno2 <- read.delim2('data/bulk/TCGA/prad_tcga_clinical_data.tsv')
pheno2$Patient.ID <- gsub("\\-", "\\.", pheno2$Patient.ID)
pheno2 <- pheno2[!duplicated(pheno2$Patient.ID), ]
rownames(pheno2) <- pheno2$Patient.ID
CommonSamples <- intersect(rownames(Pheno_tcga), rownames(pheno2))
pheno2 <- pheno2[CommonSamples, ]
all(rownames(pheno2) == rownames(CoxData_tcga))

# add the gleason score to CoxData_tcga
CoxData_tcga$gleason <- pheno2$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer
CoxData_tcga$gleason <- as.factor(CoxData_tcga$gleason)
levels(CoxData_tcga$gleason)

colnames(CoxData_tcga)[colnames(CoxData_tcga) == 'tcga_prob_logReg'] <- 'PRN signature'

# fit the multivariate COX with gleason as cofactor 
Fit_sig_tcga_PFS_coxph_logReg_withGS <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ `PRN signature` + gleason, data = CoxData_tcga)
summary(Fit_sig_tcga_PFS_coxph_logReg_withGS)

tiff('./figures/survival/multivariateCox.tiff', width = 2500, height = 3000, res = 400)
ggforest(Fit_sig_tcga_PFS_coxph_logReg_withGS, fontsize = 1, main = 'HR')
dev.off()

###################
# fit multivariate COX with gleason and tumor purity as cofactor 
Fit_sig_tcga_PFS_coxph_logReg_withGS_purity <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ `PRN signature` + gleason + purity, data = CoxData_tcga)
summary(Fit_sig_tcga_PFS_coxph_logReg_withGS_purity)

tiff('./figures/survival/multivariateCox_gleason_purity.tiff', width = 2500, height = 3000, res = 400)
ggforest(Fit_sig_tcga_PFS_coxph_logReg_withGS_purity, fontsize = 1, main = 'HR')
dev.off()

###################
# fit multivariate COX with tumor purity as cofactor 
Fit_sig_tcga_PFS_coxph_logReg_withPurity <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ `PRN signature` + purity, data = CoxData_tcga)
summary(Fit_sig_tcga_PFS_coxph_logReg_withPurity)

tiff('./figures/survival/multivariateCox_purity.tiff', width = 2500, height = 3000, res = 400)
ggforest(Fit_sig_tcga_PFS_coxph_logReg_withPurity, fontsize = 1, main = 'HR')
dev.off()



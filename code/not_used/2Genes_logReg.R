###########################################################################


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

#################
################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')

### combine in 1 dataset: Training
Data_metabric <- as.data.frame(cbind(t(Expr_metabric), group_metabric))
Data_metabric$group_metabric <- as.factor(Data_metabric$group_metabric)
levels(Data_metabric$group_metabric) <- c('0', '1')
colnames(Data_metabric)[colnames(Data_metabric) %in% c('group_metabric')] <- c('os')

# alternative
#Data_metabric <- as.data.frame(cbind(t(Expr_metabric), RFS = Pheno_metabric$Relapse.Free.Status))
#Data_metabric$RFS <- as.factor(Data_metabric$RFS)
#levels(Data_metabric$RFS) <- c('0', '1')

##########
### tcga
Data_tcga <- as.data.frame(cbind(t(Expr_tcga), group_tcga))
Data_tcga$group_tcga <- as.factor(Data_tcga$group_tcga)
levels(Data_tcga$group_tcga) <- c('0', '1')
colnames(Data_tcga)[colnames(Data_tcga) %in% c('group_tcga')] <- c('os')

# alternative
#Data_tcga <- as.data.frame(cbind(t(Expr_tcga), RFS = Pheno_tcga$Disease.Free.Status))

###########################################################################
### TRAINING using logistic regression
###########################################################################


# the model

twoGnsModel_model <- glm(os ~ KLF7 + ZNF627, data = Data_metabric, family = "binomial")
summary(twoGnsModel_model)

save(twoGnsModel_model, file = "./objs/twoGnsModel_model_logreg.rda")

###########################################################################
############################################################################
### predict in the training dataset
# Make predictions

Train_prob_twoGnsModel <- twoGnsModel_model %>% predict(Data_metabric , type = "response")

### Threshold
thr_twoGnsModel <- coords(roc(group_metabric, Train_prob_twoGnsModel, levels = c("0", "1"), direction = "<"), "best")["threshold"]
thr_twoGnsModel

### ROC Curve
ROCTrain_twoGnsModel <- roc(group_metabric, Train_prob_twoGnsModel, plot = F, print.thres=thr_twoGnsModel$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("0", "1"), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_twoGnsModel

### Get predictions based on best threshold from ROC curve
Train_predClasses_twoGnsModel <- ifelse(Train_prob_twoGnsModel >= thr_twoGnsModel$threshold, "1", "0")
table(Train_predClasses_twoGnsModel)
Train_predClasses_twoGnsModel <- factor(Train_predClasses_twoGnsModel, levels = c('0', '1'))


### Resubstitution performance in the TRAINING set
ConfusionTrain_twoGnsModel <- confusionMatrix(Train_predClasses_twoGnsModel, group_metabric, positive = "1", mode = "everything")
ConfusionTrain_twoGnsModel

## MCC
MCC_Train_twoGnsModel <- mltools::mcc(pred = Train_predClasses_twoGnsModel, actuals = group_metabric)
MCC_Train_twoGnsModel

#########################################################################
#########################################################################
### Testing
tcga_prob_twoGnsModel <-  twoGnsModel_model %>% predict(Data_tcga , type = "response")


### ROC
ROC_tcga_twoGnsModel <- roc(group_tcga, tcga_prob_twoGnsModel, plot = F, print.thres=thr_twoGnsModel$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("0", "1"), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_tcga_twoGnsModel

############################
### Get predictions based on best threshold from ROC curve
tcga_predClasses_twoGnsModel <- ifelse(tcga_prob_twoGnsModel >= thr_twoGnsModel$threshold, "1", "0")
table(tcga_predClasses_twoGnsModel)
tcga_predClasses_twoGnsModel <- factor(tcga_predClasses_twoGnsModel, levels = c('0', '1'))

##################################
### CI  in testing 1
Confusion_tcga_twoGnsModel <- confusionMatrix(tcga_predClasses_twoGnsModel, group_tcga, positive = "1", mode = "everything")
Confusion_tcga_twoGnsModel

################
## MCC
MCC_tcga_twoGnsModel <- mltools::mcc(pred = tcga_predClasses_twoGnsModel, actuals = group_tcga)
MCC_tcga_twoGnsModel


##########################
## Keep only the relevant information (Metastasis Event and Time)
Phenotype_metabric <- cbind(Pheno_metabric[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Relapse.Free.Status", "Relapse.Free.Status..Months.", "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC", "X3.Gene.classifier.subtype")], 
                            Train_prob_twoGnsModel, Train_predClasses_twoGnsModel)

Phenotype_tcga <- cbind(Pheno_tcga[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Progression.Free.Status", "Progress.Free.Survival..Months.")], 
                        tcga_prob_twoGnsModel, tcga_predClasses_twoGnsModel)

#Expr_metabric <- Expr_metabric[ClassifierGenes, ]
#Expr_tcga <- Expr_tcga[ClassifierGenes, ]


# create a merged pdata and Z-scores object
CoxData_metabric <- data.frame(Phenotype_metabric)
CoxData_tcga <- data.frame(Phenotype_tcga)

# divide the probabilities into quartiles
CoxData_metabric <- CoxData_metabric %>%
  mutate(metabric_prob_twoGnsModel_quartiles = ntile(Train_prob_twoGnsModel, 4), 
         metabric_prob_twoGnsModel_quintiles = ntile(Train_prob_twoGnsModel, 5))

CoxData_tcga <- CoxData_tcga %>%
  mutate(tcga_prob_twoGnsModel_quartiles = ntile(tcga_prob_twoGnsModel, 4),
         tcga_prob_twoGnsModel_quintiles = ntile(tcga_prob_twoGnsModel, 5))

#CutPoint <- surv_cutpoint(data = CoxData, time = "Time", event = "Event", variables = "ResidualDisease_Score")
#CutPoint

#SurvData <- surv_categorize(CutPoint)

#SurvData$ResidualDisease_Score <- factor(SurvData$ResidualDisease_Score, levels = c("low", "high"))


########################################################################  
## Fit survival curves

# Metabric: 

# init a list for classifier pairs
# THR_25_fil_list <- as.list(unlist(THR_25_fil))
# THR_50_1_fil_list <- as.list(unlist(THR_50_1_fil))
# THR_50_2_fil_list <- as.list(unlist(THR_50_2_fil))
# 
# 
# names(THR_25_fil_list) <- THR_25_fil
# names(THR_50_1_fil_list) <- THR_50_1_fil
# names(THR_50_2_fil_list) <- THR_50_2_fil
# 
# surv_func_metabric_os <- function(x){
#   f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
#   return(surv_fit(f, data = CoxData_metabric))
# }
# 
# fit_list_metabric_os <- lapply(THR_25_fil_list, surv_func_metabric_os)
# 
# # calculate the pvalue
# Pval_list_metabric_os <- surv_pvalue(fit_list_metabric_os)
# 
# Pval_df_metabric_os <- do.call(rbind.data.frame, Pval_list_metabric_os)

############
## Plot survival curves

#plot_list_metabric_os <- ggsurvplot_list(fit_list_metabric_os, CoxData_metabric, legend.title = names(fit_list_metabric_os), pval = TRUE)


#Splot_metabric_os <- arrange_ggsurvplots(plot_list_metabric_os, title = "Survival plots using the 10 pairs individually (THR25)", ncol = 2, nrow = 5)
#ggsave("./figures/sep28/THR25_Metabric_os.pdf", Splot_metabric_os, width = 40, height = 40, units = "cm")


# OS
## metabric all genes
Fit_sig_metabric_os_twoGnsModel <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_twoGnsModel, data = CoxData_metabric)


## by quartiles
Fit_sig_metabric_os_twoGnsModel_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quartiles, data = CoxData_metabric)

## by quintiles
Fit_sig_metabric_os_twoGnsModel_quintiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quintiles, data = CoxData_metabric)


# RFS
## metabric all genes
Fit_sig_metabric_RFS_twoGnsModel <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_predClasses_twoGnsModel, data = CoxData_metabric)

## by quartiles
Fit_sig_metabric_RFS_twoGnsModel_quartiles <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quartiles, data = CoxData_metabric)

## by quintiles
Fit_sig_metabric_RFS_twoGnsModel_quintiles <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quintiles, data = CoxData_metabric)

#################################
## by clinical groups

## pam50

# keep only the major pam50 subtypes
CoxData_metabric_PAM <- CoxData_metabric %>%
  filter(Pam50...Claudin.low.subtype %in% c('Basal', 'claudin-low', 'Her2', 'LumA', 'LumB'))

# keep only quartiles 1 and 4 (has to be in each of THR25, twoGnsModel, and THR50_2)
CoxData_metabric_PAM_Q1vsQ4_twoGnsModel <- CoxData_metabric_PAM %>%
  filter(metabric_prob_twoGnsModel_quartiles %in% c('1', '4'))

CoxData_metabric_PAM_Q1vsQ5_twoGnsModel <- CoxData_metabric_PAM %>%
  filter(metabric_prob_twoGnsModel_quintiles %in% c('1', '5'))


# os: all
Fit_sig_metabric_os_twoGnsModel_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_twoGnsModel + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM)

# os: quartiles: all
Fit_sig_metabric_os_twoGnsModel_quartiles_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quartiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM)

# os: quintiles: all
Fit_sig_metabric_os_twoGnsModel_quintiles_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quintiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM)

# os: quartiles: Q1 vs Q4
Fit_sig_metabric_os_twoGnsModel_Q1vsQ4_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quartiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM_Q1vsQ4_twoGnsModel)

# os: quintiles: Q1 vs Q5
Fit_sig_metabric_os_twoGnsModel_Q1vsQ5_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quintiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM_Q1vsQ5_twoGnsModel)

##########
# rfs: all
Fit_sig_metabric_rfs_twoGnsModel_PAM <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_predClasses_twoGnsModel + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM)

# os: quartiles: all
Fit_sig_metabric_rfs_twoGnsModel_quartiles_PAM <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quartiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM)

# rfs: quintiles: all
Fit_sig_metabric_rfs_twoGnsModel_quintiles_PAM <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quintiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM)

# rfs: quartiles: Q1 vs Q4
Fit_sig_metabric_rfs_twoGnsModel_Q1vsQ4_PAM <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quartiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM_Q1vsQ4_twoGnsModel)

# rfs: quintiles: Q1 vs Q5
Fit_sig_metabric_rfs_twoGnsModel_Q1vsQ5_PAM <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quintiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM_Q1vsQ5_twoGnsModel)

#############
## ER

# keep only quartiles 1 and 4 (has to be in each of THR25, twoGnsModel, and THR50_2)

# Note: we can use this for both ER and X3
CoxData_metabric_Q1vsQ4_twoGnsModel <- CoxData_metabric %>%
  filter(metabric_prob_twoGnsModel_quartiles %in% c('1', '4'))

CoxData_metabric_Q1vsQ5_twoGnsModel <- CoxData_metabric %>%
  filter(metabric_prob_twoGnsModel_quintiles %in% c('1', '5'))

# os: all
Fit_sig_metabric_os_twoGnsModel_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_twoGnsModel + ER.status.measured.by.IHC, data = CoxData_metabric)

# os: quartiles: all
Fit_sig_metabric_os_twoGnsModel_quartiles_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quartiles + ER.status.measured.by.IHC, data = CoxData_metabric)

# os: quintiles: all
Fit_sig_metabric_os_twoGnsModel_quintiles_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quintiles + ER.status.measured.by.IHC, data = CoxData_metabric)

# os: quartiles: Q1 vs Q4
Fit_sig_metabric_os_twoGnsModel_Q1vsQ4_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quartiles + ER.status.measured.by.IHC, data = CoxData_metabric_Q1vsQ4_twoGnsModel)

# os: quintiles: Q1 vs Q5
Fit_sig_metabric_os_twoGnsModel_Q1vsQ5_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quintiles + ER.status.measured.by.IHC, data = CoxData_metabric_Q1vsQ5_twoGnsModel)

##############
# rfs: all
Fit_sig_metabric_rfs_twoGnsModel_ER <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_predClasses_twoGnsModel + ER.status.measured.by.IHC, data = CoxData_metabric)

# rfs: quartiles: all
Fit_sig_metabric_rfs_twoGnsModel_quartiles_ER <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quartiles + ER.status.measured.by.IHC, data = CoxData_metabric)

# rfs: quintiles: all
Fit_sig_metabric_rfs_twoGnsModel_quintiles_ER <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quintiles + ER.status.measured.by.IHC, data = CoxData_metabric)

# rfs: quartiles: Q1 vs Q4
Fit_sig_metabric_rfs_twoGnsModel_Q1vsQ4_ER <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quartiles + ER.status.measured.by.IHC, data = CoxData_metabric_Q1vsQ4_twoGnsModel)

# rfs: quintiles: Q1 vs Q5
Fit_sig_metabric_rfs_twoGnsModel_Q1vsQ5_ER <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quintiles + ER.status.measured.by.IHC, data = CoxData_metabric_Q1vsQ5_twoGnsModel)

#######################################
## X3

# os: all
Fit_sig_metabric_os_twoGnsModel_X3 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_twoGnsModel + X3.Gene.classifier.subtype, data = CoxData_metabric)

# os: quartiles: all
Fit_sig_metabric_os_twoGnsModel_quartiles_X3 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quartiles + X3.Gene.classifier.subtype, data = CoxData_metabric)

# os: quintiles: all
Fit_sig_metabric_os_twoGnsModel_quintiles_X3 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quintiles + X3.Gene.classifier.subtype, data = CoxData_metabric)

# os: quartiles: Q1 vs Q4
Fit_sig_metabric_os_twoGnsModel_Q1vsQ4_X3 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quartiles + X3.Gene.classifier.subtype, data = CoxData_metabric_Q1vsQ4_twoGnsModel)

# os: quintiles: Q1 vs Q5
Fit_sig_metabric_os_twoGnsModel_Q1vsQ5_X3 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quintiles + X3.Gene.classifier.subtype, data = CoxData_metabric_Q1vsQ5_twoGnsModel)

###########
# rfs: all
Fit_sig_metabric_rfs_twoGnsModel_X3 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_predClasses_twoGnsModel + X3.Gene.classifier.subtype, data = CoxData_metabric)

# rfs: quartiles: all
Fit_sig_metabric_rfs_twoGnsModel_quartiles_X3 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quartiles + X3.Gene.classifier.subtype, data = CoxData_metabric)

# rfs: quintiles: all
Fit_sig_metabric_rfs_twoGnsModel_quintiles_X3 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quintiles + X3.Gene.classifier.subtype, data = CoxData_metabric)

# rfs: quartiles: Q1 vs Q4
Fit_sig_metabric_rfs_twoGnsModel_Q1vsQ4_X3 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quartiles + X3.Gene.classifier.subtype, data = CoxData_metabric_Q1vsQ4_twoGnsModel)

# rfs: quintiles: Q1 vs Q5
Fit_sig_metabric_rfs_twoGnsModel_Q1vsQ5_X3 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quintiles + X3.Gene.classifier.subtype, data = CoxData_metabric_Q1vsQ5_twoGnsModel)

############################################################################
############################################################################
# plot OS

tiff("./figures/twoGns/twoGnsModel_metabric_os_allpairs.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_twoGnsModel,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC OS'
)
dev.off()

########
# by quartiles
tiff("./figures/twoGns/twoGnsModel_metabric_os_quartiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_twoGnsModel_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC OS: quartiles'
)
dev.off()

#############
# by quintiles
tiff("./figures/twoGns/twoGnsModel_metabric_os_quintiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_twoGnsModel_quintiles,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC OS: quintiles'
)
dev.off()

######################################
# plot RFS

tiff("./figures/twoGns/twoGnsModel_metabric_RFS.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_twoGnsModel,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'KLF7+ZNF627 and METABRIC RFS'
)
dev.off()

########
# by quartiles
tiff("./figures/twoGns/twoGnsModel_metabric_RFS_quartiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_twoGnsModel_quartiles,
           risk.table = FALSE,
           pval = FALSE,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'KLF7+ZNF627 and METABRIC RFS: quartiles'
)
dev.off()

########
# by quintiles

tiff("./figures/twoGns/twoGnsModel_metabric_RFS_quintiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_twoGnsModel_quintiles,
           risk.table = FALSE,
           pval = FALSE,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC RFS: quintiles'
)
dev.off()

############################################################################
############################################################################
### by clinical group

## PAM50

# OS

tiff("./figures/twoGns/twoGnsModel_metabric_os_PAM50.tiff", width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_os_twoGnsModel_PAM,
           risk.table = FALSE,
           pval = FALSE,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'KLF7+ZNF627 and METABRIC OS by PAM50 subtypes'
)
dev.off()


#########
# RFS


tiff("./figures/twoGns/twoGnsModel_metabric_rfs_PAM50.tiff", width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_rfs_twoGnsModel_PAM,
           risk.table = FALSE,
           pval = FALSE,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'KLF7+ZNF627 and METABRIC RFS by PAM50 subtypes'
)
dev.off()

######################################################
# OS: quartiles: all


tiff("./figures/twoGns/twoGnsModel_metabric_os_PAM50_quartiles.tiff", width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_os_twoGnsModel_quartiles_PAM,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           legend.title	= 'Quartiles',
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC OS by PAM50 subtypes: quartiles'
)
dev.off()

####################################################
## OS: quintiles: all
tiff("./figures/twoGns/twoGnsModel_metabric_os_PAM50_quintiles.tiff", width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_os_twoGnsModel_quintiles_PAM,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
           legend.title	= 'Quintiles',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC OS by PAM50 subtypes: quintiles'
)
dev.off()

######################################################
# OS: quartiles: Q1 vs Q4
tiff("./figures/twoGns/twoGnsModel_metabric_os_PAM50_Q1vsQ4.tiff", width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_os_twoGnsModel_Q1vsQ4_PAM,
           risk.table = FALSE,
           pval = FALSE,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           legend.title	= 'Quartiles',
           legend.labs = c('Q1', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627  and METABRIC OS by PAM50 subtypes: Q1 vs Q4'
)
dev.off()

######################################################
# OS: quintiles: Q1 vs Q5
tiff("./figures/twoGns/twoGnsModel_metabric_os_PAM50_Q1vsQ5.tiff", width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_os_twoGnsModel_Q1vsQ5_PAM,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           legend.labs = c('Q1', 'Q5'),
           legend.title	= 'Quintiles',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC OS by PAM50 subtypes: Q1 vs Q5'
)
dev.off()

#############################
# RFS: quartiles: all

tiff("./figures/twoGns/twoGnsModel_metabric_rfs_PAM50_quartiles.tiff", width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_rfs_twoGnsModel_quartiles_PAM,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           legend.title	= 'Quartiles',
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC RFS by PAM50 subtypes: quartiles'
)
dev.off()

####################################################
## RFS: quintiles: all
tiff("./figures/twoGns/twoGnsModel_metabric_rfs_PAM50_quintiles.tiff", width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_rfs_twoGnsModel_quintiles_PAM,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
           legend.title	= 'Quintiles',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC RFS by PAM50 subtypes: quintiles'
)
dev.off()

#############################
# RFS: quartiles: Q1 vs Q4
tiff("./figures/twoGns/twoGnsModel_metabric_rfs_PAM50_Q1vsQ4.tiff", width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_rfs_twoGnsModel_Q1vsQ4_PAM,
           risk.table = FALSE,
           pval = F,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           #legend.labs = c('Basal', 'Claudin-low', 'Her2+', 'Luminal A', 'Luminal B'),
           legend.title	= 'Quartiles',
           pval.size = 15,
           #break.x.by = 20,
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q4'),
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC RFS by PAM50 subtypes: Q1 vs Q4'
)
dev.off()

#############################
# RFS: quintiles: Q1 vs Q5
tiff("./figures/twoGns/twoGnsModel_metabric_rfs_PAM50_Q1vsQ5.tiff", width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_rfs_twoGnsModel_Q1vsQ5_PAM,
           risk.table = FALSE,
           pval = T,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           #legend.labs = c('Basal', 'Claudin-low', 'Her2+', 'Luminal A', 'Luminal B'),
           legend.title	= 'Quintiles',
           pval.size = 15,
           #break.x.by = 20,
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q5'),
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC RFS by PAM50 subtypes: Q1 vs Q5'
)
dev.off()
##############################################################################################
##############################################################################################
## ER

# OS


# pdf("./figures/logreg/logistic_regression_oct11/metabric/twoGnsModel/OS/2_groups/twoGnsModel_metabric_os_ER.pdf", width = 8, height = 8, onefile = F)
# ggsurvplot(Fit_sig_metabric_os_twoGnsModel_ER,
#            risk.table = FALSE,
#            pval = TRUE,
#            short.panel.labs = T,
#            facet.by = "ER.status.measured.by.IHC",
#            ggtheme = theme_minimal(),
#            risk.table.y.text.col = FALSE,
#            risk.table.y.text = FALSE, title = 'THR 50_1 and METABRIC OS by ER status')
# dev.off()


######################
# RFS
# pdf("./figures/logreg/logistic_regression_oct11/metabric/twoGnsModel/RFS/2_groups/twoGnsModel_metabric_rfs_ER.pdf", width = 8, height = 8, onefile = F)
# ggsurvplot(Fit_sig_metabric_rfs_twoGnsModel_ER,
#            risk.table = FALSE,
#            pval = TRUE,
#            short.panel.labs = T,
#            facet.by = "ER.status.measured.by.IHC",
#            ggtheme = theme_minimal(),
#            risk.table.y.text.col = FALSE,
#            risk.table.y.text = FALSE, title = 'THR 50_1 and METABRIC RFS by ER status')
# dev.off()

####################################################
# OS: quartiles: all

# pdf("./figures/logreg/logistic_regression_oct11/metabric/twoGnsModel/OS/all_quartiles/twoGnsModel_metabric_os_ER_quartiles.pdf", width = 8, height = 8, onefile = F)
# ggsurvplot(Fit_sig_metabric_os_twoGnsModel_quartiles_ER,
#            risk.table = FALSE,
#            pval = TRUE,
#            short.panel.labs = T,
#            facet.by = "ER.status.measured.by.IHC",
#            ggtheme = theme_minimal(),
#            legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
#            risk.table.y.text.col = FALSE,
#            risk.table.y.text = FALSE, title = 'THR 50_1 and METABRIC OS by ER status: quartiles')
# dev.off()

####################################################
# OS: quartiles: Q1 vs Q4


# pdf("./figures/logreg/logistic_regression_oct11/metabric/twoGnsModel/OS/Q1_vs_Q4/twoGnsModel_metabric_os_ER_Q1vsQ4.pdf", width = 8, height = 8, onefile = F)
# ggsurvplot(Fit_sig_metabric_os_twoGnsModel_Q1vsQ4_ER,
#            risk.table = FALSE,
#            pval = TRUE,
#            short.panel.labs = T,
#            facet.by = "ER.status.measured.by.IHC",
#            ggtheme = theme_minimal(),
#            legend.labs = c('Q1', 'Q4'),
#            risk.table.y.text.col = FALSE,
#            risk.table.y.text = FALSE, title = 'THR 50_1 and METABRIC OS by ER status: Q1 vs Q4')
# dev.off()


####################################################
# RFS: quartiles: all

# pdf("./figures/logreg/logistic_regression_oct11/metabric/twoGnsModel/RFS/all_quartiles/twoGnsModel_metabric_rfs_ER_quartiles.pdf", width = 8, height = 8, onefile = F)
# ggsurvplot(Fit_sig_metabric_rfs_twoGnsModel_quartiles_ER,
#            risk.table = FALSE,
#            pval = TRUE,
#            short.panel.labs = T,
#            facet.by = "ER.status.measured.by.IHC",
#            ggtheme = theme_minimal(),
#            risk.table.y.text.col = FALSE,
#            legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
#            risk.table.y.text = FALSE, title = 'THR 50_1 and METABRIC RFS by ER status: quartiles')
# dev.off()

####################################################
# RFS: quartiles: Q1 vs Q4

# pdf("./figures/logreg/logistic_regression_oct11/metabric/twoGnsModel/RFS/Q1_vs_Q4/twoGnsModel_metabric_rfs_ER_Q1vsQ4.pdf", width = 8, height = 8, onefile = F)
# ggsurvplot(Fit_sig_metabric_rfs_twoGnsModel_Q1vsQ4_ER,
#            risk.table = FALSE,
#            pval = TRUE,
#            short.panel.labs = T,
#            facet.by = "ER.status.measured.by.IHC",
#            ggtheme = theme_minimal(),
#            risk.table.y.text.col = FALSE,
#            legend.labs = c('Q1', 'Q4'),
#            risk.table.y.text = FALSE, title = 'THR 50_1 and METABRIC RFS by ER status: Q1 vs Q4')
# dev.off()

##############################################################################################
##############################################################################################
## X3

# OS
tiff("./figures/twoGns/twoGnsModel_metabric_os_X3.tiff", width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_os_twoGnsModel_X3,
           risk.table = FALSE,
           pval = FALSE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC OS by X3 classifier subtypes'
)
dev.off()

######################
# RFS
tiff("./figures/twoGns/twoGnsModel_metabric_rfs_X3.tiff", width = 3000, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_rfs_twoGnsModel_X3,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC RFS by X3 classifier subtypes'
)
dev.off()

####################################################
# OS: quartiles: all

tiff("./figures/twoGns/twoGnsModel_metabric_os_X3_quartiles.tiff", width = 3000, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_os_twoGnsModel_quartiles_X3,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           legend.title	= 'Quartiles',
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC OS by X3 classifier subtypes: quartiles'
)
dev.off()

####################################################
## OS: quintiles: all
tiff("./figures/twoGns/twoGnsModel_metabric_os_X3_quintiles.tiff", width = 3000, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_os_twoGnsModel_quintiles_X3,
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
           title = 'KLF7+ZNF627 and METABRIC OS by X3 subtypes: quintiles'
)
dev.off()

####################################################
# OS: quartiles: Q1 vs Q4

tiff("./figures/twoGns/twoGnsModel_metabric_os_X3_Q1vsQ4.tiff", width = 3000, height = 2800, res = 300)
ggsurvplot(Fit_sig_metabric_os_twoGnsModel_Q1vsQ4_X3,
           risk.table = FALSE,
           pval = FALSE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           legend.title	= 'Quartiles',
           legend.labs = c('Q1', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           title = 'KLF7+ZNF627 and METABRIC OS by X3 classifier subtypes: Q1 vs Q4'
)
dev.off()

####################################################
## OS: quintiles: Q1 vs Q5
tiff("./figures/twoGns/twoGnsModel_metabric_os_X3_Q1vsQ5.tiff", width = 3000, height = 2800, res = 300)
ggsurvplot(Fit_sig_metabric_os_twoGnsModel_Q1vsQ5_X3,
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
           title = 'KLF7+ZNF627 and METABRIC OS by X3 subtypes: Q1 vs Q5'
)
dev.off()

###############################################################
# RFS: quartiles: all
tiff("./figures/twoGns/twoGnsModel_metabric_rfs_X3_quartiles.tiff", width = 3000, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_rfs_twoGnsModel_quartiles_X3,
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
           title = 'KLF7+ZNF627 and METABRIC RFS by X3 classifier subtypes: quartiles'
)
dev.off()

####################################################
## RFS: quintiles: all
tiff("./figures/twoGns/twoGnsModel_metabric_rfs_X3_quintiles.tiff", width = 3000, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_rfs_twoGnsModel_quintiles_X3,
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
           title = 'KLF7+ZNF627 and METABRIC RFS by X3 subtypes: quintiles'
)
dev.off()

###############################################################
# RFS: quartiles: Q1 vs Q4

tiff("./figures/twoGns/twoGnsModel_metabric_rfs_X3_Q1vsQ4.tiff", width = 3000, height = 2800, res = 300)
ggsurvplot(Fit_sig_metabric_rfs_twoGnsModel_Q1vsQ4_X3,
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
           title = 'KLF7+ZNF627 and METABRIC RFS by X3 classifier subtypes: Q1 vs Q4'
)
dev.off()

####################################################
## RFS: quintiles: Q1 vs Q5
tiff("./figures/twoGns/twoGnsModel_metabric_rfs_X3_Q1vsQ5.tiff", width = 3000, height = 2800, res = 300)
ggsurvplot(Fit_sig_metabric_rfs_twoGnsModel_Q1vsQ5_X3,
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
           title = 'KLF7+ZNF627 and METABRIC RFS by X3 subtypes: Q1 vs Q5'
)
dev.off()



##############################################################################
##############################################################################
##############################################################################
## fit coxph model:

########
## OS

# by probaility

#Fit_sig_metabric_os_coxph_THR25 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_prob_THR25, data = CoxData_metabric)
#summary(Fit_sig_metabric_os_coxph_THR25)

Fit_sig_metabric_os_coxph_twoGnsModel <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_prob_twoGnsModel, data = CoxData_metabric)
summary(Fit_sig_metabric_os_coxph_twoGnsModel)

#Fit_sig_metabric_os_coxph_THR50_2 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_prob_THR50_2, data = CoxData_metabric)
#summary(Fit_sig_metabric_os_coxph_THR50_2)

#png('./figures/logreg/THR25_HR_metabric_os.png', width = 2000, height = 2000, res = 300)
#ggforest(Fit_sig_metabric_coxph_THR25, fontsize = 0.5)
#dev.off()

#png('./figures/logreg/twoGnsModel_HR_metabric_os.png', width = 2000, height = 2000, res = 300)
#ggforest(Fit_sig_metabric_coxph_twoGnsModel, fontsize = 0.5)
#dev.off()

#png('./figures/logreg/THR50_2_HR_metabric_os.png', width = 2000, height = 2000, res = 300)
#ggforest(Fit_sig_metabric_coxph_THR50_2, fontsize = 0.5)
#dev.off()

########
## by quartiles

# make a factor with Q1 (lowest risk) being the reference
#CoxData_metabric$metabric_prob_THR25_quartiles <- factor(CoxData_metabric$metabric_prob_THR25_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_metabric$metabric_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_metabric$metabric_prob_THR25_quartiles))

CoxData_metabric$metabric_prob_twoGnsModel_quartiles <- factor(CoxData_metabric$metabric_prob_twoGnsModel_quartiles, levels = c('1', '2', '3', '4'))
levels(CoxData_metabric$metabric_prob_twoGnsModel_quartiles) <- paste0('Q', levels(CoxData_metabric$metabric_prob_twoGnsModel_quartiles))

#CoxData_metabric$metabric_prob_THR50_2_quartiles <- factor(CoxData_metabric$metabric_prob_THR50_2_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_metabric$metabric_prob_THR50_2_quartiles) <- paste0('Q', levels(CoxData_metabric$metabric_prob_THR50_2_quartiles))

# fit
#Fit_sig_metabric_os_coxph_THR25_quartiles <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR25_quartiles, data = CoxData_metabric)
#summary_metabric_os_coxph_THR25_quartiles <- summary(Fit_sig_metabric_os_coxph_THR25_quartiles)

Fit_sig_metabric_os_coxph_twoGnsModel_quartiles <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quartiles, data = CoxData_metabric)
summary_metabric_os_coxph_twoGnsModel_quartiles <- summary(Fit_sig_metabric_os_coxph_twoGnsModel_quartiles)

#Fit_sig_metabric_os_coxph_THR50_2_quartiles <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR50_2_quartiles, data = CoxData_metabric)
#summary_metabric_os_coxph_THR50_2_quartiles <- summary(Fit_sig_metabric_os_coxph_THR50_2_quartiles)

#summary_list_metabric_os_quartiles <- list(THR25 = summary_metabric_os_coxph_THR25_quartiles, twoGnsModel = summary_metabric_os_coxph_twoGnsModel_quartiles, THR50_2 = summary_metabric_os_coxph_THR50_2_quartiles) 

# get the HR
# HR_list_metabric_os_coxph_quartiles <- lapply(summary_list_metabric_os_quartiles, function(x){
#   HR <- x$conf.int[, 'exp(coef)']
#   lower_95CI <- x$conf.int[, 'lower .95']
#   upper_95CI <- x$conf.int[, 'upper .95']
#   Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
#   Pvalue_logrank_test <- x$sctest['pvalue']
#   Pvalue_wald_test <- x$waldtest['pvalue']
#   data.frame(HR = HR, lower_95CI = lower_95CI, upper_95CI = upper_95CI, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test,
#              Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
# })

# HR_df_metabric_os_coxph_quartiles_THR25 <- as.data.frame(HR_list_metabric_os_coxph_quartiles$THR25)
# HR_df_metabric_os_coxph_quartiles_THR25$quartile <- gsub('.+quartiles', '', rownames(HR_df_metabric_os_coxph_quartiles_THR25))
# 
# HR_df_metabric_os_coxph_quartiles_twoGnsModel <- as.data.frame(HR_list_metabric_os_coxph_quartiles$twoGnsModel)
# HR_df_metabric_os_coxph_quartiles_twoGnsModel$quartile <- gsub('.+quartiles', '', rownames(HR_df_metabric_os_coxph_quartiles_twoGnsModel))
# 
# HR_df_metabric_os_coxph_quartiles_THR50_2 <- as.data.frame(HR_list_metabric_os_coxph_quartiles$THR50_2)
# HR_df_metabric_os_coxph_quartiles_THR50_2$quartile <- gsub('.+quartiles', '', rownames(HR_df_metabric_os_coxph_quartiles_THR50_2))
# 
# # save the results
# write.csv(HR_df_metabric_os_coxph_quartiles_THR25, 'objs/HR/metabric/OS/THR25_quartiles_HR.csv')
# write.csv(HR_df_metabric_os_coxph_quartiles_twoGnsModel, 'objs/HR/metabric/OS/twoGnsModel_quartiles_HR.csv')
# write.csv(HR_df_metabric_os_coxph_quartiles_THR50_2, 'objs/HR/metabric/OS/THR50_2_quartiles_HR.csv')

################
# RFS

# by probaility

#Fit_sig_metabric_RFS_coxph_THR25 <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_prob_THR25, data = CoxData_metabric)
#summary(Fit_sig_metabric_RFS_coxph_THR25)

Fit_sig_metabric_RFS_coxph_twoGnsModel <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_prob_twoGnsModel, data = CoxData_metabric)
summary(Fit_sig_metabric_RFS_coxph_twoGnsModel)

#Fit_sig_metabric_RFS_coxph_THR50_2 <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_prob_THR50_2, data = CoxData_metabric)
#summary(Fit_sig_metabric_RFS_coxph_THR50_2)


########
## by quartiles

# fit
#Fit_sig_metabric_RFS_coxph_THR25_quartiles <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR25_quartiles, data = CoxData_metabric)
#summary_metabric_RFS_coxph_THR25_quartiles <- summary(Fit_sig_metabric_RFS_coxph_THR25_quartiles)

Fit_sig_metabric_RFS_coxph_twoGnsModel_quartiles <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quartiles, data = CoxData_metabric)
summary_metabric_RFS_coxph_twoGnsModel_quartiles <- summary(Fit_sig_metabric_RFS_coxph_twoGnsModel_quartiles)

#Fit_sig_metabric_RFS_coxph_THR50_2_quartiles <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR50_2_quartiles, data = CoxData_metabric)
#summary_metabric_RFS_coxph_THR50_2_quartiles <- summary(Fit_sig_metabric_RFS_coxph_THR50_2_quartiles)

#summary_list_metabric_RFS_quartiles <- list(THR25 = summary_metabric_RFS_coxph_THR25_quartiles, twoGnsModel = summary_metabric_RFS_coxph_twoGnsModel_quartiles, THR50_2 = summary_metabric_RFS_coxph_THR50_2_quartiles) 

# get the HR
# HR_list_metabric_RFS_coxph_quartiles <- lapply(summary_list_metabric_RFS_quartiles, function(x){
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
# HR_df_metabric_RFS_coxph_quartiles_THR25 <- as.data.frame(HR_list_metabric_RFS_coxph_quartiles$THR25)
# HR_df_metabric_RFS_coxph_quartiles_THR25$quartile <- gsub('.+quartiles', '', rownames(HR_df_metabric_RFS_coxph_quartiles_THR25))
# 
# HR_df_metabric_RFS_coxph_quartiles_twoGnsModel <- as.data.frame(HR_list_metabric_RFS_coxph_quartiles$twoGnsModel)
# HR_df_metabric_RFS_coxph_quartiles_twoGnsModel$quartile <- gsub('.+quartiles', '', rownames(HR_df_metabric_RFS_coxph_quartiles_twoGnsModel))
# 
# HR_df_metabric_RFS_coxph_quartiles_THR50_2 <- as.data.frame(HR_list_metabric_RFS_coxph_quartiles$THR50_2)
# HR_df_metabric_RFS_coxph_quartiles_THR50_2$quartile <- gsub('.+quartiles', '', rownames(HR_df_metabric_RFS_coxph_quartiles_THR50_2))
# 
# # save the results
# write.csv(HR_df_metabric_RFS_coxph_quartiles_THR25, 'objs/HR/metabric/RFS/THR25_quartiles_HR.csv')
# write.csv(HR_df_metabric_RFS_coxph_quartiles_twoGnsModel, 'objs/HR/metabric/RFS/twoGnsModel_quartiles_HR.csv')
# write.csv(HR_df_metabric_RFS_coxph_quartiles_THR50_2, 'objs/HR/metabric/RFS/THR50_2_quartiles_HR.csv')


########################################################################  
########################################################################  
## by clinical groups

########
## by quartiles

## make a factor with Q1 (lowest risk) being the reference

# CoxData_metabric
#CoxData_metabric$metabric_prob_THR25_quartiles <- factor(CoxData_metabric$metabric_prob_THR25_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_metabric$metabric_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_metabric$metabric_prob_THR25_quartiles))

#CoxData_metabric$metabric_prob_twoGnsModel_quartiles <- factor(CoxData_metabric$metabric_prob_twoGnsModel_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_metabric$metabric_prob_twoGnsModel_quartiles) <- paste0('Q', levels(CoxData_metabric$metabric_prob_twoGnsModel_quartiles))

#CoxData_metabric$metabric_prob_THR50_2_quartiles <- factor(CoxData_metabric$metabric_prob_THR50_2_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_metabric$metabric_prob_THR50_2_quartiles) <- paste0('Q', levels(CoxData_metabric$metabric_prob_THR50_2_quartiles))

######
# CoxData_metabric_PAM
#CoxData_metabric_PAM$metabric_prob_THR25_quartiles <- factor(CoxData_metabric_PAM$metabric_prob_THR25_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_metabric_PAM$metabric_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_metabric_PAM$metabric_prob_THR25_quartiles))

CoxData_metabric_PAM$metabric_prob_twoGnsModel_quartiles <- factor(CoxData_metabric_PAM$metabric_prob_twoGnsModel_quartiles, levels = c('1', '2', '3', '4'))
levels(CoxData_metabric_PAM$metabric_prob_twoGnsModel_quartiles) <- paste0('Q', levels(CoxData_metabric_PAM$metabric_prob_twoGnsModel_quartiles))

#CoxData_metabric_PAM$metabric_prob_THR50_2_quartiles <- factor(CoxData_metabric_PAM$metabric_prob_THR50_2_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_metabric_PAM$metabric_prob_THR50_2_quartiles) <- paste0('Q', levels(CoxData_metabric_PAM$metabric_prob_THR50_2_quartiles))

######
# CoxData_metabric_PAM_Q1vsQ4: 
#CoxData_metabric_PAM_Q1vsQ4_THR25$metabric_prob_THR25_quartiles <- factor(CoxData_metabric_PAM_Q1vsQ4_THR25$metabric_prob_THR25_quartiles, levels = c('1', '4'))
#levels(CoxData_metabric_PAM_Q1vsQ4_THR25$metabric_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_metabric_PAM_Q1vsQ4_THR25$metabric_prob_THR25_quartiles))

CoxData_metabric_PAM_Q1vsQ4_twoGnsModel$metabric_prob_twoGnsModel_quartiles <- factor(CoxData_metabric_PAM_Q1vsQ4_twoGnsModel$metabric_prob_twoGnsModel_quartiles, levels = c('1', '4'))
levels(CoxData_metabric_PAM_Q1vsQ4_twoGnsModel$metabric_prob_twoGnsModel_quartiles) <- paste0('Q', levels(CoxData_metabric_PAM_Q1vsQ4_twoGnsModel$metabric_prob_twoGnsModel_quartiles))

#CoxData_metabric_PAM_Q1vsQ4_THR50_2$metabric_prob_THR50_2_quartiles <- factor(CoxData_metabric_PAM_Q1vsQ4_THR50_2$metabric_prob_THR50_2_quartiles, levels = c('1', '4'))
#levels(CoxData_metabric_PAM_Q1vsQ4_THR50_2$metabric_prob_THR50_2_quartiles) <- paste0('Q', levels(CoxData_metabric_PAM_Q1vsQ4_THR50_2$metabric_prob_THR50_2_quartiles))


######
# CoxData_metabric_Q1vsQ4
#CoxData_metabric_Q1vsQ4_THR25$metabric_prob_THR25_quartiles <- factor(CoxData_metabric_Q1vsQ4_THR25$metabric_prob_THR25_quartiles, levels = c('1', '4'))
#levels(CoxData_metabric_Q1vsQ4_THR25$metabric_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_metabric_Q1vsQ4_THR25$metabric_prob_THR25_quartiles))

CoxData_metabric_Q1vsQ4_twoGnsModel$metabric_prob_twoGnsModel_quartiles <- factor(CoxData_metabric_Q1vsQ4_twoGnsModel$metabric_prob_twoGnsModel_quartiles, levels = c('1', '4'))
levels(CoxData_metabric_Q1vsQ4_twoGnsModel$metabric_prob_twoGnsModel_quartiles) <- paste0('Q', levels(CoxData_metabric_Q1vsQ4_twoGnsModel$metabric_prob_twoGnsModel_quartiles))

#CoxData_metabric_Q1vsQ4_THR50_2$metabric_prob_THR50_2_quartiles <- factor(CoxData_metabric_Q1vsQ4_THR50_2$metabric_prob_THR50_2_quartiles, levels = c('1', '4'))
#levels(CoxData_metabric_Q1vsQ4_THR50_2$metabric_prob_THR50_2_quartiles) <- paste0('Q', levels(CoxData_metabric_Q1vsQ4_THR50_2$metabric_prob_THR50_2_quartiles))


#CoxData_metabric_PAM$Pam50...Claudin.low.subtype <- factor(CoxData_metabric_PAM$Pam50...Claudin.low.subtype)

##############################################
## fit

## PAM50

# os: quartiles: all
# fit on each PAM50 subtype
# lapply(split(CoxData_metabric_PAM, CoxData_metabric_PAM$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR25_quartiles, data = x)))

lapply(split(CoxData_metabric_PAM, CoxData_metabric_PAM$Pam50...Claudin.low.subtype),
       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quartiles, data = x)))

#lapply(split(CoxData_metabric_PAM, CoxData_metabric_PAM$Pam50...Claudin.low.subtype),
#       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR50_2_quartiles, data = x)))

################################
# os: quartiles: Q1 vs Q4

# fit on each PAM50 subtype

# lapply(split(CoxData_metabric_PAM_Q1vsQ4_THR25, CoxData_metabric_PAM_Q1vsQ4_THR25$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR25_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_twoGnsModel, CoxData_metabric_PAM_Q1vsQ4_twoGnsModel$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_THR50_2, CoxData_metabric_PAM_Q1vsQ4_THR50_2$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR50_2_quartiles, data = x)))

################################
# rfs: quartiles: all
# fit on each PAM50 subtype
#lapply(split(CoxData_metabric_PAM, CoxData_metabric_PAM$Pam50...Claudin.low.subtype),
#       function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR25_quartiles, data = x)))

lapply(split(CoxData_metabric_PAM, CoxData_metabric_PAM$Pam50...Claudin.low.subtype),
       function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quartiles, data = x)))

#lapply(split(CoxData_metabric_PAM, CoxData_metabric_PAM$Pam50...Claudin.low.subtype),
#       function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR50_2_quartiles, data = x)))

################################
# rfs: quartiles: Q1 vs Q4

# fit on each PAM50 subtype

# lapply(split(CoxData_metabric_PAM_Q1vsQ4_THR25, CoxData_metabric_PAM_Q1vsQ4_THR25$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR25_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_twoGnsModel, CoxData_metabric_PAM_Q1vsQ4_twoGnsModel$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_THR50_2, CoxData_metabric_PAM_Q1vsQ4_THR50_2$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR50_2_quartiles, data = x)))



##############################################
## fit

## X3 classifier

# os: quartiles: all
# fit on each X3 subtypes
#lapply(split(CoxData_metabric, CoxData_metabric$X3.Gene.classifier.subtype),
#       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR25_quartiles, data = x)))

lapply(split(CoxData_metabric, CoxData_metabric$X3.Gene.classifier.subtype),
       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quartiles, data = x)))

#lapply(split(CoxData_metabric, CoxData_metabric$X3.Gene.classifier.subtype),
#       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR50_2_quartiles, data = x)))

################################
# os: quartiles: Q1 vs Q4

# fit on each X3 subtypes

# lapply(split(CoxData_metabric_PAM_Q1vsQ4_THR25, CoxData_metabric_PAM_Q1vsQ4_THR25$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR25_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_twoGnsModel, CoxData_metabric_PAM_Q1vsQ4_twoGnsModel$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_twoGnsModel_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_THR50_2, CoxData_metabric_PAM_Q1vsQ4_THR50_2$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR50_2_quartiles, data = x)))

################################
# rfs: quartiles: all
# fit on each X3 subtypes
#lapply(split(CoxData_metabric, CoxData_metabric$X3.Gene.classifier.subtype),
#       function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR25_quartiles, data = x)))

lapply(split(CoxData_metabric, CoxData_metabric$X3.Gene.classifier.subtype),
       function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quartiles, data = x)))

#lapply(split(CoxData_metabric, CoxData_metabric$X3.Gene.classifier.subtype),
#       function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR50_2_quartiles, data = x)))

################################
# rfs: quartiles: Q1 vs Q4

# fit on each X3 subtypes

# lapply(split(CoxData_metabric_PAM_Q1vsQ4_THR25, CoxData_metabric_PAM_Q1vsQ4_THR25$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR25_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_twoGnsModel, CoxData_metabric_PAM_Q1vsQ4_twoGnsModel$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_twoGnsModel_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_THR50_2, CoxData_metabric_PAM_Q1vsQ4_THR50_2$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR50_2_quartiles, data = x)))






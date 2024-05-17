
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
library(sva)

#############################
E6 <- c('LMNB2', 'CDC20', 'KIF2C', 'FAM64A', 'KIF4A', 'TPX2')

#############################################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')
##############################################


###############################################
## BATCH correction
###############################################
# Combine the expression data into a single matrix
commonGns <- intersect(rownames(Expr_tcga), rownames(Expr_metabric))
Expr_tcga <- Expr_tcga[commonGns, ]
Expr_metabric <- Expr_metabric[commonGns, ]

combined_expr <- cbind(Expr_tcga, Expr_metabric)

# Create a batch vector indicating the dataset origin for each column in the combined matrix
batch <- c(rep(1, ncol(Expr_tcga)), rep(2, ncol(Expr_metabric)))

# Perform batch correction
corrected_expr <- ComBat(dat = combined_expr, batch = batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)

# After correction, you can split the corrected expression matrix back into TCGA and METABRIC parts if necessary
corrected_Expr_tcga <- corrected_expr[, 1:ncol(Expr_tcga)]
corrected_Expr_metabric <- corrected_expr[, (ncol(Expr_tcga) + 1):ncol(corrected_expr)]


##############
# filter the signatures to include only the genes present in the expr matrices
E6_fil <- E6[E6 %in% rownames(corrected_Expr_tcga) & E6 %in% rownames(corrected_Expr_metabric)]

setdiff(E6, E6_fil)

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

corrected_Expr_tcga <- corrected_Expr_tcga[, colnames(corrected_Expr_tcga) %in% rownames(Pheno_tcga)]

##############################################
### combine in 1 dataset
##############################################

#############
# fix the survival information
Pheno_tcga$Disease.Free.Status <- gsub("\\:.+", "", Pheno_tcga$Disease.Free.Status)
Pheno_tcga$Overall.Survival.Status <- gsub("\\:.+", "", Pheno_tcga$Overall.Survival.Status)

Pheno_tcga$Disease.Free.Status <- as.numeric(Pheno_tcga$Disease.Free.Status)
Pheno_tcga$Overall.Survival.Status <- as.numeric(Pheno_tcga$Overall.Survival.Status)

table(Pheno_tcga$Disease.Free.Status)
table(Pheno_tcga$Overall.Survival.Status)

table(Pheno_tcga$Disease.Free.Status)
table(Pheno_tcga$Overall.Survival.Status)

Pheno_tcga$Disease.Free..Months. <- as.numeric(Pheno_tcga$Disease.Free..Months.)
Pheno_tcga$Overall.Survival..Months. <- as.numeric(Pheno_tcga$Overall.Survival..Months.)


group_tcga <- as.factor(Pheno_tcga$Overall.Survival.Status)
table(group_tcga)
Data_tcga <- as.data.frame(cbind(t(corrected_Expr_tcga), group_tcga))
Data_tcga$group_tcga <- as.factor(Data_tcga$group_tcga)
levels(Data_tcga$group_tcga) <- c('0', '1')
table(Data_tcga$group_tcga)
colnames(Data_tcga)[colnames(Data_tcga) %in% c('group_tcga')] <- c('os')
table(Data_tcga$os)

group_tcga <- Data_tcga$os


####################################
# get ER+HER-
####################################
table(Pheno_tcga$ER.Status.By.IHC)
table(Pheno_tcga$IHC.HER2)
table(Pheno_tcga$Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code)

table(Pheno_tcga$ER.status.measured.by.IHC, Pheno_tcga$Lymph.nodes.examined.positive)


Pheno_tcga_ERposHER2neg <- Pheno_tcga %>%
  dplyr::filter(ER.Status.By.IHC == 'Positive',
                IHC.HER2 == 'Negative') 

Expr_tcga_ERposHER2neg <- corrected_Expr_tcga[, rownames(Pheno_tcga_ERposHER2neg)]
all(colnames(Expr_tcga_ERposHER2neg) == rownames(Pheno_tcga_ERposHER2neg))

########################################
## Make a dataframe with survival time and event
########################################
os_event <- as.factor(Pheno_tcga_ERposHER2neg$Overall.Survival.Status)
os_time <- Pheno_tcga_ERposHER2neg$Overall.Survival..Months.
summary(os_time)

rfs_event <- as.factor(Pheno_tcga_ERposHER2neg$Disease.Free.Status)
rfs_time <- Pheno_tcga_ERposHER2neg$Disease.Free..Months.
summary(rfs_time)


Data_tcga_ERposHER2neg <- as.data.frame(t(Expr_tcga_ERposHER2neg))
Data_tcga_ERposHER2neg$os_event <- os_event
Data_tcga_ERposHER2neg$os_time <- os_time
summary(Data_tcga_ERposHER2neg$os_time)

Data_tcga_ERposHER2neg$rfs_event <- rfs_event
Data_tcga_ERposHER2neg$rfs_time <- rfs_time
summary(Data_tcga_ERposHER2neg$rfs_time)

all(rownames(Data_tcga_ERposHER2neg) == rownames(Pheno_tcga_ERposHER2neg))

levels(Data_tcga_ERposHER2neg$os_event) <- c('low.risk', 'high.risk')
levels(Data_tcga_ERposHER2neg$rfs_event) <- c('low.risk', 'high.risk')


################################################################
### METABRIC
################################################################
rownames(corrected_Expr_metabric) <- gsub('-', '', rownames(corrected_Expr_metabric))
Data_metabric <- as.data.frame(cbind(t(corrected_Expr_metabric), group_metabric))
Data_metabric$group_metabric <- as.factor(Data_metabric$group_metabric)
levels(Data_metabric$group_metabric) <- c('low.risk', 'high.risk')
colnames(Data_metabric)[colnames(Data_metabric) %in% c('group_metabric')] <- c('os')

##################
# get ER+HER-LN-
##################
table(Pheno_metabric$ER.status.measured.by.IHC)
table(Pheno_metabric$HER2.Status)

table(Pheno_metabric$ER.status.measured.by.IHC, Pheno_metabric$Lymph.nodes.examined.positive)


Pheno_metabric_ERposHER2negLNneg <- Pheno_metabric %>%
  dplyr::filter(ER.status.measured.by.IHC == 'Positve',
                HER2.Status == 'Negative', 
                Lymph.nodes.examined.positive == 0)  

Expr_metabric_ERposHER2negLNneg <- corrected_Expr_metabric[, rownames(Pheno_metabric_ERposHER2negLNneg)]
all(colnames(Expr_metabric_ERposHER2negLNneg) == rownames(Pheno_metabric_ERposHER2negLNneg))

########################################
## Make a dataframe with survival time and event
########################################
os_event <- as.factor(Pheno_metabric_ERposHER2negLNneg$Overall.Survival.Status)
os_time <- Pheno_metabric_ERposHER2negLNneg$Overall.Survival..Months.
summary(os_time)

rfs_event <- as.factor(Pheno_metabric_ERposHER2negLNneg$Relapse.Free.Status)
rfs_time <- Pheno_metabric_ERposHER2negLNneg$Relapse.Free.Status..Months.
summary(rfs_time)


Data_metabric_ERposHER2negLNneg <- as.data.frame(t(Expr_metabric_ERposHER2negLNneg))
Data_metabric_ERposHER2negLNneg$os_event <- os_event
Data_metabric_ERposHER2negLNneg$os_time <- os_time
summary(Data_metabric_ERposHER2negLNneg$os_time)

Data_metabric_ERposHER2negLNneg$rfs_event <- rfs_event
Data_metabric_ERposHER2negLNneg$rfs_time <- rfs_time
summary(Data_metabric_ERposHER2negLNneg$rfs_time)

all(rownames(Data_metabric_ERposHER2negLNneg) == rownames(Pheno_metabric_ERposHER2negLNneg))


levels(Data_metabric_ERposHER2negLNneg$os_event) <- c('low.risk', 'high.risk')
levels(Data_metabric_ERposHER2negLNneg$rfs_event) <- c('low.risk', 'high.risk')

#############################################################################################################
##############################################################################################################
## Train the model

#########################
## Get the best parameters

# control <- trainControl(method="repeatedcv", number=5, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = T)
# 
# # 5-fold cross validation repeated 5 times (to find the best parameters)
# set.seed(333)
# fit_init <- train(as.formula((paste("os_event ~", paste(E6_fil, collapse = "+")))), data=Data_tcga_ERposHER2neg, method="svmPoly", trControl=control, tuneLength = 5, metric = "ROC")
# fit_init

###########################
# Use the best parameters in the final model
Grid <- expand.grid(degree = 3, scale = 1, C = 0.25)

set.seed(333)

E6_svm_model <- train(as.formula((paste("os_event ~", paste(E6_fil, collapse = "+")))), data=Data_tcga_ERposHER2neg, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
E6_svm_model

save(E6_svm_model, file = 'objs/E6_svm_model.rda')

###########################################################################
# Make predictions
############################################################################

######################
# Training set: TCGA
######################

tcga_prob <- E6_svm_model %>% predict(Data_tcga_ERposHER2neg , type = "prob")

### Threshold
thr <- coords(roc(Data_tcga_ERposHER2neg$os_event, tcga_prob[,2], levels = c("low.risk", "high.risk"), direction = "<"), "best")["threshold"]
thr

### ROC Curve
roc_tcga <- roc(Data_tcga_ERposHER2neg$os_event, tcga_prob[,2], plot = F, print.thres=thr$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("low.risk", "high.risk"), direction = "<", col="blue", lwd=2, grid=TRUE)
roc_tcga

### Get predictions based on best threshold from ROC curve
tcga_classes <- ifelse(tcga_prob[,2] >= thr$threshold, "high.risk", "low.risk")
table(tcga_classes)
tcga_classes <- factor(tcga_classes, levels = c('low.risk', 'high.risk'))

###
cm_tcga <- confusionMatrix(tcga_classes, Data_tcga_ERposHER2neg$os_event, positive = "high.risk", mode = "everything")
cm_tcga

## MCC
mcc_tcga <- mltools::mcc(pred = tcga_classes, actuals = Data_tcga_ERposHER2neg$os_event)
mcc_tcga

######################
# Testing set: Metabric
######################

metabric_prob <- E6_svm_model %>% predict(Data_metabric_ERposHER2negLNneg, type = "prob")

### ROC Curve
roc_metabric <- roc(Data_metabric_ERposHER2negLNneg$os_event, metabric_prob[,2], plot = F, print.thres=thr$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("low.risk", "high.risk"), direction = "<", col="blue", lwd=2, grid=TRUE)
roc_metabric

### Get predictions based on best threshold from ROC curve
metabric_classes <- ifelse(metabric_prob[,2] >= thr$threshold, "high.risk", "low.risk")
table(metabric_classes)
metabric_classes <- factor(metabric_classes, levels = c('low.risk', 'high.risk'))

###
cm_metabric <- confusionMatrix(metabric_classes, Data_metabric_ERposHER2negLNneg$os_event, positive = "high.risk", mode = "everything")
cm_metabric

## MCC
mcc_metabric <- mltools::mcc(pred = metabric_classes, actuals = Data_metabric_ERposHER2negLNneg$os_event)
mcc_metabric

##############################################################################
## Prepare for Survival analysis
##############################################################################

###########################
# TCGA
###########################

## Keep only the relevant information (Metastasis Event and Time)
survival_tcga <- cbind(Pheno_tcga_ERposHER2neg[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Disease.Free.Status", "Disease.Free..Months.", 
                                                                 "Subtype", "ER.Status.By.IHC", "PR.status.by.ihc")], 
                                     tcga_prob, tcga_classes)

# create a merged pdata and Z-scores object
survival_tcga <- data.frame(survival_tcga)

# divide the probabilities into quartiles
survival_tcga <- survival_tcga %>%
  mutate(quartiles = ntile(high.risk, 4), 
         quintiles = ntile(high.risk, 5),
         tertiles = ntile(high.risk, 3)
  )

###########################
# Metabric
###########################

## Keep only the relevant information (Metastasis Event and Time)
survival_metabric <- cbind(Pheno_metabric_ERposHER2negLNneg[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Relapse.Free.Status", "Relapse.Free.Status..Months.", 
                                                   "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC", "PR.Status")], 
                       metabric_prob, metabric_classes)

# create a merged pdata and Z-scores object
survival_metabric <- data.frame(survival_metabric)

# divide the probabilities into quartiles
survival_metabric <- survival_metabric %>%
  mutate(quartiles = ntile(high.risk, 4), 
         quintiles = ntile(high.risk, 5),
         tertiles = ntile(high.risk, 3)
  )
###########################################################################
# new approach to quartiles using the average expression of all genes in the signature
###########################################################################

########################################################
# using the GA algorithm to divide the probabilities into three groups

library(GA)

# Define the fitness function
fitness_function <- function(cutoffs, prognostic_index, surv_time, surv_event) {
  # Create the groups based on the cut-offs
  groups <- findInterval(prognostic_index, c(-Inf, sort(cutoffs), Inf))
  
  # If there's only one unique group, return a very low fitness
  if (length(unique(groups)) < 2) return(1e8)
  
  # Compute the survival difference
  surv_diff_result <- survdiff(Surv(surv_time, surv_event) ~ groups)
  
  # Check if p-value is present, if not return a low fitness value
  if (is.null(surv_diff_result$pval)) return(1e8)
  
  # Return the negative p-value
  return(-surv_diff_result$pval)
}


###########################################################################
# Run the Genetic Algorithm to find the optimal cut-offs on prediction probabilities
###########################################################################


##############################################################################
# TCGA OS
##############################################################################
tcga_result_os_tertiles <- ga(type = "real-valued", 
                                           fitness = function(cutoffs) fitness_function(cutoffs, survival_tcga$high.risk, 
                                                                                        survival_tcga$`Overall.Survival..Months.`,
                                                                                        survival_tcga$`Overall.Survival.Status`),
                                           min = rep(min(survival_tcga$high.risk), 2), 
                                           max = rep(max(survival_tcga$high.risk), 2), 
                                           names = c("Cutoff 1", "Cutoff 2"),
                                           popSize = 100, maxiter = 200, run = 100, seed = 123456)

# Extract the optimal cut-offs
tcga_optimal_cutoffs_os_tertiles <- tcga_result_os_tertiles@solution

tcga_optimal_cutoffs_os_tertiles

tcga_best_p_value_os_tertiles <- Inf
tcga_best_cutoffs_os_tertiles <- NULL

for (i in 1:nrow(tcga_optimal_cutoffs_os_tertiles)) {
  tcga_groups_os_tertiles <- findInterval(survival_tcga$high.risk, c(-Inf, sort(tcga_optimal_cutoffs_os_tertiles[i, ]), Inf))
  tcga_surv_diff_result_os_tertiles <- survdiff(Surv(survival_tcga$`Overall.Survival..Months.`, survival_tcga$`Overall.Survival.Status`) ~ tcga_groups_os_tertiles)
  
  if (!is.null(tcga_surv_diff_result_os_tertiles$pval) && tcga_surv_diff_result_os_tertiles$pval < tcga_best_p_value_os_tertiles) {
    tcga_best_p_value_os_tertiles <- tcga_surv_diff_result_os_tertiles$pval
    tcga_best_cutoffs_os_tertiles <- tcga_optimal_cutoffs_os_tertiles[i, ]
  }
}

##############################################################################
# TCGA RFS
##############################################################################
tcga_result_rfs_tertiles <- ga(type = "real-valued", 
                              fitness = function(cutoffs) fitness_function(cutoffs, survival_tcga$high.risk, 
                                                                           survival_tcga$Disease.Free..Months.,
                                                                           survival_tcga$Disease.Free.Status),
                              min = rep(min(survival_tcga$high.risk), 2), 
                              max = rep(max(survival_tcga$high.risk), 2), 
                              names = c("Cutoff 1", "Cutoff 2"),
                              popSize = 100, maxiter = 200, run = 100, seed = 123456)

# Extract the optimal cut-offs
tcga_optimal_cutoffs_rfs_tertiles <- tcga_result_rfs_tertiles@solution

tcga_optimal_cutoffs_rfs_tertiles

tcga_best_p_value_rfs_tertiles <- Inf
tcga_best_cutoffs_rfs_tertiles <- NULL

for (i in 1:nrow(tcga_optimal_cutoffs_rfs_tertiles)) {
  tcga_groups_rfs_tertiles <- findInterval(survival_tcga$high.risk, c(-Inf, sort(tcga_optimal_cutoffs_rfs_tertiles[i, ]), Inf))
  tcga_surv_diff_result_rfs_tertiles <- survdiff(Surv(survival_tcga$Disease.Free..Months., survival_tcga$Disease.Free.Status) ~ tcga_groups_rfs_tertiles)
  
  if (!is.null(tcga_surv_diff_result_rfs_tertiles$pval) && tcga_surv_diff_result_rfs_tertiles$pval < tcga_best_p_value_rfs_tertiles) {
    tcga_best_p_value_rfs_tertiles <- tcga_surv_diff_result_rfs_tertiles$pval
    tcga_best_cutoffs_rfs_tertiles <- tcga_optimal_cutoffs_rfs_tertiles[i, ]
  }
}

##############################################################################
# METABRIC OS
##############################################################################
metabric_result_os_tertiles <- ga(type = "real-valued", 
                              fitness = function(cutoffs) fitness_function(cutoffs, survival_metabric$high.risk, 
                                                                           survival_metabric$`Overall.Survival..Months.`,
                                                                           survival_metabric$`Overall.Survival.Status`),
                              min = rep(min(survival_metabric$high.risk), 2), 
                              max = rep(max(survival_metabric$high.risk), 2), 
                              names = c("Cutoff 1", "Cutoff 2"),
                              popSize = 100, maxiter = 200, run = 100, seed = 123456)

# Extract the optimal cut-offs
metabric_optimal_cutoffs_os_tertiles <- metabric_result_os_tertiles@solution

metabric_optimal_cutoffs_os_tertiles

metabric_best_p_value_os_tertiles <- Inf
metabric_best_cutoffs_os_tertiles <- NULL

for (i in 1:nrow(metabric_optimal_cutoffs_os_tertiles)) {
  metabric_groups_os_tertiles <- findInterval(survival_metabric$high.risk, c(-Inf, sort(metabric_optimal_cutoffs_os_tertiles[i, ]), Inf))
  metabric_surv_diff_result_os_tertiles <- survdiff(Surv(survival_metabric$`Overall.Survival..Months.`, survival_metabric$`Overall.Survival.Status`) ~ metabric_groups_os_tertiles)
  
  if (!is.null(metabric_surv_diff_result_os_tertiles$pval) && metabric_surv_diff_result_os_tertiles$pval < metabric_best_p_value_os_tertiles) {
    metabric_best_p_value_os_tertiles <- metabric_surv_diff_result_os_tertiles$pval
    metabric_best_cutoffs_os_tertiles <- metabric_optimal_cutoffs_os_tertiles[i, ]
  }
}

##############################################################################
# METABRIC RFS
##############################################################################
metabric_result_rfs_tertiles <- ga(type = "real-valued", 
                               fitness = function(cutoffs) fitness_function(cutoffs, survival_metabric$high.risk, 
                                                                            survival_metabric$Relapse.Free.Status..Months.,
                                                                            survival_metabric$Relapse.Free.Status),
                               min = rep(min(survival_metabric$high.risk), 2), 
                               max = rep(max(survival_metabric$high.risk), 2), 
                               names = c("Cutoff 1", "Cutoff 2"),
                               popSize = 100, maxiter = 200, run = 100, seed = 123456)

# Extract the optimal cut-offs
metabric_optimal_cutoffs_rfs_tertiles <- metabric_result_rfs_tertiles@solution

metabric_optimal_cutoffs_rfs_tertiles

metabric_best_p_value_rfs_tertiles <- Inf
metabric_best_cutoffs_rfs_tertiles <- NULL

for (i in 1:nrow(metabric_optimal_cutoffs_rfs_tertiles)) {
  metabric_groups_rfs_tertiles <- findInterval(survival_metabric$high.risk, c(-Inf, sort(metabric_optimal_cutoffs_rfs_tertiles[i, ]), Inf))
  metabric_surv_diff_result_rfs_tertiles <- survdiff(Surv(survival_metabric$Relapse.Free.Status..Months., survival_metabric$Relapse.Free.Status) ~ metabric_groups_rfs_tertiles)
  
  if (!is.null(metabric_surv_diff_result_rfs_tertiles$pval) && metabric_surv_diff_result_rfs_tertiles$pval < metabric_best_p_value_rfs_tertiles) {
    metabric_best_p_value_rfs_tertiles <- metabric_surv_diff_result_rfs_tertiles$pval
    metabric_best_cutoffs_rfs_tertiles <- metabric_optimal_cutoffs_rfs_tertiles[i, ]
  }
}

############################################################################################################
# Fit the survival curves
############################################################################################################

# Create the survival object
tcga_os_surv_obj <- Surv(survival_tcga$`Overall.Survival..Months.`, survival_tcga$`Overall.Survival.Status`)
tcga_rfs_surv_obj <- Surv(survival_tcga$Disease.Free..Months., survival_tcga$Disease.Free.Status)

metabric_os_surv_obj <- Surv(survival_metabric$Overall.Survival..Months., survival_metabric$Overall.Survival.Status)
metabric_rfs_surv_obj <- Surv(survival_metabric$Relapse.Free.Status..Months., survival_metabric$Relapse.Free.Status)

####################################################
# two predicted classes
####################################################

##########################
## TCGA OS
##########################

tcga_fit_os_two_classes <- survfit(tcga_os_surv_obj ~ tcga_classes)

tiff("./figures/E6_tcga_prognostic/E6_tcga_svm_os_two_predClasses.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(tcga_fit_os_two_classes, 
           data = survival_tcga, 
           xlim = c(0,240),
           break.x.by = 40,
           risk.table = TRUE, pval = TRUE,
           ggtheme = theme_survminer(base_size = 16, font.x = c(16, 'bold.italic', 'black'), font.y = c(16, 'bold.italic', 'black'), font.tickslab = c(16, 'plain', 'black'), font.legend = c(16, 'bold', 'black')),
           legend.labs = c('low-risk', 'high-risk')
)
dev.off()

# COX
tcga_cox_os_two_classes <- coxph(tcga_os_surv_obj ~ as.factor(tcga_classes), data = survival_tcga)
summary(tcga_cox_os_two_classes)

##########################
## TCGA DFS
##########################

tcga_fit_rfs_two_classes <- survfit(tcga_rfs_surv_obj ~ tcga_classes)

tiff("./figures/E6_tcga_prognostic/E6_tcga_svm_rfs_two_predClasses.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(tcga_fit_rfs_two_classes, 
           data = survival_tcga, 
           xlim = c(0,240),
           break.x.by = 40,
           risk.table = TRUE, pval = TRUE,
           ggtheme = theme_survminer(base_size = 16, font.x = c(16, 'bold.italic', 'black'), font.y = c(16, 'bold.italic', 'black'), font.tickslab = c(16, 'plain', 'black'), font.legend = c(16, 'bold', 'black')),
           legend.labs = c('low-risk', 'high-risk')
)
dev.off()

# COX
tcga_cox_rfs_two_classes <- coxph(tcga_rfs_surv_obj ~ as.factor(tcga_classes), data = survival_tcga)
summary(tcga_cox_rfs_two_classes)

##########################
## METABRIC OS
##########################

metabric_fit_os_two_classes <- survfit(metabric_os_surv_obj ~ metabric_classes)

tiff("./figures/E6_metabric_prognostic/E6_metabric_svm_os_two_predClasses.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(metabric_fit_os_two_classes, 
           data = survival_metabric, 
           xlim = c(0,240),
           break.x.by = 40,
           risk.table = TRUE, pval = TRUE,
           ggtheme = theme_survminer(base_size = 16, font.x = c(16, 'bold.italic', 'black'), font.y = c(16, 'bold.italic', 'black'), font.tickslab = c(16, 'plain', 'black'), font.legend = c(16, 'bold', 'black')),
           legend.labs = c('low-risk', 'high-risk')
)
dev.off()

# COX
metabric_cox_os_two_classes <- coxph(metabric_os_surv_obj ~ as.factor(metabric_classes), data = survival_metabric)
summary(metabric_cox_os_two_classes)

##########################
## METABRIC FFS
##########################

metabric_fit_rfs_two_classes <- survfit(metabric_rfs_surv_obj ~ metabric_classes)

tiff("./figures/E6_metabric_prognostic/E6_metabric_svm_rfs_two_predClasses.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(metabric_fit_rfs_two_classes, 
           data = survival_metabric, 
           xlim = c(0,240),
           break.x.by = 40,
           risk.table = TRUE, pval = TRUE,
           ggtheme = theme_survminer(base_size = 16, font.x = c(16, 'bold.italic', 'black'), font.y = c(16, 'bold.italic', 'black'), font.tickslab = c(16, 'plain', 'black'), font.legend = c(16, 'bold', 'black')),
           legend.labs = c('low-risk', 'high-risk')
)
dev.off()

# COX
metabric_cox_rfs_two_classes <- coxph(metabric_rfs_surv_obj ~ as.factor(metabric_classes), data = survival_metabric)
summary(metabric_cox_rfs_two_classes)

####################################################
# Best cutoff tertiles
####################################################

##########################
## TCGA OS
##########################

tcga_fit_os_bestcutoff_tertiles <- survfit(os_surv_obj ~ groups_os_tertiles)
tiff("./figures/E6_tcga_prognostic/E6_tcga_os_bestCutoff_tertiles_svm.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(fit_os_bestcutoff_tertiles, 
           data = CoxData_tcga_ERposHER2neg, 
           xlim = c(0,240),
           break.x.by = 40,
           risk.table = TRUE, pval = TRUE,
           ggtheme = theme_survminer(base_size = 16, font.x = c(16, 'bold.italic', 'black'), font.y = c(16, 'bold.italic', 'black'), font.tickslab = c(16, 'plain', 'black'), font.legend = c(16, 'bold', 'black')),
           legend.labs = c('Q1', 'Q2', 'Q3')
)
dev.off()

# Fit the Cox proportional hazards model using best cut-offs
cox_fit_os_bestcutoff_tertiles <- coxph(os_surv_obj ~ as.factor(groups_os_tertiles), data = CoxData_tcga_ERposHER2neg)
summary(cox_fit_os_bestcutoff_tertiles)

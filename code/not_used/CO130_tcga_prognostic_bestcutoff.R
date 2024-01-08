
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


########################
# OS
########################

result_prob_os <- ga(type = "real-valued", 
                     fitness = function(cutoffs) fitness_function(cutoffs, CoxData_tcga$Train_prob_CO130, 
                                                                  CoxData_tcga$Overall.Survival..Months.,
                                                                  CoxData_tcga$Overall.Survival.Status),
                     min = rep(min(CoxData_tcga$Train_prob_CO130), 2), 
                     max = rep(max(CoxData_tcga$Train_prob_CO130), 2), 
                     names = c("Cutoff 1", "Cutoff 2"),
                     popSize = 50, maxiter = 100, run = 100, seed = 123456)

# Extract the optimal cut-offs
optimal_cutoffs_prob_os <- result_prob_os@solution

optimal_cutoffs_prob_os

best_p_value_os <- Inf
best_cutoffs_os <- NULL

for (i in 1:nrow(optimal_cutoffs_prob_os)) {
  groups_prob_os <- findInterval(CoxData_tcga$Train_prob_CO130, c(-Inf, sort(optimal_cutoffs_prob_os[i, ]), Inf))
  surv_diff_result_prob_os <- survdiff(Surv(CoxData_tcga$`Overall.Survival..Months.`, CoxData_tcga$`Overall.Survival.Status`) ~ groups_prob_os)
  
  if (!is.null(surv_diff_result_prob_os$pval) && surv_diff_result_prob_os$pval < best_p_value_os) {
    best_p_value_os <- surv_diff_result_prob_os$pval
    best_cutoffs_os <- optimal_cutoffs_prob_os[i, ]
  }
}

# Create the survival object
os_surv_obj <- Surv(CoxData_tcga$`Overall.Survival..Months.`, CoxData_tcga$`Overall.Survival.Status`)

# Fit the survival curves
fit_os_bestcutoff_tertiles_prob <- survfit(os_surv_obj ~ groups_prob_os)
tiff("./figures/CO130_tcga_prognostic/CO130_tcga_os_bestCutoff_tertiles_prob.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(fit_os_bestcutoff_tertiles_prob, 
           data = CoxData_tcga, 
           risk.table = TRUE, pval = TRUE,
           xlim = c(0,120),
           break.x.by = 20,
           ggtheme = theme_survminer(base_size = 16, font.x = c(16, 'bold.italic', 'black'), font.y = c(16, 'bold.italic', 'black'), font.tickslab = c(16, 'plain', 'black'), font.legend = c(16, 'bold', 'black')),
           legend.labs = c('Q1', 'Q2', 'Q3')
)
dev.off()


# Fit the Cox proportional hazards model using best cut-offs
cox_fit_os_bestcutoff_tertiles_prob <- coxph(os_surv_obj ~ as.factor(groups_prob_os), data = CoxData_tcga)
summary(cox_fit_os_bestcutoff_tertiles_prob)


########################
# DFS
########################

result_prob_dfs <- ga(type = "real-valued", 
                      fitness = function(cutoffs) fitness_function(cutoffs, CoxData_tcga$Train_prob_CO130, 
                                                                   CoxData_tcga$Disease.Free..Months.,
                                                                   CoxData_tcga$Disease.Free.Status),
                      min = rep(min(CoxData_tcga$Train_prob_CO130), 2), 
                      max = rep(max(CoxData_tcga$Train_prob_CO130), 2), 
                      names = c("Cutoff 1", "Cutoff 2"),
                      popSize = 50, maxiter = 100, run = 100, seed = 123456)

# Extract the optimal cut-offs
optimal_cutoffs_prob_dfs <- result_prob_dfs@solution

optimal_cutoffs_prob_dfs

best_p_value_dfs <- Inf
best_cutoffs_dfs <- NULL

for (i in 1:nrow(optimal_cutoffs_prob_dfs)) {
  groups_prob_dfs <- findInterval(CoxData_tcga$Train_prob_CO130, c(-Inf, sort(optimal_cutoffs_prob_dfs[i, ]), Inf))
  surv_diff_result_prob_dfs <- survdiff(Surv(CoxData_tcga$Disease.Free..Months., CoxData_tcga$Disease.Free.Status) ~ groups_prob_dfs)
  
  if (!is.null(surv_diff_result_prob_dfs$pval) && surv_diff_result_prob_dfs$pval < best_p_value_dfs) {
    best_p_value_dfs <- surv_diff_result_prob_dfs$pval
    best_cutoffs_dfs <- optimal_cutoffs_prob_dfs[i, ]
  }
}

# Create the survival object
dfs_surv_obj <- Surv(CoxData_tcga$Disease.Free..Months., CoxData_tcga$Disease.Free.Status)

# Fit the survival curves
fit_dfs_bestcutoff_tertiles_prob <- survfit(dfs_surv_obj ~ groups_prob_dfs)


# 10yrs
tiff("./figures/CO130_tcga_prognostic/CO130_tcga_dfs_bestCutoff_tertiles_prob_10yrs.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(fit_dfs_bestcutoff_tertiles_prob, 
           data = CoxData_tcga, 
           risk.table = TRUE, pval = TRUE,
           xlim = c(0,120),
           break.x.by = 20,
           ggtheme = theme_survminer(base_size = 16, font.x = c(16, 'bold.italic', 'black'), font.y = c(16, 'bold.italic', 'black'), font.tickslab = c(16, 'plain', 'black'), font.legend = c(16, 'bold', 'black')),
           legend.labs = c('Q1', 'Q2', 'Q3')
)
dev.off()


# Fit the Cox proportional hazards model using best cut-offs
cox_fit_dfs_bestcutoff_tertiles_prob <- coxph(dfs_surv_obj ~ as.factor(groups_prob_dfs), data = CoxData_tcga)
summary(cox_fit_dfs_bestcutoff_tertiles_prob)



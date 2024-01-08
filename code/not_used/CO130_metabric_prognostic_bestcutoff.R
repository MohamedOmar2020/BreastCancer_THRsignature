
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
#################
THR_signature <- readxl::read_xlsx("./data/THR_Signatures_Jan25_2023.xlsx")

# get the THR50 signature
THR_70 <- THR_signature$`THR-70`[!is.na(THR_signature$`THR-70`)]

THR_70 <- gsub('-', '', THR_70)

# combine
CO130 <- c(ET60, THR_70)
CO130 <- CO130[!CO130 %in% 'Gene']

#############################################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')
##############################################

##############################################
# fix gene names
##############################################
# Fix in TCGA
setdiff(CO130, rownames(Expr_tcga_refAll))
grep('^FAM63A', rownames(Expr_tcga_refAll), value = TRUE) # MINDY1
grep('^FAM176A', rownames(Expr_tcga_refAll), value = TRUE) # EVA1A
grep('^LEPREL1', rownames(Expr_tcga_refAll), value = TRUE) # P3H2
grep('^DULLARD', rownames(Expr_tcga_refAll), value = TRUE) # SDHAF3

rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'DULLARD'] <- 'HSA011916'
rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'GPR56'] <- 'ADGRG1'
rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'FAM116B'] <- 'DENND6B'
rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'FAM46B'] <- 'TENT5B'

###############
# Fix in metabric
setdiff(CO130, rownames(Expr_metabric_refAll))
grep('^FAM63A', rownames(Expr_metabric_refAll), value = TRUE) # MINDY1
grep('^FAM176A', rownames(Expr_metabric_refAll), value = TRUE) # EVA1A
grep('^LEPREL1', rownames(Expr_metabric_refAll), value = TRUE) # P3H2
grep('^RSNL2', rownames(Expr_metabric_refAll), value = TRUE) # SDHAF3

rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'FAM63A'] <- 'MINDY1'
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'FAM176A'] <- 'EVA1A'
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'LEPREL1'] <- 'P3H2'
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'ACN9'] <- 'SDHAF3'
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'GPR56'] <- 'ADGRG1'
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'CTDNEP1'] <- 'HSA011916'
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'FAM116B'] <- 'DENND6B'
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'FAM46B'] <- 'TENT5B'

##############
# filter the signatures to include only the genes present in the expr matrices
CO130_fil <- CO130[CO130 %in% rownames(Expr_tcga_refAll) & CO130 %in% rownames(Expr_metabric_refAll)]

setdiff(CO130, CO130_fil)

##############################################
### combine in 1 dataset
##############################################
Data_metabric <- as.data.frame(cbind(t(Expr_metabric_refAll), group_metabric))
Data_metabric$group_metabric <- as.factor(Data_metabric$group_metabric)
levels(Data_metabric$group_metabric) <- c('0', '1')
colnames(Data_metabric)[colnames(Data_metabric) %in% c('group_metabric')] <- c('os')



#############################################################################################################
##############################################################################################################
# the model
CO130_model <- glm(as.formula((paste("os ~", paste(CO130_fil, collapse = "+")))), data = Data_metabric, family = "binomial")
summary(CO130_model)

save(CO130_model, file = 'objs/CO130_model_logreg.rda')

###########################################################################
############################################################################
### predict in the training dataset
# Make predictions

Train_prob_CO130 <- CO130_model %>% predict(Data_metabric , type = "response")

### Threshold
thr_CO130 <- coords(roc(group_metabric, Train_prob_CO130, levels = c("0", "1"), direction = "<"), "best")["threshold"]
thr_CO130

### ROC Curve
ROCTrain_CO130 <- roc(group_metabric, Train_prob_CO130, plot = F, print.thres=thr_CO130$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("0", "1"), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_CO130

### Get predictions based on best threshold from ROC curve
Train_predClasses_CO130 <- ifelse(Train_prob_CO130 >= thr_CO130$threshold, "1", "0")
table(Train_predClasses_CO130)
Train_predClasses_CO130 <- factor(Train_predClasses_CO130, levels = c('0', '1'))


### Resubstitution performance in the TRAINING set
ConfusionTrain_CO130 <- confusionMatrix(Train_predClasses_CO130, group_metabric, positive = "1", mode = "everything")
ConfusionTrain_CO130

## MCC
MCC_Train_CO130 <- mltools::mcc(pred = Train_predClasses_CO130, actuals = group_metabric)
MCC_Train_CO130

##########################
## Keep only the relevant information (Metastasis Event and Time)
Phenotype_metabric <- cbind(Pheno_metabric[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Relapse.Free.Status", "Relapse.Free.Status..Months.", "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC", "X3.Gene.classifier.subtype")], 
                            Train_prob_CO130, Train_predClasses_CO130)

# create a merged pdata and Z-scores object
CoxData_metabric <- data.frame(Phenotype_metabric)

# divide the probabilities into quartiles
CoxData_metabric <- CoxData_metabric %>%
  mutate(metabric_prob_CO130_quartiles = ntile(Train_prob_CO130, 4), 
         metabric_prob_CO130_quintiles = ntile(Train_prob_CO130, 5),
         metabric_prob_CO130_tertiles = ntile(Train_prob_CO130, 3)
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


########################
# OS
########################

result_prob_os <- ga(type = "real-valued", 
                  fitness = function(cutoffs) fitness_function(cutoffs, CoxData_metabric$Train_prob_CO130, 
                                                               CoxData_metabric$`Overall.Survival..Months.`,
                                                               CoxData_metabric$`Overall.Survival.Status`),
                  min = rep(min(CoxData_metabric$Train_prob_CO130), 2), 
                  max = rep(max(CoxData_metabric$Train_prob_CO130), 2), 
                  names = c("Cutoff 1", "Cutoff 2"),
                  popSize = 50, maxiter = 100, run = 100, seed = 123456)

# Extract the optimal cut-offs
optimal_cutoffs_prob_os <- result_prob_os@solution

optimal_cutoffs_prob_os

best_p_value_os <- Inf
best_cutoffs_os <- NULL

for (i in 1:nrow(optimal_cutoffs_prob_os)) {
  groups_prob_os <- findInterval(CoxData_metabric$Train_prob_CO130, c(-Inf, sort(optimal_cutoffs_prob_os[i, ]), Inf))
  surv_diff_result_prob_os <- survdiff(Surv(CoxData_metabric$`Overall.Survival..Months.`, CoxData_metabric$`Overall.Survival.Status`) ~ groups_prob_os)
  
  if (!is.null(surv_diff_result_prob_os$pval) && surv_diff_result_prob_os$pval < best_p_value_os) {
    best_p_value_os <- surv_diff_result_prob_os$pval
    best_cutoffs_os <- optimal_cutoffs_prob_os[i, ]
  }
}

# Create the survival object
os_surv_obj <- Surv(CoxData_metabric$`Overall.Survival..Months.`, CoxData_metabric$`Overall.Survival.Status`)

# Fit the survival curves
fit_os_bestcutoff_tertiles_prob <- survfit(os_surv_obj ~ groups_prob_os)
tiff("./figures/CO130_metabric_prognostic/CO130_metabric_os_bestCutoff_tertiles_prob.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(fit_os_bestcutoff_tertiles_prob, 
           data = CoxData_metabric, 
           xlim = c(0,240),
           break.x.by = 40,
           risk.table = TRUE, pval = TRUE,
           ggtheme = theme_survminer(base_size = 16, font.x = c(16, 'bold.italic', 'black'), font.y = c(16, 'bold.italic', 'black'), font.tickslab = c(16, 'plain', 'black'), font.legend = c(16, 'bold', 'black')),
           legend.labs = c('Q1', 'Q2', 'Q3')
           )
dev.off()

# Fit the Cox proportional hazards model using best cut-offs
cox_fit_os_bestcutoff_tertiles_prob <- coxph(os_surv_obj ~ as.factor(groups_prob_os), data = CoxData_metabric)
summary(cox_fit_os_bestcutoff_tertiles_prob)


########################
# RFS
########################

result_prob_rfs <- ga(type = "real-valued", 
                     fitness = function(cutoffs) fitness_function(cutoffs, CoxData_metabric$Train_prob_CO130, 
                                                                  CoxData_metabric$Relapse.Free.Status..Months.,
                                                                  CoxData_metabric$Relapse.Free.Status),
                     min = rep(min(CoxData_metabric$Train_prob_CO130), 2), 
                     max = rep(max(CoxData_metabric$Train_prob_CO130), 2), 
                     names = c("Cutoff 1", "Cutoff 2"),
                     popSize = 50, maxiter = 100, run = 100, seed = 123456)

# Extract the optimal cut-offs
optimal_cutoffs_prob_rfs <- result_prob_rfs@solution

optimal_cutoffs_prob_rfs

best_p_value_rfs <- Inf
best_cutoffs_rfs <- NULL

for (i in 1:nrow(optimal_cutoffs_prob_rfs)) {
  groups_prob_rfs <- findInterval(CoxData_metabric$Train_prob_CO130, c(-Inf, sort(optimal_cutoffs_prob_rfs[i, ]), Inf))
  surv_diff_result_prob_rfs <- survdiff(Surv(CoxData_metabric$Relapse.Free.Status..Months., CoxData_metabric$Relapse.Free.Status) ~ groups_prob_rfs)
  
  if (!is.null(surv_diff_result_prob_rfs$pval) && surv_diff_result_prob_rfs$pval < best_p_value_rfs) {
    best_p_value_rfs <- surv_diff_result_prob_rfs$pval
    best_cutoffs_rfs <- optimal_cutoffs_prob_rfs[i, ]
  }
}

# Create the survival object
rfs_surv_obj <- Surv(CoxData_metabric$Relapse.Free.Status..Months., CoxData_metabric$Relapse.Free.Status)

# Fit the survival curves
fit_rfs_bestcutoff_tertiles_prob <- survfit(rfs_surv_obj ~ groups_prob_rfs)
tiff("./figures/CO130_metabric_prognostic/CO130_metabric_rfs_bestCutoff_tertiles_prob.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(fit_rfs_bestcutoff_tertiles_prob, 
           data = CoxData_metabric, 
           xlim = c(0,240),
           break.x.by = 40,
           risk.table = TRUE, pval = TRUE,
           ggtheme = theme_survminer(base_size = 16, font.x = c(16, 'bold.italic', 'black'), font.y = c(16, 'bold.italic', 'black'), font.tickslab = c(16, 'plain', 'black'), font.legend = c(16, 'bold', 'black')),
           legend.labs = c('Q1', 'Q2', 'Q3')
)
dev.off()

# Fit the Cox proportional hazards model using best cut-offs
cox_fit_rfs_bestcutoff_tertiles_prob <- coxph(rfs_surv_obj ~ as.factor(groups_prob_rfs), data = CoxData_metabric)
summary(cox_fit_rfs_bestcutoff_tertiles_prob)




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
library(boot)
library(gridExtra)

#############################################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')
##############################################

# Fix gene names
# Replace special characters with an underscore
orig_gene_names <- rownames(Expr_metabric)
rownames(Expr_metabric) <- gsub("[^[:alnum:]_]", "_", rownames(Expr_metabric))

# Ensure names are valid R variable names
rownames(Expr_metabric) <- make.names(rownames(Expr_metabric), unique = TRUE)

# check if any names were changed
if (!all(orig_gene_names == rownames(Expr_metabric))) {
  warning("Some gene names were modified to ensure validity")
}

#############################
#############################
THR_6h <- c('CCDC148', 'MOAP1', 'RHOB', 'XBP1', 'AGGF1', 'TBCK')
THR_6b <- c('CDC20', 'PARP12', 'TRANK1', 'TAP2', 'LYN', 'IRF1')
THR_6e <- c('LMNB2', 'CDC20', 'KIF2C', 'FAM64A', 'KIF4A', 'TPX2')
THR_12t <- c('RARRES3', 'PARP12', 'IFIT3', 'B2M', 'TNFAIP8', 'GBP1', 'PSMB9', 'FYN', 'LYN', 'MSN', 'IL15RA', 'IL15')

##############################################
### combine in 1 dataset
##############################################
rownames(Expr_metabric) <- gsub('-', '', rownames(Expr_metabric))

os_event <- as.factor(Pheno_metabric$Overall.Survival.Status)
os_time <- Pheno_metabric$Overall.Survival..Months.
summary(os_time)

rfs_event <- as.factor(Pheno_metabric$Relapse.Free.Status)
rfs_time <- Pheno_metabric$Relapse.Free.Status..Months.
summary(rfs_time)

# make a dataframe
Data_metabric <- as.data.frame(t(Expr_metabric))

# add os info
Data_metabric$os_event <- os_event
Data_metabric$os_time <- os_time
summary(Data_metabric$os_time)

# add rfs info
Data_metabric$rfs_event <- rfs_event
Data_metabric$rfs_time <- rfs_time
summary(Data_metabric$rfs_time)

# add stage
table(Pheno_metabric$Tumor.Stage)
Data_metabric$stage <- as.factor(Pheno_metabric$Tumor.Stage)
Data_metabric$stage <- factor(Data_metabric$stage, levels= c('0', '1', '2', '3', '4'))
table(Data_metabric$stage)

# add grade
table(Pheno_metabric$Neoplasm.Histologic.Grade)
Data_metabric$grade <- as.factor(Pheno_metabric$Neoplasm.Histologic.Grade)
Data_metabric$grade <- factor(Data_metabric$grade, levels= c('1', '2', '3'))
table(Data_metabric$grade)

# add age
summary(Pheno_metabric$Age.at.Diagnosis)
Data_metabric$age <- as.numeric(Pheno_metabric$Age.at.Diagnosis)
summary(Data_metabric$age)

# add HER2 status
table(Pheno_metabric$HER2.Status)
Data_metabric$HER2.Status <- Pheno_metabric$HER2.Status
table(Data_metabric$HER2.Status)

# add ER status
table(Pheno_metabric$ER.status.measured.by.IHC)
Data_metabric$ER.Status <- Pheno_metabric$ER.status.measured.by.IHC
table(Data_metabric$ER.Status)

# add PAM50 subtypes
table(Pheno_metabric$Pam50...Claudin.low.subtype)
Data_metabric$PAM50 <- Pheno_metabric$Pam50...Claudin.low.subtype
table(Data_metabric$PAM50)

# check
all(rownames(Data_metabric) == rownames(Pheno_metabric))

# Ensure that 'os_even' and 'rfs_even' are numeric
Data_metabric$os_event <- as.numeric(as.character(Data_metabric$os_event))
Data_metabric$rfs_event <- as.numeric(as.character(Data_metabric$rfs_event))

####################################
# get HER2+
####################################
table(Data_metabric$HER2.Status)

Data_metabric_HER2_pos <- Data_metabric %>%
  dplyr::filter(HER2.Status == 'Positive')  

####################################
# get ER+
####################################
table(Data_metabric$ER.Status)

Data_metabric_ER_pos <- Data_metabric %>%
  dplyr::filter(ER.Status == 'Positve')  

####################################
# get ER-
####################################
table(Data_metabric$ER.Status)

Data_metabric_ER_neg <- Data_metabric %>%
  dplyr::filter(ER.Status == 'Negative')  

####################################
# get Luminal-B
####################################
table(Data_metabric$PAM50)

Data_metabric_LumB <- Data_metabric %>%
  dplyr::filter(PAM50 == 'LumB')  

#############################################################################################################
# Risk score and survival analysis
##############################################################################################################
generate_glm_risk_scores <- function(expr_data, signature_genes) {
  # Ensure signature genes are present in the expr_data columns
  if (!all(signature_genes %in% colnames(expr_data))) {
    stop("Not all signature genes are present in the expression dataset.")
  }
  
  if (!"os_event" %in% colnames(expr_data)) {
    stop("os_event variable not found in the data.")
  }
  
  formula_str <- paste("os_event ~", paste(signature_genes, collapse = " + "))
  
  # Fit the GLM with binomial family
  glm_fit <- glm(as.formula(formula_str), data = expr_data, family = binomial())
  
  # Check if glm fitting was successful
  if (!exists("glm_fit")) {
    stop("GLM fitting was unsuccessful.")
  }
  
  # Compute risk scores: Here we use the model's coefficients to compute a linear predictor
  risk_scores <- as.matrix(expr_data[, signature_genes, drop = FALSE]) %*% coef(glm_fit)[-1] + coef(glm_fit)[1]
  
  # Calculate the median of the risk scores
  median_score <- median(risk_scores)
  
  # Dichotomize patients based on the median risk score
  risk_categories <- ifelse(risk_scores > median_score, "High risk", "Low risk")
  risk_categories <- factor(risk_categories, levels = c("Low risk", "High risk"))
  # Returning risk categories
  return(risk_categories)
}

#########################################################################
# THR-6h in HER2-positive
######################################################################
THR6h_risk_categories <- generate_glm_risk_scores(Data_metabric_HER2_pos, THR_6h)

# Add risk scores to the dataframe
Data_metabric_HER2_pos$THR6h_risk_categories <- THR6h_risk_categories

# get the os survival estimates using km
THR6h_km_os_fit <- survfit(Surv(os_time, os_event) ~ THR6h_risk_categories, data = Data_metabric_HER2_pos)

# get the rfs survival estimates using km
THR6h_km_rfs_fit <- survfit(Surv(rfs_time, rfs_event) ~ THR6h_risk_categories, data = Data_metabric_HER2_pos)

# plot os
tiff("./figures/THR_G/THR_6h/THR6h_metabric_os.tiff", width = 2000, height = 2000, res = 350)
ggsurvplot(THR6h_km_os_fit,
           risk.table = TRUE,
           pval = TRUE,
           legend.labs = levels(Data_metabric_HER2_pos$THR6h_risk_categories),
           legend.title = '',
           xlim = c(0,240),
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
           risk.table.y.text.col = TRUE,
           palette = 'jco',
           risk.table.y.text = TRUE, 
           #title = 'THR 50 and METABRIC OS'
)
dev.off()

# plot rfs
tiff("./figures/THR_G/THR_6h/THR6h_metabric_rfs.tiff", width = 2000, height = 2000, res = 350)
ggsurvplot(THR6h_km_rfs_fit,
           risk.table = TRUE,
           pval = TRUE,
           legend.labs = levels(Data_metabric_HER2_pos$THR6h_risk_categories),
           legend.title = '',
           xlim = c(0,240),
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
           risk.table.y.text.col = TRUE,
           palette = 'jco',
           risk.table.y.text = TRUE, 
           #title = 'THR 50 and METABRIC OS'
)
dev.off()

# multivariate COX: os 
THR6h_cox_multi_os <- coxph(Surv(os_time, os_event) ~ THR6h_risk_categories + age + stage + grade, data = Data_metabric_HER2_pos)
summary(THR6h_cox_multi_os)

# multivariate COX: rfs 
THR6h_cox_multi_rfs <- coxph(Surv(rfs_time, rfs_event) ~ THR6h_risk_categories + age + stage + grade, data = Data_metabric_HER2_pos)
summary(THR6h_cox_multi_rfs)

#########################################################################
# THR-6e in ER-positive
######################################################################
THR6e_risk_categories <- generate_glm_risk_scores(Data_metabric_ER_pos, THR_6e)

# Add risk scores to the dataframe
Data_metabric_ER_pos$THR6e_risk_categories <- THR6e_risk_categories

# get the survival estimates using km
THR6e_km_os_fit <- survfit(Surv(os_time, os_event) ~ THR6e_risk_categories, data = Data_metabric_ER_pos)

# get the rfs survival estimates using km
THR6e_km_rfs_fit <- survfit(Surv(rfs_time, rfs_event) ~ THR6e_risk_categories, data = Data_metabric_ER_pos)

# plot os
tiff("./figures/THR_G/THR_6e/THR6e_metabric_os.tiff", width = 2000, height = 2000, res = 350)
ggsurvplot(THR6e_km_os_fit,
           risk.table = TRUE,
           pval = TRUE,
           legend.labs = levels(Data_metabric_ER_pos$THR6e_risk_categories),
           legend.title = '',
           xlim = c(0,240),
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
           risk.table.y.text.col = TRUE,
           palette = 'jco',
           risk.table.y.text = TRUE, 
           #title = 'THR 50 and METABRIC OS'
)
dev.off()

# plot rfs
tiff("./figures/THR_G/THR_6e/THR6e_metabric_rfs.tiff", width = 2000, height = 2000, res = 350)
ggsurvplot(THR6e_km_rfs_fit,
           risk.table = TRUE,
           pval = TRUE,
           legend.labs = levels(Data_metabric_ER_pos$THR6e_risk_categories),
           legend.title = '',
           xlim = c(0,240),
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
           risk.table.y.text.col = TRUE,
           palette = 'jco',
           risk.table.y.text = TRUE, 
           #title = 'THR 50 and METABRIC OS'
)
dev.off()

# multivariate COX: os
THR6e_cox_multi_os <- coxph(Surv(os_time, os_event) ~ THR6e_risk_categories + age + stage + grade, data = Data_metabric_ER_pos)
summary(THR6e_cox_multi_os)

# multivariate COX: rfs
THR6e_cox_multi_rfs <- coxph(Surv(rfs_time, rfs_event) ~ THR6e_risk_categories + age + stage + grade, data = Data_metabric_ER_pos)
summary(THR6e_cox_multi_rfs)


#########################################################################
# THR-6b in Luminal B
######################################################################
THR6b_risk_categories <- generate_glm_risk_scores(Data_metabric_LumB, THR_6b)

# Add risk scores to the dataframe
Data_metabric_LumB$THR6b_risk_categories <- THR6b_risk_categories

# get the survival estimates using km
THR6b_km_os_fit <- survfit(Surv(os_time, os_event) ~ THR6b_risk_categories, data = Data_metabric_LumB)

# get the rfs survival estimates using km
THR6b_km_rfs_fit <- survfit(Surv(rfs_time, rfs_event) ~ THR6b_risk_categories, data = Data_metabric_LumB)

# plot os
tiff("./figures/THR_G/THR_6b/THR6b_metabric_os.tiff", width = 2000, height = 2000, res = 350)
ggsurvplot(THR6b_km_os_fit,
           risk.table = TRUE,
           pval = TRUE,
           legend.labs = levels(Data_metabric_LumB$THR6b_risk_categories),
           legend.title = '',
           xlim = c(0,240),
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
           risk.table.y.text.col = TRUE,
           palette = 'jco',
           risk.table.y.text = TRUE, 
           #title = 'THR 50 and METABRIC OS'
)
dev.off()

# plot rfs
tiff("./figures/THR_G/THR_6b/THR6b_metabric_rfs.tiff", width = 2000, height = 2000, res = 350)
ggsurvplot(THR6b_km_rfs_fit,
           risk.table = TRUE,
           pval = TRUE,
           legend.labs = levels(Data_metabric_LumB$THR6b_risk_categories),
           legend.title = '',
           xlim = c(0,240),
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
           risk.table.y.text.col = TRUE,
           palette = 'jco',
           risk.table.y.text = TRUE, 
           #title = 'THR 50 and METABRIC OS'
)
dev.off()

# multivariate COX: os
THR6b_cox_multi_os <- coxph(Surv(os_time, os_event) ~ THR6b_risk_categories + age + stage + grade, data = Data_metabric_LumB)
summary(THR6b_cox_multi_os)

# multivariate COX: rfs
THR6b_cox_multi_rfs <- coxph(Surv(rfs_time, rfs_event) ~ THR6b_risk_categories + age + stage + grade, data = Data_metabric_LumB)
summary(THR6b_cox_multi_rfs)

#########################################################################
# THR-12t in ER-negative
######################################################################
THR12t_risk_categories <- generate_glm_risk_scores(Data_metabric_ER_neg, THR_12t)

# Add risk scores to the dataframe
Data_metabric_ER_neg$THR12t_risk_categories <- THR12t_risk_categories

# get the survival estimates using km
THR12t_km_os_fit <- survfit(Surv(os_time, os_event) ~ THR12t_risk_categories, data = Data_metabric_ER_neg)

# get the rfs survival estimates using km
THR12t_km_rfs_fit <- survfit(Surv(rfs_time, rfs_event) ~ THR12t_risk_categories, data = Data_metabric_ER_neg)

# plot os
tiff("./figures/THR_G/THR_12t/THR12t_metabric_os.tiff", width = 2000, height = 2000, res = 350)
ggsurvplot(THR12t_km_os_fit,
           risk.table = TRUE,
           pval = TRUE,
           legend.labs = levels(Data_metabric_ER_neg$THR12t_risk_categories),
           legend.title = '',
           xlim = c(0,240),
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
           risk.table.y.text.col = TRUE,
           palette = 'jco',
           risk.table.y.text = TRUE, 
           #title = 'THR 50 and METABRIC OS'
)
dev.off()

# plot rfs
tiff("./figures/THR_G/THR_12t/THR12t_metabric_rfs.tiff", width = 2000, height = 2000, res = 350)
ggsurvplot(THR12t_km_rfs_fit,
           risk.table = TRUE,
           pval = TRUE,
           legend.labs = levels(Data_metabric_ER_neg$THR12t_risk_categories),
           legend.title = '',
           xlim = c(0,240),
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
           risk.table.y.text.col = TRUE,
           palette = 'jco',
           risk.table.y.text = TRUE, 
           #title = 'THR 50 and METABRIC OS'
)
dev.off()

# multivariate COX: os
THR12t_cox_multi_os <- coxph(Surv(os_time, os_event) ~ THR12t_risk_categories + age + stage + grade, data = Data_metabric_ER_neg)
summary(THR12t_cox_multi_os)

# multivariate COX: rfs
THR12t_cox_multi_rfs <- coxph(Surv(rfs_time, rfs_event) ~ THR12t_risk_categories + age + stage + grade, data = Data_metabric_ER_neg)
summary(THR12t_cox_multi_rfs)

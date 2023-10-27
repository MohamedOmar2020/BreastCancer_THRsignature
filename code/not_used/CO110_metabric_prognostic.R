
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
load('./objs/forKTSP.rda')
##############################################

##############################################
# fix gene names
##############################################

# Fix in TCGA
setdiff(CO110, rownames(Expr_tcga_refAll))
grep('^FAM63A', rownames(Expr_tcga_refAll), value = TRUE) # MINDY1
grep('^FAM176A', rownames(Expr_tcga_refAll), value = TRUE) # EVA1A
grep('^LEPREL1', rownames(Expr_tcga_refAll), value = TRUE) # P3H2
grep('^ACN9', rownames(Expr_tcga_refAll), value = TRUE) # SDHAF3

rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'FAM63A'] <- 'MINDY1'
rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'FAM176A'] <- 'EVA1A'
rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'LEPREL1'] <- 'P3H2'
rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'ACN9'] <- 'SDHAF3'
rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'CTDNEP1'] <- 'HSA011916'
rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'GPR56'] <- 'ADGRG1'
rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'FAM116B'] <- 'DENND6B'
rownames(Expr_tcga_refAll)[rownames(Expr_tcga_refAll) == 'FAM46B'] <- 'TENT5B'

###############
# Fix in metabric
setdiff(CO110, rownames(Expr_metabric_refAll))
grep('^FAM63A', rownames(Expr_metabric_refAll), value = TRUE) # MINDY1
grep('^FAM176A', rownames(Expr_metabric_refAll), value = TRUE) # EVA1A
grep('^LEPREL1', rownames(Expr_metabric_refAll), value = TRUE) # P3H2
grep('^ACN9', rownames(Expr_metabric_refAll), value = TRUE) # SDHAF3

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
CO110_fil <- CO110[CO110 %in% rownames(Expr_tcga_refAll) & CO110 %in% rownames(Expr_metabric_refAll)]

setdiff(CO110, CO110_fil)

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
CO110_model <- glm(as.formula((paste("os ~", paste(CO110_fil, collapse = "+")))), data = Data_metabric, family = "binomial")
summary(CO110_model)

save(CO110_model, file = 'objs/CO110_model_logreg.rda')
###########################################################################
############################################################################
### predict in the training dataset
# Make predictions

Train_prob_CO110 <- CO110_model %>% predict(Data_metabric , type = "response")

### Threshold
thr_CO110 <- coords(roc(group_metabric, Train_prob_CO110, levels = c("0", "1"), direction = "<"), "best")["threshold"]
thr_CO110

### ROC Curve
ROCTrain_CO110 <- roc(group_metabric, Train_prob_CO110, plot = F, print.thres=thr_CO110$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("0", "1"), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_CO110

### Get predictions based on best threshold from ROC curve
Train_predClasses_CO110 <- ifelse(Train_prob_CO110 >= thr_CO110$threshold, "1", "0")
table(Train_predClasses_CO110)
Train_predClasses_CO110 <- factor(Train_predClasses_CO110, levels = c('0', '1'))


### Resubstitution performance in the TRAINING set
ConfusionTrain_CO110 <- confusionMatrix(Train_predClasses_CO110, group_metabric, positive = "1", mode = "everything")
ConfusionTrain_CO110

## MCC
MCC_Train_CO110 <- mltools::mcc(pred = Train_predClasses_CO110, actuals = group_metabric)
MCC_Train_CO110

##########################
## Keep only the relevant information (Metastasis Event and Time)
Phenotype_metabric <- cbind(Pheno_metabric[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Relapse.Free.Status", "Relapse.Free.Status..Months.", "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC", "X3.Gene.classifier.subtype")], 
                            Train_prob_CO110, Train_predClasses_CO110)

# create a merged pdata and Z-scores object
CoxData_metabric <- data.frame(Phenotype_metabric)

# divide the probabilities into quartiles
CoxData_metabric <- CoxData_metabric %>%
  mutate(metabric_prob_CO110_quartiles = ntile(Train_prob_CO110, 4), 
         metabric_prob_CO110_quintiles = ntile(Train_prob_CO110, 5),
         metabric_prob_CO110_tertiles = ntile(Train_prob_CO110, 3)
         )


########################################################################  
## Fit survival curves
########################################################################  

# OS
## metabric all genes
Fit_sig_metabric_os_CO110 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_CO110, data = CoxData_metabric)


## by quartiles
Fit_sig_metabric_os_CO110_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_quartiles, data = CoxData_metabric)

## by quintiles
Fit_sig_metabric_os_CO110_quintiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_quintiles, data = CoxData_metabric)

## by tertiles
Fit_sig_metabric_os_CO110_tertiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_tertiles, data = CoxData_metabric)

#############
# RFS
## metabric all genes
Fit_sig_metabric_RFS_CO110 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_predClasses_CO110, data = CoxData_metabric)

## by quartiles
Fit_sig_metabric_RFS_CO110_quartiles <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_quartiles, data = CoxData_metabric)

## by quintiles
Fit_sig_metabric_RFS_CO110_quintiles <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_quintiles, data = CoxData_metabric)

## by tertiles
Fit_sig_metabric_RFS_CO110_tertiles <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_tertiles, data = CoxData_metabric)

#################################
## by clinical groups

## pam50

# keep only the major pam50 subtypes
CoxData_metabric_PAM <- CoxData_metabric %>%
  filter(Pam50...Claudin.low.subtype %in% c('Basal', 'claudin-low', 'Her2', 'LumA', 'LumB'))

# keep only quartiles 1 and 4 (has to be in each of THR25, CO110_1, and CO110_2)

CoxData_metabric_PAM_Q1vsQ4_CO110 <- CoxData_metabric_PAM %>%
  filter(metabric_prob_CO110_quartiles %in% c('1', '4'))

CoxData_metabric_PAM_Q1vsQ5_CO110 <- CoxData_metabric_PAM %>%
  filter(metabric_prob_CO110_quintiles %in% c('1', '5'))


# os: all
Fit_sig_metabric_os_CO110_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_CO110 + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM)

# os: quartiles: all
Fit_sig_metabric_os_CO110_quartiles_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_quartiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM)

# os: quintiles: all
Fit_sig_metabric_os_CO110_quintiles_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_quintiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM)

# os: quartiles: Q1 vs Q4
Fit_sig_metabric_os_CO110_Q1vsQ4_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_quartiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM_Q1vsQ4_CO110)

# os: quintiles: Q1 vs Q5
Fit_sig_metabric_os_CO110_Q1vsQ5_PAM <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_quintiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM_Q1vsQ5_CO110)

##########
# rfs: all
Fit_sig_metabric_rfs_CO110_PAM <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_predClasses_CO110 + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM)

# rfs: quartiles: all
Fit_sig_metabric_rfs_CO110_quartiles_PAM <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_quartiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM)

# rfs: quintiles: all
Fit_sig_metabric_rfs_CO110_quintiles_PAM <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_quintiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM)

# rfs: quartiles: Q1 vs Q4
Fit_sig_metabric_rfs_CO110_Q1vsQ4_PAM <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_quartiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM_Q1vsQ4_CO110)

# rfs: quintiles: Q1 vs Q5
Fit_sig_metabric_rfs_CO110_Q1vsQ5_PAM <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_quintiles + Pam50...Claudin.low.subtype, data = CoxData_metabric_PAM_Q1vsQ5_CO110)

#############
## ER

# keep only quartiles 1 and 4 (has to be in each of THR25, CO110_1, and CO110_2)

# Note: we can use this for both ER and X3
CoxData_metabric_Q1vsQ4_CO110 <- CoxData_metabric %>%
  filter(metabric_prob_CO110_quartiles %in% c('1', '4'))

CoxData_metabric_Q1vsQ5_CO110 <- CoxData_metabric %>%
  filter(metabric_prob_CO110_quintiles %in% c('1', '5'))

# os: all
Fit_sig_metabric_os_CO110_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_CO110 + ER.status.measured.by.IHC, data = CoxData_metabric)

# os: quartiles: all
Fit_sig_metabric_os_CO110_quartiles_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_quartiles + ER.status.measured.by.IHC, data = CoxData_metabric)

# os: quintiles: all
Fit_sig_metabric_os_CO110_quintiles_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_quintiles + ER.status.measured.by.IHC, data = CoxData_metabric)

# os: quartiles: Q1 vs Q4
Fit_sig_metabric_os_CO110_Q1vsQ4_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_quartiles + ER.status.measured.by.IHC, data = CoxData_metabric_Q1vsQ4_CO110)

# os: quintiles: Q1 vs Q5
Fit_sig_metabric_os_CO110_Q1vsQ5_ER <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_quintiles + ER.status.measured.by.IHC, data = CoxData_metabric_Q1vsQ5_CO110)

##############
# rfs: all
Fit_sig_metabric_rfs_CO110_ER <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_predClasses_CO110 + ER.status.measured.by.IHC, data = CoxData_metabric)

# rfs: quartiles: all
Fit_sig_metabric_rfs_CO110_quartiles_ER <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_quartiles + ER.status.measured.by.IHC, data = CoxData_metabric)

# rfs: quintiles: all
Fit_sig_metabric_rfs_CO110_quintiles_ER <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_quintiles + ER.status.measured.by.IHC, data = CoxData_metabric)

# rfs: quartiles: Q1 vs Q4
Fit_sig_metabric_rfs_CO110_Q1vsQ4_ER <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_quartiles + ER.status.measured.by.IHC, data = CoxData_metabric_Q1vsQ4_CO110)

# rfs: quintiles: Q1 vs Q5
Fit_sig_metabric_rfs_CO110_Q1vsQ5_ER <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_quintiles + ER.status.measured.by.IHC, data = CoxData_metabric_Q1vsQ5_CO110)

#######################################
## X3

# os: all
Fit_sig_metabric_os_CO110_X3 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_CO110 + X3.Gene.classifier.subtype, data = CoxData_metabric)

# os: quartiles: all
Fit_sig_metabric_os_CO110_quartiles_X3 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_quartiles + X3.Gene.classifier.subtype, data = CoxData_metabric)

# os: quintiles: all
Fit_sig_metabric_os_CO110_quintiles_X3 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_quintiles + X3.Gene.classifier.subtype, data = CoxData_metabric)

# os: quartiles: Q1 vs Q4
Fit_sig_metabric_os_CO110_Q1vsQ4_X3 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_quartiles + X3.Gene.classifier.subtype, data = CoxData_metabric_Q1vsQ4_CO110)

# os: quintiles: Q1 vs Q5
Fit_sig_metabric_os_CO110_Q1vsQ5_X3 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_quintiles + X3.Gene.classifier.subtype, data = CoxData_metabric_Q1vsQ5_CO110)

###########
# rfs: all
Fit_sig_metabric_rfs_CO110_X3 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_predClasses_CO110 + X3.Gene.classifier.subtype, data = CoxData_metabric)

# rfs: quartiles: all
Fit_sig_metabric_rfs_CO110_quartiles_X3 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_quartiles + X3.Gene.classifier.subtype, data = CoxData_metabric)

# rfs: quintiles: all
Fit_sig_metabric_rfs_CO110_quintiles_X3 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_quintiles + X3.Gene.classifier.subtype, data = CoxData_metabric)

# rfs: quartiles: Q1 vs Q4
Fit_sig_metabric_rfs_CO110_Q1vsQ4_X3 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_quartiles + X3.Gene.classifier.subtype, data = CoxData_metabric_Q1vsQ4_CO110)

# rfs: quintiles: Q1 vs Q5
Fit_sig_metabric_rfs_CO110_Q1vsQ5_X3 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_quintiles + X3.Gene.classifier.subtype, data = CoxData_metabric_Q1vsQ5_CO110)

############################################################################
# plot OS
############################################################################

tiff("./figures/CO110_metabric/logreg/CO110_metabric_os_all.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_CO110,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           #title = 'THR 50 and METABRIC OS'
)
dev.off()

########
# by quartiles
tiff("./figures/CO110_metabric/logreg/CO110_metabric_os_quartiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_CO110_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles'
)
dev.off()

#############
# by quintiles
tiff("./figures/CO110_metabric/logreg/CO110_metabric_os_quintiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_CO110_quintiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles'
)
dev.off()

#############
# by tertiles
tiff("./figures/CO110_metabric/logreg/CO110_metabric_os_tertiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_CO110_tertiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('Q1', 'Q2', 'Q3'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles'
)
dev.off()

######################################
# plot RFS
######################################
tiff("./figures/CO110_metabric/logreg/CO110_metabric_RFS_allpairs.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_CO110,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           #title = 'THR 50_1 (logistic regression) and METABRIC RFS'
)
dev.off()


########
# by quartiles

pdf("./figures/CO110_metabric/logreg/CO110_metabric_RFS_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_RFS_CO110_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           #risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and METABRIC RFS: quartiles'
           )
dev.off()

########
# by quintiles

tiff("./figures/CO110_metabric/logreg/CO110_metabric_rfs_quintiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_CO110_quintiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles'
)
dev.off()


########
# by tertiles

tiff("./figures/CO110_metabric/logreg/CO110_metabric_rfs_tertiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_CO110_tertiles,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           legend.labs = c('Q1', 'Q2', 'Q3'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles'
)
dev.off()

############################################################################
############################################################################
### by clinical group

## PAM50

# OS

pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_os_PAM50.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_os_CO110_PAM,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           ggtheme = theme_minimal(),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'CO110 and METABRIC OS by PAM50 subtypes')
dev.off()


#########
# RFS

pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_RFS_PAM50.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_rfs_CO110_PAM,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           ggtheme = theme_minimal(),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'CO110 and METABRIC RFS by PAM50 subtypes')
dev.off()


######################################################
# OS: quartiles: all


pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_os_PAM50_quartiles.pdf", width = 12, height = 10, onefile = F)
ggsurvplot(Fit_sig_metabric_os_CO110_quartiles_PAM,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           ggtheme = theme_survminer(base_size = 25, font.x = c(25, 'bold.italic', 'black'), font.y = c(25, 'bold.italic', 'black'), font.tickslab = c(25, 'plain', 'black'), font.legend = c(25, 'bold', 'black')),
           palette = 'jco',
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'CO110 and METABRIC OS by PAM50 subtypes: quartiles')
dev.off()

####################################################
## OS: quintiles: all

tiff("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_os_PAM50_quintiles.tiff", width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_os_CO110_quintiles_PAM,
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
           #title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles'
)
dev.off()

######################################################
# OS: quartiles: Q1 vs Q4

tiff("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_os_PAM50_Q1vsQ4.tiff",  width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_os_CO110_Q1vsQ4_PAM,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           legend.title	= 'Quartiles',
           pval.size = 15,
           #break.x.by = 20,
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           legend.labs = c('Q1', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC OS by PAM50 subtypes: Q1 vs Q4'
)
dev.off()

######################################################
# OS: quintiles: Q1 vs Q5
tiff("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_os_PAM50_Q1vsQ5.tiff",  width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_os_CO110_Q1vsQ5_PAM,
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
           #title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles'
)
dev.off()

#############################
# RFS: quartiles: all

tiff("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_rfs_PAM50_quartiles.tiff",  width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_rfs_CO110_quartiles_PAM,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           legend.title	= 'Quartiles',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text = FALSE, title = 'CO110 and METABRIC RFS by PAM50 subtypes: quartiles')
dev.off()

####################################################
## RFS: quintiles: all

tiff("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_rfs_PAM50_quintiles.tiff",  width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_rfs_CO110_quintiles_PAM,
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
           #title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles'
)
dev.off()

#############################
# RFS: quartiles: Q1 vs Q4
tiff("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_rfs_PAM50_Q1vsQ4.tiff",  width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_rfs_CO110_Q1vsQ4_PAM,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 12,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           #legend.labs = c('Basal', 'Claudin-low', 'Her2+', 'Luminal A', 'Luminal B'),
           legend.title	= 'Quartiles',
           #break.x.by = 20,
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q4'),
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC RFS by PAM50 subtypes: Q1 vs Q4'
)
dev.off()

#############################
# RFS: quintiles: Q1 vs Q5
tiff("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_rfs_PAM50_Q1vsQ5.tiff",  width = 3200, height = 2200, res = 300)
ggsurvplot(Fit_sig_metabric_rfs_CO110_Q1vsQ5_PAM,
           risk.table = FALSE,
           pval = T,
           short.panel.labs = T,
           facet.by = "Pam50...Claudin.low.subtype",
           #legend.labs = c('Basal', 'Claudin-low', 'Her2+', 'Luminal A', 'Luminal B'),
           legend.title	= 'Quintiles',
           pval.size = 12,
           #break.x.by = 20,
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q5'),
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC RFS by PAM50 subtypes: Q1 vs Q4'
)
dev.off()
##############################################################################################
##############################################################################################
## ER

# OS

tiff("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_os_ER.tiff", width = 3000, height = 2000, res =300)
ggsurvplot(Fit_sig_metabric_os_CO110_ER,
           risk.table = FALSE,
           pval = FALSE,
           short.panel.labs = T,
           facet.by = "ER.status.measured.by.IHC",
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC OS by ER status'
           )
dev.off()

# Os COXPH by class
lapply(split(CoxData_metabric, CoxData_metabric$ER.status.measured.by.IHC),
       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_CO110, data = x)))

# by prob
lapply(split(CoxData_metabric, CoxData_metabric$ER.status.measured.by.IHC),
       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_prob_CO110, data = x)))

######################
# RFS

tiff("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_RFS_ER.tiff", width = 3000, height = 2000, res =300)
ggsurvplot(Fit_sig_metabric_rfs_CO110_ER,
           risk.table = FALSE,
           pval = FALSE,
           short.panel.labs = T,
           facet.by = "ER.status.measured.by.IHC",
           palette = 'jco',
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(17, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC OS by ER status'
)
dev.off()


# RFS COXPH by class
lapply(split(CoxData_metabric, CoxData_metabric$ER.status.measured.by.IHC),
       function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_predClasses_CO110, data = x)))

# by prob
lapply(split(CoxData_metabric, CoxData_metabric$ER.status.measured.by.IHC),
       function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_prob_CO110, data = x)))

####################################################
# OS: quartiles: all

pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_os_ER_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_os_CO110_quartiles_ER,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "ER.status.measured.by.IHC",
           ggtheme = theme_minimal(),
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC OS by ER status: quartiles'
           )
dev.off()



####################################################
# OS: quartiles: Q1 vs Q4
pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_os_ER_Q1vsQ4.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_os_CO110_Q1vsQ4_ER,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "ER.status.measured.by.IHC",
           ggtheme = theme_minimal(),
           legend.labs = c('Q1', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = THR 50_1 and METABRIC OS by ER status: Q1 vs Q4'
           )
dev.off()

####################################################
# RFS: quartiles: all
pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_rfs_ER_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_rfs_CO110_quartiles_ER,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "ER.status.measured.by.IHC",
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC RFS by ER status: quartiles'
           )
dev.off()

####################################################
# RFS: quartiles: Q1 vs Q4
pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_rfs_ER_Q1vsQ4.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_rfs_CO110_Q1vsQ4_ER,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "ER.status.measured.by.IHC",
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           legend.labs = c('Q1', 'Q4'),
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC RFS by ER status: Q1 vs Q4'
           )
dev.off()

##############################################################################################
##############################################################################################
## X3

# OS
pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_os_X3.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_os_CO110_X3,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC OS by X3 classifier subtypes'
           )
dev.off()


######################
# RFS
pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_rfs_X3.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_rfs_CO110_X3,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC RFS by X3 classifier subtypes'
           )
dev.off()


####################################################
# OS: quartiles: all

pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_os_X3_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_os_CO110_quartiles_X3,
           risk.table = FALSE,
           pval = TRUE,
           short.panel.labs = T,
           facet.by = "X3.Gene.classifier.subtype",
           ggtheme = theme_minimal(),
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC OS by X3 classifier subtypes: quartiles'
           )
dev.off()

####################################################
## OS: quintiles: all
pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_os_X3_quintiles.tiff", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_os_CO110_quintiles_X3,
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
           #title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles'
)
dev.off()

####################################################
# OS: quartiles: Q1 vs Q4
pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_os_X3_Q1vsQ4.tiff", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_os_CO110_Q1vsQ4_X3,
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
           #title = 'THR 50_1 and METABRIC OS by PAM50 subtypes: Q1 vs Q4'
)
dev.off()

####################################################
## OS: quintiles: Q1 vs Q5
pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_os_X3_Q1vsQ5.tiff", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_os_CO110_Q1vsQ5_X3,
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
           #title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles'
)
dev.off()


###############################################################
# RFS: quartiles: all
pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_rfs_X3_quartiles.tiff", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_rfs_CO110_quartiles_X3,
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
           #title = 'THR 50_1 and METABRIC RFS by X3 classifier subtypes: quartiles'
)
dev.off()

####################################################
## RFS: quintiles: all
pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_rfs_X3_quintiles.tiff", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_rfs_CO110_quintiles_X3,
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
           #title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles'
)
dev.off()


###############################################################
# RFS: quartiles: Q1 vs Q4

pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_rfs_X3_Q1vsQ4.tiff", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_rfs_CO110_Q1vsQ4_X3,
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
           #title = 'THR 50_1 and METABRIC RFS by X3 classifier subtypes: Q1 vs Q4'
)
dev.off()


####################################################
## RFS: quintiles: Q1 vs Q5
pdf("./figures/CO110_metabric/logreg/byClinicalGroup/CO110_metabric_rfs_X3_Q1vsQ5.tiff", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_rfs_CO110_Q1vsQ5_X3,
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
           #title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles'
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

Fit_sig_metabric_os_coxph_CO110_1 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_prob_CO110_1, data = CoxData_metabric)
summary(Fit_sig_metabric_os_coxph_CO110_1)

#Fit_sig_metabric_os_coxph_CO110_2 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_prob_CO110_2, data = CoxData_metabric)
#summary(Fit_sig_metabric_os_coxph_CO110_2)

#png('./figures/logreg/THR25_HR_metabric_os.png', width = 2000, height = 2000, res = 300)
#ggforest(Fit_sig_metabric_coxph_THR25, fontsize = 0.5)
#dev.off()

#png('./figures/logreg/CO110_1_HR_metabric_os.png', width = 2000, height = 2000, res = 300)
#ggforest(Fit_sig_metabric_coxph_CO110_1, fontsize = 0.5)
#dev.off()

#png('./figures/logreg/CO110_2_HR_metabric_os.png', width = 2000, height = 2000, res = 300)
#ggforest(Fit_sig_metabric_coxph_CO110_2, fontsize = 0.5)
#dev.off()

########
## by quartiles

# make a factor with Q1 (lowest risk) being the reference
#CoxData_metabric$metabric_prob_THR25_quartiles <- factor(CoxData_metabric$metabric_prob_THR25_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_metabric$metabric_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_metabric$metabric_prob_THR25_quartiles))

CoxData_metabric$metabric_prob_CO110_1_quartiles <- factor(CoxData_metabric$metabric_prob_CO110_1_quartiles, levels = c('1', '2', '3', '4'))
levels(CoxData_metabric$metabric_prob_CO110_1_quartiles) <- paste0('Q', levels(CoxData_metabric$metabric_prob_CO110_1_quartiles))

#CoxData_metabric$metabric_prob_CO110_2_quartiles <- factor(CoxData_metabric$metabric_prob_CO110_2_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_metabric$metabric_prob_CO110_2_quartiles) <- paste0('Q', levels(CoxData_metabric$metabric_prob_CO110_2_quartiles))

# fit
#Fit_sig_metabric_os_coxph_THR25_quartiles <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR25_quartiles, data = CoxData_metabric)
#summary_metabric_os_coxph_THR25_quartiles <- summary(Fit_sig_metabric_os_coxph_THR25_quartiles)

Fit_sig_metabric_os_coxph_CO110_1_quartiles <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_1_quartiles, data = CoxData_metabric)
summary_metabric_os_coxph_CO110_1_quartiles <- summary(Fit_sig_metabric_os_coxph_CO110_1_quartiles)

#Fit_sig_metabric_os_coxph_CO110_2_quartiles <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_2_quartiles, data = CoxData_metabric)
#summary_metabric_os_coxph_CO110_2_quartiles <- summary(Fit_sig_metabric_os_coxph_CO110_2_quartiles)

#summary_list_metabric_os_quartiles <- list(THR25 = summary_metabric_os_coxph_THR25_quartiles, CO110_1 = summary_metabric_os_coxph_CO110_1_quartiles, CO110_2 = summary_metabric_os_coxph_CO110_2_quartiles) 

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
# HR_df_metabric_os_coxph_quartiles_CO110_1 <- as.data.frame(HR_list_metabric_os_coxph_quartiles$CO110_1)
# HR_df_metabric_os_coxph_quartiles_CO110_1$quartile <- gsub('.+quartiles', '', rownames(HR_df_metabric_os_coxph_quartiles_CO110_1))
# 
# HR_df_metabric_os_coxph_quartiles_CO110_2 <- as.data.frame(HR_list_metabric_os_coxph_quartiles$CO110_2)
# HR_df_metabric_os_coxph_quartiles_CO110_2$quartile <- gsub('.+quartiles', '', rownames(HR_df_metabric_os_coxph_quartiles_CO110_2))
# 
# # save the results
# write.csv(HR_df_metabric_os_coxph_quartiles_THR25, 'objs/HR/metabric/OS/THR25_quartiles_HR.csv')
# write.csv(HR_df_metabric_os_coxph_quartiles_CO110_1, 'objs/HR/metabric/OS/CO110_1_quartiles_HR.csv')
# write.csv(HR_df_metabric_os_coxph_quartiles_CO110_2, 'objs/HR/metabric/OS/CO110_2_quartiles_HR.csv')

################
# RFS

# by probaility

#Fit_sig_metabric_RFS_coxph_THR25 <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_prob_THR25, data = CoxData_metabric)
#summary(Fit_sig_metabric_RFS_coxph_THR25)

Fit_sig_metabric_RFS_coxph_CO110_1 <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_prob_CO110_1, data = CoxData_metabric)
summary(Fit_sig_metabric_RFS_coxph_CO110_1)

#Fit_sig_metabric_RFS_coxph_CO110_2 <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_prob_CO110_2, data = CoxData_metabric)
#summary(Fit_sig_metabric_RFS_coxph_CO110_2)


########
## by quartiles

# fit
#Fit_sig_metabric_RFS_coxph_THR25_quartiles <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR25_quartiles, data = CoxData_metabric)
#summary_metabric_RFS_coxph_THR25_quartiles <- summary(Fit_sig_metabric_RFS_coxph_THR25_quartiles)

Fit_sig_metabric_RFS_coxph_CO110_1_quartiles <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_1_quartiles, data = CoxData_metabric)
summary_metabric_RFS_coxph_CO110_1_quartiles <- summary(Fit_sig_metabric_RFS_coxph_CO110_1_quartiles)

#Fit_sig_metabric_RFS_coxph_CO110_2_quartiles <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_2_quartiles, data = CoxData_metabric)
#summary_metabric_RFS_coxph_CO110_2_quartiles <- summary(Fit_sig_metabric_RFS_coxph_CO110_2_quartiles)

#summary_list_metabric_RFS_quartiles <- list(THR25 = summary_metabric_RFS_coxph_THR25_quartiles, CO110_1 = summary_metabric_RFS_coxph_CO110_1_quartiles, CO110_2 = summary_metabric_RFS_coxph_CO110_2_quartiles) 

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
# HR_df_metabric_RFS_coxph_quartiles_CO110_1 <- as.data.frame(HR_list_metabric_RFS_coxph_quartiles$CO110_1)
# HR_df_metabric_RFS_coxph_quartiles_CO110_1$quartile <- gsub('.+quartiles', '', rownames(HR_df_metabric_RFS_coxph_quartiles_CO110_1))
# 
# HR_df_metabric_RFS_coxph_quartiles_CO110_2 <- as.data.frame(HR_list_metabric_RFS_coxph_quartiles$CO110_2)
# HR_df_metabric_RFS_coxph_quartiles_CO110_2$quartile <- gsub('.+quartiles', '', rownames(HR_df_metabric_RFS_coxph_quartiles_CO110_2))
# 
# # save the results
# write.csv(HR_df_metabric_RFS_coxph_quartiles_THR25, 'objs/HR/metabric/RFS/THR25_quartiles_HR.csv')
# write.csv(HR_df_metabric_RFS_coxph_quartiles_CO110_1, 'objs/HR/metabric/RFS/CO110_1_quartiles_HR.csv')
# write.csv(HR_df_metabric_RFS_coxph_quartiles_CO110_2, 'objs/HR/metabric/RFS/CO110_2_quartiles_HR.csv')


########################################################################  
########################################################################  
## by clinical groups

########
## by quartiles

## make a factor with Q1 (lowest risk) being the reference

# CoxData_metabric
#CoxData_metabric$metabric_prob_THR25_quartiles <- factor(CoxData_metabric$metabric_prob_THR25_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_metabric$metabric_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_metabric$metabric_prob_THR25_quartiles))

#CoxData_metabric$metabric_prob_CO110_1_quartiles <- factor(CoxData_metabric$metabric_prob_CO110_1_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_metabric$metabric_prob_CO110_1_quartiles) <- paste0('Q', levels(CoxData_metabric$metabric_prob_CO110_1_quartiles))

#CoxData_metabric$metabric_prob_CO110_2_quartiles <- factor(CoxData_metabric$metabric_prob_CO110_2_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_metabric$metabric_prob_CO110_2_quartiles) <- paste0('Q', levels(CoxData_metabric$metabric_prob_CO110_2_quartiles))

######
# CoxData_metabric_PAM
#CoxData_metabric_PAM$metabric_prob_THR25_quartiles <- factor(CoxData_metabric_PAM$metabric_prob_THR25_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_metabric_PAM$metabric_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_metabric_PAM$metabric_prob_THR25_quartiles))

CoxData_metabric_PAM$metabric_prob_CO110_1_quartiles <- factor(CoxData_metabric_PAM$metabric_prob_CO110_1_quartiles, levels = c('1', '2', '3', '4'))
levels(CoxData_metabric_PAM$metabric_prob_CO110_1_quartiles) <- paste0('Q', levels(CoxData_metabric_PAM$metabric_prob_CO110_1_quartiles))

#CoxData_metabric_PAM$metabric_prob_CO110_2_quartiles <- factor(CoxData_metabric_PAM$metabric_prob_CO110_2_quartiles, levels = c('1', '2', '3', '4'))
#levels(CoxData_metabric_PAM$metabric_prob_CO110_2_quartiles) <- paste0('Q', levels(CoxData_metabric_PAM$metabric_prob_CO110_2_quartiles))

######
# CoxData_metabric_PAM_Q1vsQ4: 
#CoxData_metabric_PAM_Q1vsQ4_THR25$metabric_prob_THR25_quartiles <- factor(CoxData_metabric_PAM_Q1vsQ4_THR25$metabric_prob_THR25_quartiles, levels = c('1', '4'))
#levels(CoxData_metabric_PAM_Q1vsQ4_THR25$metabric_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_metabric_PAM_Q1vsQ4_THR25$metabric_prob_THR25_quartiles))

CoxData_metabric_PAM_Q1vsQ4_CO110_1$metabric_prob_CO110_1_quartiles <- factor(CoxData_metabric_PAM_Q1vsQ4_CO110_1$metabric_prob_CO110_1_quartiles, levels = c('1', '4'))
levels(CoxData_metabric_PAM_Q1vsQ4_CO110_1$metabric_prob_CO110_1_quartiles) <- paste0('Q', levels(CoxData_metabric_PAM_Q1vsQ4_CO110_1$metabric_prob_CO110_1_quartiles))

#CoxData_metabric_PAM_Q1vsQ4_CO110_2$metabric_prob_CO110_2_quartiles <- factor(CoxData_metabric_PAM_Q1vsQ4_CO110_2$metabric_prob_CO110_2_quartiles, levels = c('1', '4'))
#levels(CoxData_metabric_PAM_Q1vsQ4_CO110_2$metabric_prob_CO110_2_quartiles) <- paste0('Q', levels(CoxData_metabric_PAM_Q1vsQ4_CO110_2$metabric_prob_CO110_2_quartiles))


######
# CoxData_metabric_Q1vsQ4
#CoxData_metabric_Q1vsQ4_THR25$metabric_prob_THR25_quartiles <- factor(CoxData_metabric_Q1vsQ4_THR25$metabric_prob_THR25_quartiles, levels = c('1', '4'))
#levels(CoxData_metabric_Q1vsQ4_THR25$metabric_prob_THR25_quartiles) <- paste0('Q', levels(CoxData_metabric_Q1vsQ4_THR25$metabric_prob_THR25_quartiles))

CoxData_metabric_Q1vsQ4_CO110_1$metabric_prob_CO110_1_quartiles <- factor(CoxData_metabric_Q1vsQ4_CO110_1$metabric_prob_CO110_1_quartiles, levels = c('1', '4'))
levels(CoxData_metabric_Q1vsQ4_CO110_1$metabric_prob_CO110_1_quartiles) <- paste0('Q', levels(CoxData_metabric_Q1vsQ4_CO110_1$metabric_prob_CO110_1_quartiles))

#CoxData_metabric_Q1vsQ4_CO110_2$metabric_prob_CO110_2_quartiles <- factor(CoxData_metabric_Q1vsQ4_CO110_2$metabric_prob_CO110_2_quartiles, levels = c('1', '4'))
#levels(CoxData_metabric_Q1vsQ4_CO110_2$metabric_prob_CO110_2_quartiles) <- paste0('Q', levels(CoxData_metabric_Q1vsQ4_CO110_2$metabric_prob_CO110_2_quartiles))


#CoxData_metabric_PAM$Pam50...Claudin.low.subtype <- factor(CoxData_metabric_PAM$Pam50...Claudin.low.subtype)

##############################################
## fit

## PAM50

# os: quartiles: all
# fit on each PAM50 subtype
# lapply(split(CoxData_metabric_PAM, CoxData_metabric_PAM$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR25_quartiles, data = x)))

lapply(split(CoxData_metabric_PAM, CoxData_metabric_PAM$Pam50...Claudin.low.subtype),
       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_1_quartiles, data = x)))

#lapply(split(CoxData_metabric_PAM, CoxData_metabric_PAM$Pam50...Claudin.low.subtype),
#       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_2_quartiles, data = x)))

################################
# os: quartiles: Q1 vs Q4

# fit on each PAM50 subtype

# lapply(split(CoxData_metabric_PAM_Q1vsQ4_THR25, CoxData_metabric_PAM_Q1vsQ4_THR25$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR25_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_CO110_1, CoxData_metabric_PAM_Q1vsQ4_CO110_1$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_1_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_CO110_2, CoxData_metabric_PAM_Q1vsQ4_CO110_2$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_2_quartiles, data = x)))

################################
# rfs: quartiles: all
# fit on each PAM50 subtype
#lapply(split(CoxData_metabric_PAM, CoxData_metabric_PAM$Pam50...Claudin.low.subtype),
#       function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR25_quartiles, data = x)))

lapply(split(CoxData_metabric_PAM, CoxData_metabric_PAM$Pam50...Claudin.low.subtype),
       function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_1_quartiles, data = x)))

#lapply(split(CoxData_metabric_PAM, CoxData_metabric_PAM$Pam50...Claudin.low.subtype),
#       function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_2_quartiles, data = x)))

################################
# rfs: quartiles: Q1 vs Q4

# fit on each PAM50 subtype

# lapply(split(CoxData_metabric_PAM_Q1vsQ4_THR25, CoxData_metabric_PAM_Q1vsQ4_THR25$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR25_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_CO110_1, CoxData_metabric_PAM_Q1vsQ4_CO110_1$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_1_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_CO110_2, CoxData_metabric_PAM_Q1vsQ4_CO110_2$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_2_quartiles, data = x)))



##############################################
## fit

## X3 classifier

# os: quartiles: all
# fit on each X3 subtypes
#lapply(split(CoxData_metabric, CoxData_metabric$X3.Gene.classifier.subtype),
#       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR25_quartiles, data = x)))

lapply(split(CoxData_metabric, CoxData_metabric$X3.Gene.classifier.subtype),
       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_1_quartiles, data = x)))

#lapply(split(CoxData_metabric, CoxData_metabric$X3.Gene.classifier.subtype),
#       function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_2_quartiles, data = x)))

################################
# os: quartiles: Q1 vs Q4

# fit on each X3 subtypes

# lapply(split(CoxData_metabric_PAM_Q1vsQ4_THR25, CoxData_metabric_PAM_Q1vsQ4_THR25$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR25_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_CO110_1, CoxData_metabric_PAM_Q1vsQ4_CO110_1$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_1_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_CO110_2, CoxData_metabric_PAM_Q1vsQ4_CO110_2$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_CO110_2_quartiles, data = x)))

################################
# rfs: quartiles: all
# fit on each X3 subtypes
#lapply(split(CoxData_metabric, CoxData_metabric$X3.Gene.classifier.subtype),
#       function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR25_quartiles, data = x)))

lapply(split(CoxData_metabric, CoxData_metabric$X3.Gene.classifier.subtype),
       function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_1_quartiles, data = x)))

#lapply(split(CoxData_metabric, CoxData_metabric$X3.Gene.classifier.subtype),
#       function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_2_quartiles, data = x)))

################################
# rfs: quartiles: Q1 vs Q4

# fit on each X3 subtypes

# lapply(split(CoxData_metabric_PAM_Q1vsQ4_THR25, CoxData_metabric_PAM_Q1vsQ4_THR25$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_THR25_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_CO110_1, CoxData_metabric_PAM_Q1vsQ4_CO110_1$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_1_quartiles, data = x)))
# 
# lapply(split(CoxData_metabric_PAM_Q1vsQ4_CO110_2, CoxData_metabric_PAM_Q1vsQ4_CO110_2$Pam50...Claudin.low.subtype),
#        function(x) summary(coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_CO110_2_quartiles, data = x)))





















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
Fit_sig_TCGA_os_CO110_1 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_predClasses_CO110_1, data = CoxData_tcga)
Fit_sig_TCGA_os_CO110_2 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_predClasses_CO110_2, data = CoxData_tcga)

## by quartiles
Fit_sig_TCGA_os_THR25_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_THR25_quartiles, data = CoxData_tcga)
Fit_sig_TCGA_os_CO110_1_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO110_1_quartiles, data = CoxData_tcga)
Fit_sig_TCGA_os_CO110_2_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO110_2_quartiles, data = CoxData_tcga)

pdf("./figures/logreg/THR25_Allpairs_TCGA_os.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_THR25,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 25 (logistic regression) and TCGA OS')
dev.off()


pdf("./figures/logreg/CO110_1_Allpairs_TCGA_os.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_CO110_1,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and TCGA OS')
dev.off()


pdf("./figures/logreg/CO110_2_Allpairs_TCGA_os.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_CO110_1,
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

pdf("./figures/logreg/CO110_1_tcga_os_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_CO110_1_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and TCGA OS: quartiles')
dev.off()

pdf("./figures/logreg/CO110_2_tcga_os_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_CO110_2_quartiles,
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

Fit_sig_TCGA_coxph_CO110_1 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO110_1, data = CoxData_tcga)
summary(Fit_sig_TCGA_coxph_CO110_1)

Fit_sig_TCGA_coxph_CO110_2 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_CO110_2, data = CoxData_tcga)
summary(Fit_sig_TCGA_coxph_CO110_2)


png('./figures/logreg/THR25_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_coxph_THR25, fontsize = 0.5)
dev.off()

png('./figures/logreg/CO110_1_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_coxph_CO110_1, fontsize = 0.5)
dev.off()

png('./figures/logreg/CO110_2_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_coxph_CO110_2, fontsize = 0.5)
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
Fit_sig_TCGA_pfs_CO110_1 <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_predClasses_CO110_1, data = CoxData_tcga)
Fit_sig_TCGA_pfs_CO110_2 <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_predClasses_CO110_2, data = CoxData_tcga)

## by quartiles
Fit_sig_TCGA_pfs_THR25_quartiles <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_THR25_quartiles, data = CoxData_tcga)
Fit_sig_TCGA_pfs_CO110_1_quartiles <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_CO110_1_quartiles, data = CoxData_tcga)
Fit_sig_TCGA_pfs_CO110_2_quartiles <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_CO110_2_quartiles, data = CoxData_tcga)

pdf("./figures/logreg/THR25_Allpairs_TCGA_PFS.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_THR25,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 25 (logistic regression) and TCGA PFS')
dev.off()

pdf("./figures/logreg/CO110_1_Allpairs_TCGA_PFS.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_CO110_1,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and TCGA PFS')
dev.off()

pdf("./figures/logreg/CO110_2_Allpairs_TCGA_PFS.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_CO110_2,
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

pdf("./figures/logreg/CO110_1_tcga_pfs_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_CO110_1_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and TCGA PFS: quartiles')
dev.off()

pdf("./figures/logreg/CO110_2_tcga_pfs_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_CO110_2_quartiles,
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

Fit_sig_TCGA_pfs_coxph_CO110_1 <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_CO110_1, data = CoxData_tcga)
summary(Fit_sig_TCGA_pfs_coxph_CO110_1)

Fit_sig_TCGA_pfs_coxph_CO110_2 <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_CO110_2, data = CoxData_tcga)
summary(Fit_sig_TCGA_pfs_coxph_CO110_2)

png('./figures/logreg/THR25_HR_tcga_pfs.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_pfs_coxph_THR25, fontsize = 0.5)
dev.off()

png('./figures/logreg/CO110_1_HR_tcga_pfs.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_pfs_coxph_CO110_1, fontsize = 0.5)
dev.off()

png('./figures/logreg/CO110_2_HR_tcga_pfs.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_pfs_coxph_CO110_2, fontsize = 0.5)
dev.off()



###########################
## heatmap

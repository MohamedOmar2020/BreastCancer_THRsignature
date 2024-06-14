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


#############################################
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

# add X3 subtypes
table(Pheno_metabric$X3.Gene.classifier.subtype)
Data_metabric$X3 <- Pheno_metabric$X3.Gene.classifier.subtype
table(Data_metabric$X3)

# check
all(rownames(Data_metabric) == rownames(Pheno_metabric))

# Ensure that 'os_even' and 'rfs_event' are numeric
Data_metabric$os_event <- as.numeric(as.character(Data_metabric$os_event))
Data_metabric$rfs_event <- as.numeric(as.character(Data_metabric$rfs_event))

#############################################################################################################
##############################################################################################################

## Keep only the relevant information 
survival_metabric <- Data_metabric[, c("os_event", "os_time", 
                                       "rfs_event", "rfs_time", 
                                       'age', 'stage', 'grade',
                                       "PAM50", "ER.Status", "HER2.Status",
                                       "X3")] 

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Divide Basal and Claudin-low using i20------------------------
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## Load the i20 signature ---------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
i20 <- read_xlsx("./figures/c3_DE_THR50_RFS/THR50_c3_longVSshortSurv_DE.xlsx")$gene

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## get TNBC----------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
table(survival_metabric$X3)
TNBC_pheno <- survival_metabric[survival_metabric$X3 == 'ER-/HER2-', ]
TNBC_pheno <- TNBC_pheno[!is.na(TNBC_pheno$X3), ]

TNBC_expr <- Expr_metabric[, rownames(TNBC_pheno)]

all(rownames(TNBC_pheno) == colnames(TNBC_expr))


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## use I20 to separate TNBC into two clusters
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

i20_filt <- intersect(i20, rownames(Expr_metabric))
TNBC_i20_expr <- TNBC_expr[i20_filt, ]


# using logistic regression:
Data_TNBC <- data.frame(cbind(t(TNBC_expr), 'rfs' = TNBC_pheno$os_event))
Data_TNBC$rfs <- as.factor(Data_TNBC$rfs)
table(Data_TNBC$rfs)
i20_logReg_TNBC <- glm(as.formula((paste("rfs ~", paste(i20_filt, collapse = "+")))), data = Data_TNBC, family = "binomial")
TNBC_pheno$i20_logReg_score <- i20_logReg_TNBC %>% predict(Data_TNBC , type = "response")

#&&&&&&&&&&&&&&&&&&
### determine the best threshold-------
#&&&&&&&&&&&&&&&&&&
ROC_thr_TNBC <- coords(roc(TNBC_pheno$rfs_event, TNBC_pheno$i20_logReg_score, direction = "<"), "best")["threshold"]
ROC_thr_TNBC

#&&&&&&&&&&&&&&&&&&
### Assign samples to subclass based on i20 expression-------
#&&&&&&&&&&&&&&&&&&
TNBC_pheno$immune_clusters <- ifelse(TNBC_pheno$i20_logReg_score >= ROC_thr_TNBC$threshold, "ER-/HER2-.i-", "ER-/HER2-.i+")
table(TNBC_pheno$immune_clusters)
TNBC_pheno$immune_clusters <- factor(TNBC_pheno$immune_clusters, levels = c('ER-/HER2-.i-', 'ER-/HER2-.i+'))

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

survival_metabric$Sample.ID <- rownames(survival_metabric)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## recombine with the rest---------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

TNBC <- data.frame(TNBC_immune_clusters = TNBC_pheno$immune_clusters,
                    `Sample.ID` = rownames(TNBC_pheno))

rownames(TNBC) <- rownames(TNBC_pheno)

# merge
survival_metabric$X3[survival_metabric$X3 %in% c('ER-/HER2-')] <- NA


survival_metabric2 <- merge(x = TNBC, y = survival_metabric, by = "Sample.ID", all.y = TRUE)

survival_metabric2 <- survival_metabric2 %>% 
  mutate(X3_merged = coalesce(TNBC_immune_clusters, X3))


table(survival_metabric2$X3_merged)

# keep only the major PAM50 subtypes
survival_metabric2 <- survival_metabric2[!is.na(survival_metabric2$X3_merged), ]


survival_metabric2$X3_merged <- as.factor(survival_metabric2$X3_merged)
table(survival_metabric2$X3_merged)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## survival analysis---------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
### OS---------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# fix the levels

cluster_colors <- c("#4DAF4A", "#E41A1C", "#1960b2", "#984EA3", "#FF7F00")


##################
Fit_metabric_X3_os <- survfit(Surv(os_time, os_event) ~ X3_merged, data = survival_metabric2)

png("./figures/revision/X3_metabric_os_clusters_i20.png", width = 2000, height = 2000, res = 300)
gg <- ggsurvplot(Fit_metabric_X3_os,
                 risk.table = FALSE,
                 pval = TRUE,
                 palette = cluster_colors,
                 xlim = c(0,240),
                 legend.labs = levels(survival_metabric2$X3_merged),
                 legend.title	= '',
                 pval.size = 10,
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
                 #title = 'THR70 clusters and RFS'
) 

gg$plot + guides(
  colour = guide_legend(ncol = 2))
dev.off()

# Cox
cox_metabric$merged_THR_clusters_i20_logReg <- factor(cox_metabric$merged_THR_clusters_i20_logReg, levels = c('E3', 'E2', 'E1', 'HER2+', 'PQNBC_i+', 'PQNBC_i-')) 
cox_Fit_metabric_os_THR70_with_i20_logReg <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ merged_THR_clusters_i20_logReg, data = cox_metabric)
summary(cox_Fit_metabric_os_THR70_with_i20_logReg)



#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
### RFS---------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

##################
Fit_metabric_X3_rfs <- survfit(Surv(rfs_time, rfs_event) ~ X3_merged, data = survival_metabric2)

png("./figures/revision/X3_metabric_rfs_clusters_i20.png", width = 2000, height = 2000, res = 300)
gg <- ggsurvplot(Fit_metabric_X3_rfs,
                 risk.table = FALSE,
                 pval = TRUE,
                 palette = cluster_colors,
                 xlim = c(0,240),
                 legend.labs = levels(survival_metabric2$X3_merged),
                 legend.title	= '',
                 pval.size = 10,
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
                 #title = 'THR70 clusters and RFS'
) 

gg$plot + guides(
  colour = guide_legend(ncol = 3))
dev.off()

# Cox
cox_metabric <- survival_metabric2
table(survival_metabric2$X3_merged)
cox_metabric$X3_merged <- factor(cox_metabric$X3_merged, levels = c('ER-/HER2-.i+','ER-/HER2-.i-', 'ER+/HER2- High Prolif', 'ER+/HER2- Low Prolif', 'HER2+')) 
table(cox_metabric$X3_merged)

cox_Fit_metabric_X3_rfs <- coxph(Surv(rfs_time, rfs_event) ~ X3_merged, data = cox_metabric)
summary(cox_Fit_metabric_X3_rfs)





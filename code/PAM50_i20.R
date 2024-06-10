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


###################################################
# load PAM50
###################################################
PAM50_signature <- readxl::read_xlsx("./data/THR_Signatures_Jan25_2023.xlsx")

# get the THR50 signature
PAM50 <- PAM50_signature$PAM50[!is.na(PAM50_signature$PAM50)]

PAM50 <- gsub('-', '', PAM50)

################
# fix gene names
rownames(Expr_metabric)[grep('^ZNF652', rownames(Expr_metabric))]

# filter the THR signatures to include only the genes present in the expr matrices
PAM50_fil <- PAM50[PAM50 %in% rownames(Expr_metabric)]

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
## get basal----------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
basal_pheno <- survival_metabric[survival_metabric$PAM50 == 'Basal', ]
basal_expr <- Expr_metabric[, rownames(basal_pheno)]

all(rownames(basal_pheno) == colnames(basal_expr))


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## get claudin-low----------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
claudin_pheno <- survival_metabric[survival_metabric$PAM50 == 'claudin-low', ]
claudin_expr <- Expr_metabric[, rownames(claudin_pheno)]

all(rownames(claudin_pheno) == colnames(claudin_expr))

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## use I20 to separate basal into two clusters
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

i20_filt <- intersect(i20, rownames(Expr_metabric))
basal_i20_expr <- basal_expr[i20_filt, ]


# using logistic regression:
Data_basal <- data.frame(cbind(t(basal_expr), 'RFS' = basal_pheno$rfs_event))
Data_basal$RFS <- as.factor(Data_basal$RFS)
table(Data_basal$RFS)
i20_logReg_basal <- glm(as.formula((paste("RFS ~", paste(i20_filt, collapse = "+")))), data = Data_basal, family = "binomial")
basal_pheno$i20_logReg_score <- i20_logReg_basal %>% predict(Data_basal , type = "response")

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## use I20 to separate claudin-low into two clusters
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
claudin_i20_expr <- claudin_expr[i20_filt, ]


# using logistic regression:
Data_claudin <- data.frame(cbind(t(claudin_expr), 'RFS' = claudin_pheno$rfs_event))
Data_claudin$RFS <- as.factor(Data_claudin$RFS)
table(Data_claudin$RFS)
i20_logReg_claudin <- glm(as.formula((paste("RFS ~", paste(i20_filt, collapse = "+")))), data = Data_claudin, family = "binomial")
claudin_pheno$i20_logReg_score <- i20_logReg_claudin %>% predict(Data_claudin , type = "response")

#&&&&&&&&&&&&&&&&&&
### determine the best threshold-------
#&&&&&&&&&&&&&&&&&&
ROC_thr_basal <- coords(roc(basal_pheno$rfs_event, basal_pheno$i20_logReg_score, direction = "<"), "best")["threshold"]
ROC_thr_basal

ROC_thr_claudin <- coords(roc(claudin_pheno$rfs_event, claudin_pheno$i20_logReg_score, direction = "<"), "best")["threshold"]
ROC_thr_claudin



#&&&&&&&&&&&&&&&&&&
### Assign samples to subclass based on i20 expression-------
#&&&&&&&&&&&&&&&&&&
basal_pheno$immune_clusters <- ifelse(basal_pheno$i20_logReg_score >= ROC_thr_basal$threshold, "Basal_i-", "Basal_i+")
table(basal_pheno$immune_clusters)
basal_pheno$immune_clusters <- factor(basal_pheno$immune_clusters, levels = c('Basal_i-', 'Basal_i+'))


claudin_pheno$immune_clusters <- ifelse(claudin_pheno$i20_logReg_score >= ROC_thr_claudin$threshold, "Claudin-low_i-", "Claudin-low_i+")
table(claudin_pheno$immune_clusters)
claudin_pheno$immune_clusters <- factor(claudin_pheno$immune_clusters, levels = c('Claudin-low_i-', 'Claudin-low_i+'))

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

survival_metabric$Sample.ID <- rownames(survival_metabric)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## recombine with the rest---------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

basal <- data.frame(basal_immune_clusters = basal_pheno$immune_clusters,
                 `Sample.ID` = rownames(basal_pheno))

rownames(basal) <- rownames(basal_pheno)


Claudin <- data.frame(claudin_immune_clusters = claudin_pheno$immune_clusters,
                    `Sample.ID` = rownames(claudin_pheno))

rownames(Claudin) <- rownames(claudin_pheno)


# merge
survival_metabric$PAM50[survival_metabric$PAM50 %in% c('Basal', 'claudin-low')] <- NA


survival_metabric2 <- merge(x = basal, y = survival_metabric, by = "Sample.ID", all.y = TRUE)
survival_metabric2 <- merge(x = Claudin, y = survival_metabric2, by = "Sample.ID", all.y = TRUE)

survival_metabric2 <- survival_metabric2 %>% 
  mutate(PAM50_merged = coalesce(basal_immune_clusters, claudin_immune_clusters, PAM50))


table(survival_metabric2$PAM50_merged)

# keep only the major PAM50 subtypes
survival_metabric2 <- survival_metabric2 %>%
  filter(PAM50_merged %in% c('Basal_i-', 'Basal_i+', 'Claudin-low_i-', 'Claudin-low_i+', 'Her2', 'LumA', 'LumB'))
  

survival_metabric2$PAM50_merged <- as.factor(survival_metabric2$PAM50_merged)
table(survival_metabric2$PAM50_merged)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## survival analysis---------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
### OS---------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# fix the levels

cluster_colors <- c("#1960b2", "#DF271F", "#66A61E", "#CD4174", "#A253A2", "#9D696C", "#1B9E77")


##################
Fit_metabric_PAM50_os <- survfit(Surv(os_time, os_event) ~ PAM50_merged, data = survival_metabric2)

png("./figures/revision/PAM50_metabric_os_clusters_i20.png", width = 2000, height = 2000, res = 300)
gg <- ggsurvplot(Fit_metabric_PAM50_os,
                 risk.table = FALSE,
                 pval = TRUE,
                 palette = cluster_colors,
                 xlim = c(0,240),
                 legend.labs = levels(survival_metabric2$PAM50_merged),
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
Fit_metabric_PAM50_rfs <- survfit(Surv(rfs_time, rfs_event) ~ PAM50_merged, data = survival_metabric2)

png("./figures/revision/PAM50_metabric_rfs_clusters_i20.png", width = 2000, height = 2000, res = 300)
gg <- ggsurvplot(Fit_metabric_PAM50_rfs,
                 risk.table = FALSE,
                 pval = TRUE,
                 palette = cluster_colors,
                 xlim = c(0,240),
                 legend.labs = levels(survival_metabric2$PAM50_merged),
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
cox_metabric <- survival_metabric
cox_metabric$merged_THR_clusters_i20_logReg <- factor(cox_metabric$merged_THR_clusters_i20_logReg, levels = c('PQNBC_i+','PQNBC_i-', 'E3', 'E2', 'E1', 'HER2+')) 

cox_Fit_metabric_rfs_THR70_with_i20_logReg <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ merged_THR_clusters_i20_logReg, data = cox_metabric)
summary(cox_Fit_metabric_rfs_THR70_with_i20_logReg)


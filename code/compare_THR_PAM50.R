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

#################
THR_signature <- readxl::read_xlsx("./data/THR Signatures_sep23.xlsx")

# get THR 25 and 50
THR50 <- THR_signature$`THR-50.1`[!is.na(THR_signature$`THR-50.1`)]

THR50 <- gsub('-', '', THR50)

################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')

##############################################
# fix gene names
##############################################
# Fix in TCGA : ALL GOOD
setdiff(THR50, rownames(Expr_tcga))
grep('^FAM63A', rownames(Expr_tcga), value = TRUE) # MINDY1
grep('^FAM176A', rownames(Expr_tcga), value = TRUE) # EVA1A
grep('^LEPREL1', rownames(Expr_tcga), value = TRUE) # P3H2
grep('^DULLARD', rownames(Expr_tcga), value = TRUE) # SDHAF3

rownames(Expr_tcga)[rownames(Expr_tcga) == 'ACN9'] <- 'SDHAF3'
rownames(Expr_tcga)[rownames(Expr_tcga) == 'FAM176A'] <- 'EVA1A'
rownames(Expr_tcga)[rownames(Expr_tcga) == 'LEPREL1'] <- 'P3H2'
rownames(Expr_tcga)[rownames(Expr_tcga) == 'FAM63A'] <- 'MINDY1'

setdiff(THR50, rownames(Expr_tcga))

###############
# Fix in metabric: 4 missing
setdiff(THR50, rownames(Expr_metabric))

rownames(Expr_metabric)[rownames(Expr_metabric) == 'FAM63A'] <- 'MINDY1'
rownames(Expr_metabric)[rownames(Expr_metabric) == 'FAM176A'] <- 'EVA1A'
rownames(Expr_metabric)[rownames(Expr_metabric) == 'ACN9'] <- 'SDHAF3'
rownames(Expr_metabric)[rownames(Expr_metabric) == 'LEPREL1'] <- 'P3H2'

setdiff(THR50, rownames(Expr_metabric_refAll))


##############
# filter the signatures to include only the genes present in the expr matrices
THR50_fil <- THR50[THR50 %in% rownames(Expr_tcga) & THR50 %in% rownames(Expr_metabric)]

setdiff(THR50, THR50_fil)

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

# check
all(rownames(Data_metabric) == rownames(Pheno_metabric))

# Ensure that 'os_even' and 'rfs_event' are numeric
Data_metabric$os_event <- as.numeric(as.character(Data_metabric$os_event))
Data_metabric$rfs_event <- as.numeric(as.character(Data_metabric$rfs_event))


#############################################################
# compare THR50 to PAM50
#############################################################

# Calculate PAM50 centroids
pam50_centroids <- Data_metabric %>%
  filter(PAM50 %in% c('Basal', 'Her2', 'LumA', 'LumB', 'claudin-low')) %>%
  group_by(PAM50) %>%
  summarize(across(all_of(PAM50_fil), median, na.rm = TRUE), .groups = 'drop')

pam50_centroids <- as.data.frame(pam50_centroids)

rownames(pam50_centroids) <- c('Basal', 'Her2', 'LumA', 'LumB', 'Claudin_low')

# Calculate THR50 centroids
thr50_centroids <- Data_metabric %>%
  filter(PAM50 %in% c('Basal', 'Her2', 'LumA', 'LumB', 'Claudin_low')) %>%  # Using PAM50 groups for simplicity
  group_by(PAM50) %>%
  summarize(across(all_of(THR50_fil), median, na.rm = TRUE), .groups = 'drop')


#######################
# classify samples

genes_in_centroid <- colnames(pam50_centroids)[-1]  # assuming first column is not a gene
genes_to_use <- genes_in_centroid[genes_in_centroid %in% colnames(Data_metabric)]

# Update classify_samples function to handle gene filtering and reordering internally
classify_samples <- function(data, centroids, genes) {
  apply(data[, genes, drop = FALSE], 1, function(x) {
    # Calculate correlations with each centroid
    corrs <- sapply(seq_len(nrow(centroids)), function(y) {
      centroid_genes <- as.numeric(centroids[y, genes, drop = FALSE])  # Ensure numeric data for correlation
      cor(x, centroid_genes, method = "spearman", use = "complete.obs")  # Handling NAs
    })
    # Return the row name (subtype) of the centroid with the highest correlation
    rownames(centroids)[which.max(corrs)]
  })
}


# Classify each sample in Data_metabric for PAM50
Data_metabric$pam50_class <- classify_samples(
  Data_metabric,
  pam50_centroids,
  genes_to_use
)

table(Data_metabric$pam50_class)


# Classify each sample in Data_metabric for THR50
Data_metabric$THR50_class <- classify_samples(
  Data_metabric,
  thr50_centroids,
  THR50_fil
)

table(Data_metabric$THR50_class)

#############################################################################
# Kaplan-Meier for PAM50 and THR50: OS
#############################################################################

# PAM50
km_os_fit_pam50 <- survfit(Surv(os_time, os_event) ~ pam50_class, data = Data_metabric)
tiff("./figures/revision/THR_vs_PAM50/PAM50_metabric_os_centroids.tiff", width = 3000, height = 2500, res = 350)
ggsurvplot(km_os_fit_pam50,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 14,
           short.panel.labs = T,
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
                           legend.title = element_text(size = 16, face = 'bold.italic', color = 'black'),
                           strip.text = element_text(size = 16, face = 'bold.italic', color = 'black')
           ),
           palette = 'jco',
           legend.title	= 'PAM50 Groups',
           legend.labs = c('Basal', 'Claudin-low', 'HER2', 'LumA', 'LumB'),
           xlim = c(0,240),
           #break.x.by = 40,
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC OS by PAM50 subtypes: quartiles'
)
dev.off()

#####################
# THR50
km_os_fit_thr50 <- survfit(Surv(os_time, os_event) ~ THR50_class, data = Data_metabric)
tiff("./figures/revision/THR_vs_PAM50/THR50_metabric_os_centroids.tiff", width = 3000, height = 2500, res = 350)
ggsurvplot(km_os_fit_thr50,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 14,
           short.panel.labs = T,
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
                           legend.title = element_text(size = 16, face = 'bold.italic', color = 'black'),
                           strip.text = element_text(size = 16, face = 'bold.italic', color = 'black')
           ),
           palette = 'jco',
           legend.title	= 'THR50 Groups',
           legend.labs = c('Group 1', 'Group 2', 'Group 3', 'Group 4'),
           xlim = c(0,240),
           #break.x.by = 40,
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC OS by PAM50 subtypes: quartiles'
)
dev.off()



#############################################################################
# Kaplan-Meier for PAM50 and THR50: RFS
#############################################################################

# PAM50
km_rfs_fit_pam50 <- survfit(Surv(rfs_time, rfs_event) ~ pam50_class, data = Data_metabric)
tiff("./figures/revision/THR_vs_PAM50/PAM50_metabric_rfs_centroids.tiff", width = 3000, height = 2500, res = 350)
ggsurvplot(km_rfs_fit_pam50,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 14,
           short.panel.labs = T,
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
                           legend.title = element_text(size = 16, face = 'bold.italic', color = 'black'),
                           strip.text = element_text(size = 16, face = 'bold.italic', color = 'black')
           ),
           palette = 'jco',
           legend.title	= 'PAM50 Groups',
           legend.labs = c('Basal', 'Claudin-low', 'HER2', 'LumA', 'LumB'),
           xlim = c(0,240),
           #break.x.by = 40,
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC OS by PAM50 subtypes: quartiles'
)
dev.off()

#####################
# THR50
km_rfs_fit_thr50 <- survfit(Surv(rfs_time, rfs_event) ~ THR50_class, data = Data_metabric)
tiff("./figures/revision/THR_vs_PAM50/THR50_metabric_rfs_centroids.tiff", width = 3000, height = 2500, res = 350)
ggsurvplot(km_rfs_fit_thr50,
           risk.table = FALSE,
           pval = TRUE,
           pval.size = 14,
           short.panel.labs = T,
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
                           legend.title = element_text(size = 16, face = 'bold.italic', color = 'black'),
                           strip.text = element_text(size = 16, face = 'bold.italic', color = 'black')
           ),
           palette = 'jco',
           legend.title	= 'THR50 Groups',
           legend.labs = c('Group 1', 'Group 2', 'Group 3', 'Group 4'),
           xlim = c(0,240),
           #break.x.by = 40,
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 and METABRIC OS by PAM50 subtypes: quartiles'
)
dev.off()

###############################################################################
# Cox Proportional Hazards model for PAM50
cox_pam50 <- coxph(Surv(os_time, os_event) ~ pam50_class, data = Data_metabric)
summary(cox_pam50)

# Cox Proportional Hazards model for THR50
cox_thr50 <- coxph(Surv(os_time, os_event) ~ THR50_class, data = Data_metabric)
summary(cox_thr50)




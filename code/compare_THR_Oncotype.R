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
grep('^BAG', rownames(Expr_tcga), value = TRUE) # MINDY1
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

grep('^BAG', rownames(Expr_metabric), value = TRUE) # MINDY1
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
# load oncotype Dx
###################################################
oncotype <-  c("MKI67", "AURKA", "BIRC5", "CCNB1", "MYBL2", "MMP11", "CTSL2", 
               "GRB7", "ERBB2", 
               "ESR1", "PGR", "BCL2", "SCUBE2", "GSTM1", "BAG1", "CD68",
               "ACTB", "GAPDH", "GUSB", "RPLP0", "TFRC")

##############
# filter the signatures to include only the genes present in the expr matrices
oncotype_fil <- oncotype[oncotype %in% rownames(Expr_tcga) & oncotype %in% rownames(Expr_metabric)]

setdiff(oncotype, oncotype_fil)

oncotype_fil_clean <- oncotype_fil

# Replace special characters with an underscore
oncotype_fil_clean <- gsub("[^[:alnum:]_]", "_", oncotype_fil_clean)

# Ensure names are valid R variable names
oncotype_fil_clean <- make.names(oncotype_fil_clean, unique = TRUE)

# Optional: Inspect the cleaned gene names
print(oncotype_fil_clean)

# check if any names were changed
if (!all(oncotype_fil_clean == oncotype_fil)) {
  warning("Some gene names were modified to ensure validity")
}
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
# Normalization
#############################################################

# Calculate the median expression of the reference genes
reference_expr <- rowMeans(Data_metabric[, c("ACTB", "GAPDH", "GUSB", "RPLP0", "TFRC")], na.rm = TRUE)

# Normalize the expression of each gene
Data_metabric[names(oncotype_fil_clean)] <- sweep(Data_metabric[names(oncotype_fil_clean)], 1, reference_expr, "/")


# Compute group scores
Data_metabric$GRB7_group_score <- pmax(0.9 * Data_metabric$GRB7 + 0.1 * Data_metabric$ERBB2, 8)
Data_metabric$ER_group_score <- (0.8 * Data_metabric$ESR1 + 1.2 * Data_metabric$PGR + Data_metabric$BCL2 + Data_metabric$SCUBE2) / 4
Data_metabric$proliferation_group_score <- pmax((Data_metabric$MKI67 + Data_metabric$AURKA + Data_metabric$BIRC5 + Data_metabric$CCNB1 + Data_metabric$MYBL2) / 5, 6.5)
Data_metabric$invasion_group_score <- (Data_metabric$CTSL2 + Data_metabric$MMP11) / 2

summary(Data_metabric$GRB7_group_score)
summary(Data_metabric$ER_group_score)
summary(Data_metabric$proliferation_group_score)
summary(Data_metabric$invasion_group_score)
summary(Data_metabric$CD68)
summary(Data_metabric$GSTM1)
summary(Data_metabric$BAG1)

# Compute the Unscaled Recurrence Score (RSU)
Data_metabric$RSU <- 0.47 * Data_metabric$GRB7_group_score - 
  0.34 * Data_metabric$ER_group_score + 
  1.04 * Data_metabric$proliferation_group_score + 
  0.10 * Data_metabric$invasion_group_score +
  0.05 * Data_metabric$CD68 -
  0.08 * Data_metabric$GSTM1 



# Rescale the Recurrence Score (RS)
Data_metabric$RS <- pmin(100, pmax(0, round(20 * (Data_metabric$RSU - 6.7), 0)))
summary(Data_metabric$RS)

# Classify the Recurrence Scores into Risk Categories
Data_metabric$risk_group <- cut(Data_metabric$RS, 
                                breaks = c(-Inf, 18, 30, Inf), 
                                labels = c("Low Risk", "Intermediate Risk", "High Risk"),
                                include.lowest = TRUE)

#############################################################################
# Kaplan-Meier for oncotype Dx
#############################################################################

# Fit a Kaplan-Meier survival curve
km_fit <- survfit(Surv(time = Data_metabric$os_time, event = Data_metabric$os_event) ~ risk_group, data = Data_metabric)

# Plot the survival curves
ggsurvplot(km_fit, 
           data = Data_metabric, 
           palette = "jco",
           ggtheme = theme_minimal(),
           #pval = TRUE,                     # Show p-value of the log-rank test
           #risk.table = TRUE,               # Add risk table
           legend.title = "Risk Group",     # Legend title
           xlab = "Time (months)",          # X-axis label
           ylab = "Survival Probability"    # Y-axis label
)





###### 
# Clean Work space
rm(list = ls())

############################################################################
### Load library
require(switchBox)
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
library(readxl)

################
# Load the METABRIC expr and pheno data
Expr_metabric <- read.delim("./data/brca_metabric_cbioportal/data_mrna_agilent_microarray_zscores_ref_all_samples.txt")
Pheno_metabric <- read.delim("./data/brca_metabric_cbioportal/brca_metabric_clinical_data.tsv")

################
# Load the TCGA expression and Phenotype data
Expr_tcga <- read.delim("./data/brca_tcga/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")
Pheno_tcga <- read.delim("./data/brca_tcga/data_clinical_patient.txt")
Pheno_tcga <- Pheno_tcga[-c(1:4), ]

###########################
## Annotation: metabric
head(rownames(Expr_metabric))
Expr_metabric <- Expr_metabric[!duplicated(Expr_metabric$Hugo_Symbol), ]
rownames(Expr_metabric) <- Expr_metabric$Hugo_Symbol
Expr_metabric$Hugo_Symbol <- NULL
Expr_metabric$Entrez_Gene_Id <- NULL
summary(is.na(rownames(Expr_metabric)))
sel <- which(apply(Expr_metabric, 1, function(x) all(is.finite(x)) ))
Expr_metabric <- Expr_metabric[sel, ]
Expr_metabric <- Expr_metabric[!is.na(rownames(Expr_metabric)),]
dim(Expr_metabric)
range(Expr_metabric)

#############################
## Annotation: TCGA
head(rownames(Expr_tcga))
Expr_tcga <- Expr_tcga[!duplicated(Expr_tcga$Hugo_Symbol), ]
rownames(Expr_tcga) <- Expr_tcga$Hugo_Symbol
Expr_tcga$Hugo_Symbol <- NULL
Expr_tcga$Entrez_Gene_Id <- NULL
summary(is.na(rownames(Expr_tcga)))
sel <- which(apply(Expr_tcga, 1, function(x) all(is.finite(x)) ))
Expr_tcga <- Expr_tcga[sel, ]
Expr_tcga <- Expr_tcga[!is.na(rownames(Expr_tcga)),]
#rownames(Expr_tcga[Expr_tcga < -30, ])
Expr_tcga <- Expr_tcga[!(rownames(Expr_tcga) %in% c('CEACAM18', 'KRTAP12-4', 'RPS4Y2', 'TTTY4C', 'DEFB134')), ]
dim(Expr_tcga)
range(Expr_tcga)

# logscale
#Expr_tcga <- log2(Expr_tcga + 6)

# fix the column names
colnames(Expr_tcga) <- gsub('\\.01', '', colnames(Expr_tcga))

####################################
### Modify the Phenotype data: metabric

Pheno_metabric$Sample.ID <- gsub("\\-", "\\.", Pheno_metabric$Sample.ID)
rownames(Pheno_metabric) <- Pheno_metabric$Sample.ID
CommonSamples <- intersect(colnames(Expr_metabric), rownames(Pheno_metabric))
Pheno_metabric <- Pheno_metabric[CommonSamples, ]

Pheno_metabric$Relapse.Free.Status <- gsub("\\:.+", "", Pheno_metabric$Relapse.Free.Status)
Pheno_metabric$Overall.Survival.Status <- gsub("\\:.+", "", Pheno_metabric$Overall.Survival.Status)

all(rownames(Pheno_metabric) == colnames(Expr_metabric))

Expr_metabric <- as.matrix(Expr_metabric)

############
## Keep only the relevant information 
Pheno_metabric$Relapse.Free.Status <- as.numeric(Pheno_metabric$Relapse.Free.Status)
Pheno_metabric$Overall.Survival.Status <- as.numeric(Pheno_metabric$Overall.Survival.Status)

table(Pheno_metabric$Relapse.Free.Status)
table(Pheno_metabric$Overall.Survival.Status)

Pheno_metabric$Relapse.Free.Status..Months. <- as.numeric(Pheno_metabric$Relapse.Free.Status..Months.)
Pheno_metabric$Overall.Survival..Months. <- as.numeric(Pheno_metabric$Overall.Survival..Months.)

group_metabric <- as.factor(Pheno_metabric$Overall.Survival.Status)
table(group_metabric)

################
## Modify the Phenotype data: tcga

Pheno_tcga$X.Patient.Identifier <- gsub("\\-", "\\.", Pheno_tcga$X.Patient.Identifier)
rownames(Pheno_tcga) <- Pheno_tcga$X.Patient.Identifier
CommonSamples <- intersect(colnames(Expr_tcga), rownames(Pheno_tcga))
Pheno_tcga <- Pheno_tcga[CommonSamples, ]

Pheno_tcga$Progression.Free.Status <- gsub("\\:.+", "", Pheno_tcga$Progression.Free.Status)
Pheno_tcga$Overall.Survival.Status <- gsub("\\:.+", "", Pheno_tcga$Overall.Survival.Status)

all(rownames(Pheno_tcga) == colnames(Expr_tcga))

Expr_tcga <- as.matrix(Expr_tcga)

################
## Keep only the relevant information (Metastasis Event and Time)
Pheno_tcga$Progression.Free.Status <- as.numeric(Pheno_tcga$Progression.Free.Status)
Pheno_tcga$Overall.Survival.Status <- as.numeric(Pheno_tcga$Overall.Survival.Status)

table(Pheno_tcga$Progression.Free.Status)
table(Pheno_tcga$Overall.Survival.Status)

Pheno_tcga$Progress.Free.Survival..Months. <- as.numeric(Pheno_tcga$Progress.Free.Survival..Months.)
Pheno_tcga$Overall.Survival..Months. <- as.numeric(Pheno_tcga$Overall.Survival..Months.)

group_tcga <- as.factor(Pheno_tcga$Overall.Survival.Status)
table(group_tcga)

####################################################################################
# save
save(Expr_metabric, Expr_tcga, group_metabric, group_tcga, Pheno_metabric, Pheno_tcga, file = './objs/forKTSP.rda')









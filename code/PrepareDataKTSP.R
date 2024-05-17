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
Expr_metabric_refAll <- read.delim("./data/brca_metabric_cbioportal/data_mrna_agilent_microarray_zscores_ref_all_samples.txt")
Expr_metabric_refDip <- read.delim("./data/brca_metabric_cbioportal/data_mrna_agilent_microarray_zscores_ref_diploid_samples.txt")
Expr_metabric <- read.delim("./data/brca_metabric_cbioportal/data_mrna_agilent_microarray.txt")

#cancer_study_identifier: brca_metabric
#genetic_alteration_type: MRNA_EXPRESSION
#datatype: CONTINUOUS
#stable_id: mrna
#show_profile_in_analysis_tab: false
#profile_name: mRNA expression (microarray)
#profile_description: Expression log intensity levels (Illumina Human v3 microarray)
#data_filename: data_mrna_agilent_microarray.txt


Pheno_metabric <- read.delim("./data/brca_metabric_cbioportal/brca_metabric_clinical_data.tsv")

################
# Load the TCGA expression and Phenotype data
Expr_tcga_refAll <- read.delim("./data/brca_tcga/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")
Expr_tcga_refDip <- read.delim("./data/brca_tcga/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt")
#Expr_tcga_refNorm <- read.delim("./data/brca_tcga/data_mrna_seq_v2_rsem_zscores_ref_normal_samples.txt")
Expr_tcga <- read.delim("./data/brca_tcga/data_mrna_seq_v2_rsem.txt")

#cancer_study_identifier: brca_tcga_pan_can_atlas_2018
#stable_id: rna_seq_v2_mrna
#genetic_alteration_type: MRNA_EXPRESSION
#datatype: CONTINUOUS
#show_profile_in_analysis_tab: FALSE
#profile_name: mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)
#profile_description:mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)
#data_filename: data_mrna_seq_v2_rsem.txt

Pheno_tcga <- read.delim("./data/brca_tcga/brca_tcga_clinical_data.tsv")
#Pheno_tcga <- Pheno_tcga[-c(1:4), ]


summary(Pheno_metabric$Overall.Survival..Months.)
summary(Pheno_metabric$Relapse.Free.Status..Months.)

summary(Pheno_tcga$Overall.Survival..Months.)
summary(Pheno_tcga$Disease.Free..Months)


######################################################
## Annotation: metabric
######################################################
# ref all
head(rownames(Expr_metabric_refAll))
Expr_metabric_refAll <- Expr_metabric_refAll[!duplicated(Expr_metabric_refAll$Hugo_Symbol), ]
rownames(Expr_metabric_refAll) <- Expr_metabric_refAll$Hugo_Symbol
Expr_metabric_refAll$Hugo_Symbol <- NULL
Expr_metabric_refAll$Entrez_Gene_Id <- NULL
summary(is.na(rownames(Expr_metabric_refAll)))
sel <- which(apply(Expr_metabric_refAll, 1, function(x) all(is.finite(x)) ))
Expr_metabric_refAll <- Expr_metabric_refAll[sel, ]
Expr_metabric_refAll <- Expr_metabric_refAll[!is.na(rownames(Expr_metabric_refAll)),]
dim(Expr_metabric_refAll)
range(Expr_metabric_refAll)

##############
## ref diploid
head(rownames(Expr_metabric_refDip))
Expr_metabric_refDip <- Expr_metabric_refDip[!duplicated(Expr_metabric_refDip$Hugo_Symbol), ]
rownames(Expr_metabric_refDip) <- Expr_metabric_refDip$Hugo_Symbol
Expr_metabric_refDip$Hugo_Symbol <- NULL
Expr_metabric_refDip$Entrez_Gene_Id <- NULL
summary(is.na(rownames(Expr_metabric_refDip)))
sel <- which(apply(Expr_metabric_refDip, 1, function(x) all(is.finite(x)) ))
Expr_metabric_refDip <- Expr_metabric_refDip[sel, ]
Expr_metabric_refDip <- Expr_metabric_refDip[!is.na(rownames(Expr_metabric_refDip)),]
dim(Expr_metabric_refDip)
range(Expr_metabric_refDip)

##############
## raw
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
range(Expr_metabric) # log-scaled

############################################################
## Annotation: TCGA
##############################################################
# ref all
head(rownames(Expr_tcga_refAll))
Expr_tcga_refAll <- Expr_tcga_refAll[!duplicated(Expr_tcga_refAll$Hugo_Symbol), ]
rownames(Expr_tcga_refAll) <- Expr_tcga_refAll$Hugo_Symbol
Expr_tcga_refAll$Hugo_Symbol <- NULL
Expr_tcga_refAll$Entrez_Gene_Id <- NULL
summary(is.na(rownames(Expr_tcga_refAll)))
sel <- which(apply(Expr_tcga_refAll, 1, function(x) all(is.finite(x)) ))
Expr_tcga_refAll <- Expr_tcga_refAll[sel, ]
Expr_tcga_refAll <- Expr_tcga_refAll[!is.na(rownames(Expr_tcga_refAll)),]
#rownames(Expr_tcga_refAll[Expr_tcga_refAll < -30, ])
Expr_tcga_refAll <- Expr_tcga_refAll[!(rownames(Expr_tcga_refAll) %in% c('CEACAM18', 'KRTAP12-4', 'RPS4Y2', 'TTTY4C', 'DEFB134')), ]
dim(Expr_tcga_refAll)
range(Expr_tcga_refAll)
# logscale
#Expr_tcga_refAll <- log2(Expr_tcga_refAll + 6)
# fix the column names
colnames(Expr_tcga_refAll) <- gsub('\\.01', '', colnames(Expr_tcga_refAll))

###############
# ref diploid
head(rownames(Expr_tcga_refDip))
Expr_tcga_refDip <- Expr_tcga_refDip[!duplicated(Expr_tcga_refDip$Hugo_Symbol), ]
rownames(Expr_tcga_refDip) <- Expr_tcga_refDip$Hugo_Symbol
Expr_tcga_refDip$Hugo_Symbol <- NULL
Expr_tcga_refDip$Entrez_Gene_Id <- NULL
summary(is.na(rownames(Expr_tcga_refDip)))
sel <- which(apply(Expr_tcga_refDip, 1, function(x) all(is.finite(x)) ))
Expr_tcga_refDip <- Expr_tcga_refDip[sel, ]
Expr_tcga_refDip <- Expr_tcga_refDip[!is.na(rownames(Expr_tcga_refDip)),]
#rownames(Expr_tcga_refDip[Expr_tcga_refDip < -30, ])
Expr_tcga_refDip <- Expr_tcga_refDip[!(rownames(Expr_tcga_refDip) %in% c('CEACAM18', 'KRTAP12-4', 'RPS4Y2', 'TTTY4C', 'DEFB134')), ]
dim(Expr_tcga_refDip)
range(Expr_tcga_refDip)
# logscale
#Expr_tcga_refDip <- log2(Expr_tcga_refDip + 6)
# fix the column names
colnames(Expr_tcga_refDip) <- gsub('\\.01', '', colnames(Expr_tcga_refDip))

###############
# ref normal
head(rownames(Expr_tcga_refNorm))
Expr_tcga_refNorm <- Expr_tcga_refNorm[!duplicated(Expr_tcga_refNorm$Hugo_Symbol), ]
rownames(Expr_tcga_refNorm) <- Expr_tcga_refNorm$Hugo_Symbol
Expr_tcga_refNorm$Hugo_Symbol <- NULL
Expr_tcga_refNorm$Entrez_Gene_Id <- NULL
summary(is.na(rownames(Expr_tcga_refNorm)))
sel <- which(apply(Expr_tcga_refNorm, 1, function(x) all(is.finite(x)) ))
Expr_tcga_refNorm <- Expr_tcga_refNorm[sel, ]
Expr_tcga_refNorm <- Expr_tcga_refNorm[!is.na(rownames(Expr_tcga_refNorm)),]
#rownames(Expr_tcga_refNorm[Expr_tcga_refNorm < -30, ])
Expr_tcga_refNorm <- Expr_tcga_refNorm[!(rownames(Expr_tcga_refNorm) %in% c('CEACAM18', 'KRTAP12-4', 'RPS4Y2', 'TTTY4C', 'DEFB134')), ]
dim(Expr_tcga_refNorm)
range(Expr_tcga_refNorm)
# logscale
#Expr_tcga_refNorm <- log2(Expr_tcga_refNorm + 6)
# fix the column names
colnames(Expr_tcga_refNorm) <- gsub('\\.01', '', colnames(Expr_tcga_refNorm))

###############
# raw
head(rownames(Expr_tcga))
Expr_tcga <- Expr_tcga[!duplicated(Expr_tcga$Hugo_Symbol), ]
rownames(Expr_tcga) <- Expr_tcga$Hugo_Symbol
Expr_tcga$Hugo_Symbol <- NULL
Expr_tcga$Entrez_Gene_Id <- NULL
summary(is.na(rownames(Expr_tcga)))
Expr_tcga <- Expr_tcga[-1,]
sel <- which(apply(Expr_tcga, 1, function(x) all(is.finite(x)) ))
Expr_tcga <- Expr_tcga[sel, ]
Expr_tcga <- Expr_tcga[!is.na(rownames(Expr_tcga)),]
#rownames(Expr_tcga[Expr_tcga < -30, ])
Expr_tcga <- Expr_tcga[!(rownames(Expr_tcga) %in% c('CEACAM18', 'KRTAP12-4', 'RPS4Y2', 'TTTY4C', 'DEFB134')), ]
dim(Expr_tcga)
range(Expr_tcga)
# logscale
Expr_tcga <- log2(Expr_tcga + 1)
# fix the column names
colnames(Expr_tcga) <- gsub('\\.01', '', colnames(Expr_tcga))

########################################################################
### Modify the Phenotype data: metabric
########################################################################
Pheno_metabric$Sample.ID <- gsub("\\-", "\\.", Pheno_metabric$Sample.ID)
rownames(Pheno_metabric) <- Pheno_metabric$Sample.ID
CommonSamples <- intersect(colnames(Expr_metabric_refAll), rownames(Pheno_metabric))
Pheno_metabric <- Pheno_metabric[CommonSamples, ]

Pheno_metabric$Relapse.Free.Status <- gsub("\\:.+", "", Pheno_metabric$Relapse.Free.Status)
Pheno_metabric$Overall.Survival.Status <- gsub("\\:.+", "", Pheno_metabric$Overall.Survival.Status)

all(rownames(Pheno_metabric) == colnames(Expr_metabric_refAll))
all(rownames(Pheno_metabric) == colnames(Expr_metabric_refDip))
all(rownames(Pheno_metabric) == colnames(Expr_metabric))

#Expr_metabric_refAll <- as.matrix(Expr_metabric_refAll)
#Expr_metabric_refDip <- as.matrix(Expr_metabric_refDip)

## Keep only the relevant information 
Pheno_metabric$Relapse.Free.Status <- as.numeric(Pheno_metabric$Relapse.Free.Status)
Pheno_metabric$Overall.Survival.Status <- as.numeric(Pheno_metabric$Overall.Survival.Status)

table(Pheno_metabric$Relapse.Free.Status)
table(Pheno_metabric$Overall.Survival.Status)

Pheno_metabric$Relapse.Free.Status..Months. <- as.numeric(Pheno_metabric$Relapse.Free.Status..Months.)
Pheno_metabric$Overall.Survival..Months. <- as.numeric(Pheno_metabric$Overall.Survival..Months.)

group_metabric <- as.factor(Pheno_metabric$Overall.Survival.Status)
table(group_metabric)

########################################################################################
## Modify the Phenotype data: tcga
########################################################################################
Pheno_tcga$Patient.ID <- gsub("\\-", "\\.", Pheno_tcga$Patient.ID)
#rownames(Pheno_tcga) <- Pheno_tcga$Patient.ID
Pheno_tcga <- Pheno_tcga[!duplicated(Pheno_tcga$Patient.ID), ]
rownames(Pheno_tcga) <- Pheno_tcga$Patient.ID
CommonSamples <- intersect(colnames(Expr_tcga_refAll), rownames(Pheno_tcga))
Pheno_tcga <- Pheno_tcga[CommonSamples, ]

Pheno_tcga$Disease.Free.Status <- gsub("\\:.+", "", Pheno_tcga$Disease.Free.Status)
Pheno_tcga$Overall.Survival.Status <- gsub("\\:.+", "", Pheno_tcga$Overall.Survival.Status)

all(rownames(Pheno_tcga) == colnames(Expr_tcga_refAll))
all(rownames(Pheno_tcga) == colnames(Expr_tcga_refDip))
all(rownames(Pheno_tcga) == colnames(Expr_tcga))

#Expr_tcga_refAll <- as.matrix(Expr_tcga_refAll)
#Expr_tcga_refDip <- as.matrix(Expr_tcga_refDip)
#Expr_tcga <- as.matrix(Expr_tcga)

## Keep only the relevant information (Metastasis Event and Time)
table(Pheno_tcga$Disease.Free.Status)
table(Pheno_tcga$Overall.Survival.Status)

Pheno_tcga$Disease.Free.Status <- as.numeric(Pheno_tcga$Disease.Free.Status)
Pheno_tcga$Overall.Survival.Status <- as.numeric(Pheno_tcga$Overall.Survival.Status)

table(Pheno_tcga$Disease.Free.Status)
table(Pheno_tcga$Overall.Survival.Status)

Pheno_tcga$Disease.Free..Months. <- as.numeric(Pheno_tcga$Disease.Free..Months.)
Pheno_tcga$Overall.Survival..Months. <- as.numeric(Pheno_tcga$Overall.Survival..Months.)

group_tcga <- as.factor(Pheno_tcga$Overall.Survival.Status)
table(group_tcga)

####################################################################################
# save
save(Expr_metabric, Expr_metabric_refAll, Expr_metabric_refDip, Expr_tcga, Expr_tcga_refAll, Expr_tcga_refDip, group_metabric, group_tcga, Pheno_metabric, Pheno_tcga, file = './objs/forKTSP.rda')










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
library(glmnet)
library(RTCGA)
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


Pheno_tcga <- Pheno_tcga[!(Pheno_tcga$Cancer.Type.Detailed %in% c('Adenoid Cystic Breast Cancer', 'Basal Cell Carcinoma', 
                                                                  'Malignant Phyllodes Tumor of the Breast', 
                                                                  'Metaplastic Breast Cancer', 'Paget Disease of the Nipple', 
                                                                  'Solid Papillary Carcinoma of the Breast')), 
]

#############################################################################################
#############################################################################################
## heatmap (CO130)
Expr_tcga_refAll_heatmap <- Expr_tcga_refAll[CO130_fil, ] 

# Create annotation for columns/samples based on some clinical variables:
Pheno_tcga_forHeatmap <- Pheno_tcga
rownames(Pheno_tcga_forHeatmap) <- NULL

AnnAll_tcga <- Pheno_tcga_forHeatmap %>% 
  as.data.frame() %>%
  dplyr::select(Patient.ID, ER.Status.By.IHC, PR.status.by.ihc, HER2.fish.status, IHC.HER2, Subtype) %>%
  #filter(Subtype %in% c('BRCA_Basal', 'BRCA_Her2', 'BRCA_LumA', 'BRCA_LumB')) %>%
  filter(ER.Status.By.IHC %in% c('Negative', 'Positive')) %>%
  filter(PR.status.by.ihc %in% c('Negative', 'Positive')) %>%
  dplyr:: mutate(Subtype = gsub('BRCA_', '', Subtype),
                 Pam50_subtypes = as.factor(Subtype),
                 IHC.HER2 = as.factor(IHC.HER2), 
                 HER2.fish.status = as.factor(HER2.fish.status),
                 ER.Status.By.IHC = as.factor(ER.Status.By.IHC),
                 PR.status.by.ihc = as.factor(PR.status.by.ihc)
  ) %>%
  column_to_rownames(var = "Patient.ID") 


table(AnnAll_tcga$Pam50_subtypes)
table(AnnAll_tcga$ER.Status.By.IHC)
table(AnnAll_tcga$PR.status.by.ihc)


# filter and transpose the expression matrix
Expr_tcga_refAll_heatmap <- Expr_tcga_refAll_heatmap[, rownames(AnnAll_tcga)]
Expr_tcga_refAll_heatmap_t <- t(Expr_tcga_refAll_heatmap)

# filter pheno (above we remove normal and NC)
Pheno_tcga <- Pheno_tcga[rownames(AnnAll_tcga), ]

# colors
ann_colors = list()
ann_colors$Pam50_subtypes <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(6)
names(ann_colors$Pam50_subtypes) <- levels(AnnAll_tcga$Pam50_subtypes)

ann_colors$ER.Status.By.IHC <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(2)
names(ann_colors$ER.Status.By.IHC) <- levels(AnnAll_tcga$ER.Status.By.IHC)

#ann_colors$X3.Gene.classifier.subtype <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(4)
#names(ann_colors$X3.Gene.classifier.subtype) <- levels(AnnAll_tcga$X3.Gene.classifier.subtype)

#ann_colors$Neoplasm.Histologic.Grade <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(3)
#names(ann_colors$Neoplasm.Histologic.Grade) <- levels(AnnAll_tcga$Neoplasm.Histologic.Grade)


breaksList = seq(-4, 4, by = 1)
ColPal <- colorRampPalette(colors = rev(brewer.pal(11,"RdYlBu")))(20)
ColPal2 <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(20))


#######################################################
# get the 5 groups
heat_tcga <- pheatmap(Expr_tcga_refAll_heatmap, 
                      scale = "none",
                      #color = rev(heat.colors(20)),
                      color =ColPal,
                      annotation_colors = ann_colors,
                      cluster_cols = T, 
                      cluster_rows = T, 
                      clustering_distance_cols = 'correlation',
                      clustering_distance_rows = 'correlation',
                      clustering_method = 'ward.D',
                      show_colnames = F,
                      show_rownames = T,
                      annotation_col = AnnAll_tcga,
                      annotation_names_col = T,
                      #annotation_row = AnnAll_tcga,
                      annotation_names_row = T,
                      fontsize = 7,
                      #fontsize_col = 3,
                      fontsize_row = 10,
                      silent = TRUE,
                      cex = 1,
                      cutree_cols = 5,
                      cutree_rows = 5,
                      breaks = seq(-1, 1, by = 0.1),
                      main = "")

clusters_tcga <- as.data.frame(cbind(t(Expr_tcga_refAll_heatmap), 
                                     'THR clusters' = cutree(heat_tcga$tree_col, 
                                                             k = 5)))

table(clusters_tcga$`THR clusters`)

# add the cluster info to the phenotype table
all(rownames(clusters_tcga) == rownames(Pheno_tcga))
Pheno_tcga$`THR clusters` <- clusters_tcga$`THR clusters`

# add the cluster info to the Ann dataframe and re-plot the heatmap
all(rownames(clusters_tcga) == rownames(AnnAll_tcga))
AnnAll_tcga$`THR clusters` <- as.factor(paste0('c', clusters_tcga$`THR clusters`))
table(AnnAll_tcga$`THR clusters`)

# re-order the annotation dataframe then the expression matrix by cluster
#AnnAll_tcga <- AnnAll_tcga[order(AnnAll_tcga$cluster, decreasing = FALSE), ]
#Expr_tcga_refAll_heatmap <- Expr_tcga_refAll_heatmap[, rownames(AnnAll_tcga)]


ann_colors$`THR clusters` <- colorRampPalette(colors = rev(brewer.pal(5,"Dark2")))(5)
names(ann_colors$`THR clusters`) <- levels(AnnAll_tcga$`THR clusters`)

# heatmap with cluster annotation
tiff('./figures/CO130_tcga_clusters/CO130_heatmap_tcga_clusters_refAll.tiff', width=3000, height=2000, res = 300)
pheatmap(Expr_tcga_refAll_heatmap, 
         scale = "none",
         #color = rev(heat.colors(20)),
         color = ColPal,
         annotation_colors = ann_colors,
         cluster_cols = T, 
         cluster_rows = T, 
         clustering_distance_cols = 'correlation',
         clustering_distance_rows = 'correlation',
         clustering_method = 'ward.D',
         show_colnames = F,
         show_rownames = T,
         annotation_col = AnnAll_tcga,
         annotation_names_col = T,
         #annotation_row = AnnAll_metabric,
         annotation_names_row = T,
         fontsize = 7,
         #fontsize_col = 3,
         fontsize_row = 3,
         cex = 1,
         cutree_cols = 5,
         cutree_rows = 5,
         breaks = seq(-1, 1, by = 0.1),
         main = "")
dev.off()

# another one with T and E annotation
levels(AnnAll_tcga$`THR clusters`) <- c('E1', 'E2', 'E4', 'E3', 'T1')
#ann_colors$`THR clusters` <- c('#66a61e', "#E7298A", '#D95F02', "#1B9E77", '#7570b3')

names(ann_colors$`THR clusters`) <- levels(AnnAll_tcga$`THR clusters`)

tiff('./figures/CO130_tcga_clusters/CO130_heatmap_tcga_clusters_E_refAll.tiff', width=3000, height=2000, res = 300)
pheatmap(Expr_tcga_refAll_heatmap, 
         scale = "none",
         #color = rev(heat.colors(20)),
         color =ColPal,
         annotation_colors = ann_colors,
         cluster_cols = T, 
         cluster_rows = T, 
         clustering_distance_cols = 'correlation',
         clustering_distance_rows = 'correlation',
         clustering_method = 'ward.D',
         show_colnames = F,
         show_rownames = T,
         annotation_col = AnnAll_tcga,
         annotation_names_col = T,
         #annotation_row = AnnAll_metabric,
         annotation_names_row = T,
         fontsize = 7,
         #fontsize_col = 3,
         fontsize_row = 3,
         cex = 1,
         cutree_cols = 5,
         cutree_rows = 5,
         breaks = seq(-1, 1, by = 0.1),
         main = "")
dev.off()

# fix the cluster names in the pheno table
table(Pheno_tcga$`THR clusters`)
Pheno_tcga$`THR clusters` <- as.factor(Pheno_tcga$`THR clusters`)
levels(Pheno_tcga$`THR clusters`) <- c('E1', 'E2', 'E4', 'E3', 'T1')
#AnnAll_tcga$`THR clusters` <- as.factor(paste0('c', clusters_tcga$`THR clusters`))
table(Pheno_tcga$`THR clusters`)
table(AnnAll_tcga$`THR clusters`)

# fix the pam50 subtypes
Pheno_tcga$Pam50_subtypes <- as.factor(Pheno_tcga$Subtype)
Pheno_tcga$Pam50_subtypes <- gsub('BRCA_', '', Pheno_tcga$Pam50_subtypes)
table(Pheno_tcga$Pam50_subtypes)

#############################################################################################################
##############################################################################################################

## Keep only the relevant information (Metastasis Event and Time)
survival_tcga <- Pheno_tcga[, c("Overall.Survival.Status", "Overall.Survival..Months.", 
                                "Disease.Free.Status", "Disease.Free..Months.", 
                                "Pam50_subtypes", "ER.Status.By.IHC",
                                "PR.status.by.ihc", "HER2.fish.status", "IHC.HER2", "THR clusters")] 

table(survival_tcga$`THR clusters`)
survival_tcga$`THR clusters` <- as.factor(survival_tcga$`THR clusters`)

survival_tcga$Disease.Free.Status <- gsub("\\:.+", "", survival_tcga$Disease.Free.Status)
survival_tcga$Overall.Survival.Status <- gsub("\\:.+", "", survival_tcga$Overall.Survival.Status)

survival_tcga$Disease.Free.Status <- as.numeric(survival_tcga$Disease.Free.Status)
survival_tcga$Overall.Survival.Status <- as.numeric(survival_tcga$Overall.Survival.Status)

table(survival_tcga$Disease.Free.Status)
table(survival_tcga$Overall.Survival.Status)

table(Pheno_tcga$Disease.Free.Status)
table(Pheno_tcga$Overall.Survival.Status)



survival_tcga$Disease.Free..Months. <- as.numeric(survival_tcga$Disease.Free..Months.)
survival_tcga$Overall.Survival..Months. <- as.numeric(survival_tcga$Overall.Survival..Months.)


summary(survival_tcga$Overall.Survival..Months.)
summary(survival_tcga$Disease.Free..Months.)

summary(as.numeric(Pheno_tcga$Overall.Survival..Months.))
summary(as.numeric(Pheno_tcga$Disease.Free..Months.))


# OS
Fit_tcga_os <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ as.factor(`THR clusters`), data = survival_tcga)

# DFS
Fit_tcga_DFS <- survfit(Surv(Disease.Free..Months., Disease.Free.Status) ~ as.factor(`THR clusters`), data = survival_tcga)

############################################################################
############################################################################
# plot OS
cluster_colors <- as.vector(ann_colors$`THR clusters`)

pdf("./figures/CO130_tcga_clusters/tcga_os_5clusters_refAll.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_tcga_os,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = levels(survival_tcga$`THR clusters`),
           legend.title	= '',
           pval.size = 12,
           break.x.by = 40,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR50 clusters and OS in TCGA'
           )
dev.off()


png("./figures/CO130_tcga_clusters/CO130_dfs_5clusters_refAll.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_tcga_DFS,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = levels(survival_tcga$`THR clusters`),
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
) + guides(
  colour = guide_legend(ncol = 3))
dev.off()



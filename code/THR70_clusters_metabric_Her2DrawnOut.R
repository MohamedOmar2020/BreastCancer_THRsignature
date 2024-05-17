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
library(readxl)
library(GSVA)
#################
# THR signature --------------
#################
THR_signature <- readxl::read_xlsx("./data/THR_Signatures_Jan25_2023.xlsx")

# get the THR70 signature
THR70 <- THR_signature$`THR-70`[!is.na(THR_signature$`THR-70`)]

THR70 <- gsub('-', '', THR70)

#############################################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')
##############################################

##############################################
# fix gene names
##############################################
# Fix in TCGA
setdiff(THR70, rownames(Expr_tcga_refAll))
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
setdiff(THR70, rownames(Expr_metabric_refAll))
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
THR70_fil <- THR70[THR70 %in% rownames(Expr_tcga_refAll) & THR70 %in% rownames(Expr_metabric_refAll)]

setdiff(THR70, THR70_fil)

#############################################################################################
#############################################################################################
## heatmap
Expr_metabric_refAll_heatmap <- Expr_metabric_refAll[THR70_fil, ] 

# Create annotation for columns/samples based on some clinical variables:
Pheno_metabric_forHeatmap <- Pheno_metabric
rownames(Pheno_metabric_forHeatmap) <- NULL

AnnAll_metabric <- Pheno_metabric_forHeatmap %>% 
  as.data.frame() %>%
  dplyr::select(Sample.ID, Pam50...Claudin.low.subtype, X3.Gene.classifier.subtype, HER2.Status, PR.Status, ER.status.measured.by.IHC, Neoplasm.Histologic.Grade) %>%
  column_to_rownames(var = "Sample.ID") %>%
  filter(Pam50...Claudin.low.subtype %in% c('Basal', 'claudin-low', 'Her2', 'LumA', 'LumB')) %>%
  dplyr::mutate(X3.Gene.classifier.subtype = as.factor(X3.Gene.classifier.subtype),
                ER.status.measured.by.IHC = as.factor(ER.status.measured.by.IHC),
                Pam50...Claudin.low.subtype = as.factor(Pam50...Claudin.low.subtype),
                HER2.Status = as.factor(HER2.Status), 
                PR.Status = as.factor(PR.Status),
                #Overall.Survival.Status = as.factor(Overall.Survival.Status),
                #Relapse.Free.Status = as.factor(Relapse.Free.Status), 
                #Tumor.Stage = as.factor(Tumor.Stage), 
                Neoplasm.Histologic.Grade = as.factor(Neoplasm.Histologic.Grade))


# filter and transpose the expression matrix
Expr_metabric_refAll_heatmap <- Expr_metabric_refAll_heatmap[, rownames(AnnAll_metabric)]
Expr_metabric_refAll_heatmap_t <- t(Expr_metabric_refAll_heatmap)

# filter pheno (above we remove normal and NC)
Pheno_metabric <- Pheno_metabric[rownames(AnnAll_metabric), ]

# colors
ann_colors = list()
ann_colors$Pam50...Claudin.low.subtype <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(5)
names(ann_colors$Pam50...Claudin.low.subtype) <- levels(AnnAll_metabric$Pam50...Claudin.low.subtype)

ann_colors$ER.status.measured.by.IHC <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(2)
names(ann_colors$ER.status.measured.by.IHC) <- levels(AnnAll_metabric$ER.status.measured.by.IHC)

ann_colors$X3.Gene.classifier.subtype <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(4)
names(ann_colors$X3.Gene.classifier.subtype) <- levels(AnnAll_metabric$X3.Gene.classifier.subtype)

ann_colors$Neoplasm.Histologic.Grade <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(3)
names(ann_colors$Neoplasm.Histologic.Grade) <- levels(AnnAll_metabric$Neoplasm.Histologic.Grade)


breaksList = seq(-4, 4, by = 1)
ColPal <- colorRampPalette(colors = rev(brewer.pal(11,"RdYlBu")))(20)
ColPal2 <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(20))

#######################################################
# get the 5 groups
heat_metabric <- pheatmap(Expr_metabric_refAll_heatmap, 
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
                          annotation_col = AnnAll_metabric,
                          annotation_names_col = T,
                          #annotation_row = AnnAll_metabric,
                          annotation_names_row = T,
                          fontsize = 7,
                          #fontsize_col = 3,
                          fontsize_row = 10,
                          cex = 1,
                          cutree_cols = 5,
                          cutree_rows = 5,
                          breaks = seq(-1, 1, by = 0.1),
                          silent = TRUE,
                          main = "")

clusters_metabric <- as.data.frame(cbind(t(Expr_metabric_refAll_heatmap), 
                                         'THR clusters' = cutree(heat_metabric$tree_col, 
                                                                 k = 5)))
table(clusters_metabric$`THR clusters`)


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# Modify the clusters by Extracting HER2+ samples----------------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
table(Pheno_metabric$X3.Gene.classifier.subtype)

HER2_samples <- which(Pheno_metabric$X3.Gene.classifier.subtype == "HER2+")
clusters_metabric$`THR clusters`[HER2_samples] <- max(clusters_metabric$`THR clusters`) + 1  # Add 1 to the maximum cluster number to ensure a new distinct cluster

table(clusters_metabric$`THR clusters`)

#######
# rename the levels
clusters_metabric$`THR clusters` <- as.factor(clusters_metabric$`THR clusters`)
levels(clusters_metabric$`THR clusters`) <- c('E2', 'E2', 'E1', 'E3', 'PQNBC', 'HER2+')

table(clusters_metabric$`THR clusters`)

##############################
## add the cluster info to the phenotype table

# Reorder the rows of Pheno_metabric
Pheno_metabric <- Pheno_metabric[match(rownames(clusters_metabric), rownames(Pheno_metabric)), ]
all(rownames(clusters_metabric) == rownames(Pheno_metabric))

# Add the 'THR clusters' column
Pheno_metabric$`THR clusters` <- clusters_metabric$`THR clusters`
table(Pheno_metabric$`THR clusters`)

#################
# add the cluster info to the Ann dataframe and re-plot the heatmap
AnnAll_metabric <- AnnAll_metabric[match(rownames(clusters_metabric), rownames(AnnAll_metabric)), ]
all(rownames(clusters_metabric) == rownames(AnnAll_metabric))
AnnAll_metabric$`THR clusters` <- as.factor(clusters_metabric$`THR clusters`)
table(AnnAll_metabric$`THR clusters`)

ann_colors <- list()
ann_colors$Pam50...Claudin.low.subtype <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(5)
names(ann_colors$Pam50...Claudin.low.subtype) <- levels(AnnAll_metabric$Pam50...Claudin.low.subtype)

ann_colors$ER.status.measured.by.IHC <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(2)
names(ann_colors$ER.status.measured.by.IHC) <- levels(AnnAll_metabric$ER.status.measured.by.IHC)

ann_colors$X3.Gene.classifier.subtype <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(4)
names(ann_colors$X3.Gene.classifier.subtype) <- levels(AnnAll_metabric$X3.Gene.classifier.subtype)

ann_colors$Neoplasm.Histologic.Grade <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(3)
names(ann_colors$Neoplasm.Histologic.Grade) <- levels(AnnAll_metabric$Neoplasm.Histologic.Grade)

ann_colors$`THR clusters` <- colorRampPalette(colors = rev(brewer.pal(5,"Dark2")))(5)
#levels(AnnAll_metabric$`THR clusters`) <- c('E1', 'E2', 'E3', 'T1', 'T1', 'HER2+')
names(ann_colors$`THR clusters`) <- levels(AnnAll_metabric$`THR clusters`)
table(AnnAll_metabric$`THR clusters`)

breaksList = seq(-4, 4, by = 1)
ColPal <- colorRampPalette(colors = rev(brewer.pal(11,"RdYlBu")))(20)
ColPal2 <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(20))

#####################################################
## Fix the order of samples in Expr_metabric_refAll_heatmap
Expr_metabric_refAll_heatmap <- Expr_metabric_refAll_heatmap[, match(rownames(clusters_metabric), colnames(Expr_metabric_refAll_heatmap))]
all(rownames(clusters_metabric) == colnames(Expr_metabric_refAll_heatmap))

###################################
# Create an ordering based on the `THR clusters` column in `AnnAll_metabric`
order_indices <- order(AnnAll_metabric$`THR clusters`)

# Use the order to rearrange both `Expr_metabric_refAll_heatmap` and `AnnAll_metabric`
Expr_metabric_refAll_heatmap_ordered <- Expr_metabric_refAll_heatmap[, order_indices]
AnnAll_metabric_ordered <- AnnAll_metabric[order_indices, ]

table(AnnAll_metabric_ordered$`THR clusters`)
##############
cluster_ends <- cumsum(table(AnnAll_metabric_ordered$`THR clusters`))

#########################################
# heatmap
tiff('./figures/THR70_metabric_clusters/THR70_heatmap_metabric_clusters_HER2_DrawnOut.tiff', width=3000, height=2000, res = 300)
pheatmap(Expr_metabric_refAll_heatmap_ordered, 
         scale = "none",
         #color = rev(heat.colors(20)),
         color =ColPal,
         annotation_colors = ann_colors,
         cluster_cols = F, 
         cluster_rows = T, 
         clustering_distance_cols = 'correlation',
         clustering_distance_rows = 'correlation',
         clustering_method = 'ward.D',
         show_colnames = F,
         show_rownames = T,
         annotation_col = AnnAll_metabric_ordered,
         annotation_names_col = T,
         #annotation_row = AnnAll_metabric,
         annotation_names_row = T,
         fontsize = 7,
         #fontsize_col = 3,
         fontsize_row = 3,
         cex = 1,
         cutree_cols = 5,
         cutree_rows = 5,
         gaps_col = cluster_ends,  # Add gaps
         breaks = seq(-1, 1, by = 0.1),
         main = "")
dev.off()

#############################################################################################################
##############################################################################################################

## Keep only the relevant information (Metastasis Event and Time)
survival_metabric <- Pheno_metabric[, c("Overall.Survival.Status", "Overall.Survival..Months.", 
                                        "Relapse.Free.Status", "Relapse.Free.Status..Months.", 
                                        "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC",
                                        "X3.Gene.classifier.subtype", "THR clusters")] 

survival_metabric$`THR clusters` <- as.factor(survival_metabric$`THR clusters`)
#levels(survival_metabric$`THR clusters`) <- c('E1', 'E2', 'E3', 'T1', 'T1', 'HER2+')
#survival_metabric$`THR clusters` <- factor(survival_metabric$`THR clusters`, levels = c('E1', 'E2a', 'E2b', 'E3', 'PNBC'))
# OS
Fit_metabric_os <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ as.factor(`THR clusters`), data = survival_metabric)

# RFS
Fit_metabric_RFS <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ as.factor(`THR clusters`), data = survival_metabric)


############################################################################
############################################################################
# plot OS
cluster_colors <- as.vector(ann_colors$`THR clusters`)
#cluster_colors <- levels(survival_metabric$`THR clusters`)

#cluster_colors <- c("#1B9E77", "#E7298A", "#66A61E" , "#D95F02", "#7570B3") 
#cluster_colors2 <- c('#1B9E77', "#c1b026", "#D95F02", '#7570B3')

#names(ann_colors$`THR clusters`) <- levels(survival_metabric$`THR clusters`)
#cluster_colors <- as.vector(ann_colors$`THR clusters`)
png("./figures/THR70_metabric_clusters/THR70_metabric_os_5clusters_20yrs_Her2_DrawnOut_2.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_metabric_os,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = levels(survival_metabric$`THR clusters`),
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

# plot RFS
png("./figures/THR70_metabric_clusters/THR70_metabric_rfs_5clusters_20yrs_Her2_DrawnOut_2.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_metabric_RFS,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = levels(survival_metabric$`THR clusters`),
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

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Divide T1 using immune signatures ------------------------
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## get T1----------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
T1_pheno <- Pheno_metabric[Pheno_metabric$`THR clusters` == 'PQNBC', ]
T1_expr <- Expr_metabric_refAll[, rownames(T1_pheno)]

all(rownames(T1_pheno) == colnames(T1_expr))

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## Load the i20 signature ---------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
i20 <- read_xlsx("./figures/c3_DE_THR50_RFS/THR50_c3_longVSshortSurv_DE.xlsx")$gene

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## Load the i45 signature ---------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
load('objs/THR56T1sep_i_genesets.rda')

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## use I20 to separate T1 into two clusters
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#&&&&&&&&&&&&&&&&&&
### calculate signature score---------
#&&&&&&&&&&&&&&&&&&
i20_filt <- intersect(i20, rownames(Expr_metabric_refAll))
T1_i20_expr <- T1_expr[i20_filt, ]

i45_filt <- intersect(i30p15, rownames(Expr_metabric_refAll))
T1_i45_expr <- T1_expr[i45_filt, ]

# using the average expression
#T1_signature_scores <- colMeans(T1_i20_expr, na.rm = TRUE)  # Use colMeans since samples are in columns

# using GSVA
gene_sets <- list(i20 = rownames(T1_i20_expr), i45 = rownames(T1_i45_expr))
T1_signature_scores <- gsva(T1_expr, gene_sets, method = "ssgsea")
T1_signature_scores_i20 <- as.vector(T1_signature_scores[1, ])
T1_signature_scores_i45 <- as.vector(T1_signature_scores[2, ])

T1_pheno$i20_score <- T1_signature_scores_i20
T1_pheno$i45_score <- T1_signature_scores_i45

# Perform PCA on the expression data of the i20 genes
#pca_result <- prcomp(t(T1_i20_expr), scale. = TRUE)
# Use the first principal component as the signature score
#T1_signature_scores <- pca_result$x[, 1]

# using z-score
#T1_signature_scores <- colMeans(scale(T1_i20_expr), na.rm = TRUE)

# using logistic regression:
Data_T1 <- data.frame(cbind(t(T1_expr), 'RFS' = T1_pheno$Relapse.Free.Status))
Data_T1$RFS <- as.factor(Data_T1$RFS)
#levels(Data_T1$RFS) <- c('longSurv', 'shortSurv')
table(Data_T1$RFS)
i20_logReg_model <- glm(as.formula((paste("RFS ~", paste(i20_filt, collapse = "+")))), data = Data_T1, family = "binomial")
T1_pheno$i20_logReg_score <- i20_logReg_model %>% predict(Data_T1 , type = "response")

#&&&&&&&&&&&&&&&&&&
### determine the best threshold-------
#&&&&&&&&&&&&&&&&&&
median_threshold_i20 <- median(T1_signature_scores_i20)
median_threshold_i45 <- median(T1_signature_scores_i45)

quantile_threshold_i20 <- quantile(T1_signature_scores_i20, probs = 0.75)  # 75th percentile as an example
quantile_threshold_i45 <- quantile(T1_signature_scores_i45, probs = 0.75)  # 75th percentile as an example

ROC_thr_I20 <- coords(roc(T1_pheno$Relapse.Free.Status, T1_pheno$i20_logReg_score, direction = "<"), "best")["threshold"]

#&&&&&&&&&&&&&&&&&&&&&&
# save the logistic regression model and the threshold-------
#&&&&&&&&&&&&&&&&&&&&&&
save(i20_logReg_model, ROC_thr_I20, file = './objs/metabric_THR70_derived_T1_THR50_i20_logReg_model_Her2_DrawnOut.rda')


#&&&&&&&&&&&&&&&&&&
### Assign samples to subclass based on i20 expression-------
#&&&&&&&&&&&&&&&&&&
T1_pheno$immune_clusters_i20 <- ifelse(T1_signature_scores_i20 >= median_threshold_i20, "PQNBC_i+", "PQNBC_i-")
table(T1_pheno$immune_clusters_i20)
T1_pheno$immune_clusters_i20 <- factor(T1_pheno$immune_clusters_i20, levels = c('PQNBC_i-', 'PQNBC_i+'))

T1_pheno$immune_clusters_i45 <- ifelse(T1_signature_scores_i45 >= median_threshold_i45, "PQNBC_i+", "PQNBC_i-")
table(T1_pheno$immune_clusters_i45)
T1_pheno$immune_clusters_i45 <- factor(T1_pheno$immune_clusters_i45, levels = c('PQNBC_i-', 'PQNBC_i+'))

T1_pheno$immune_clusters_i20_logReg <- ifelse(T1_pheno$i20_logReg_score >= ROC_thr_I20$threshold, "PQNBC_i-", "PQNBC_i+")
table(T1_pheno$immune_clusters_i20_logReg)
T1_pheno$immune_clusters_i20_logReg <- factor(T1_pheno$immune_clusters_i20_logReg, levels = c('PQNBC_i-', 'PQNBC_i+'))


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## Keep only the relevant information (Metastasis Event and Time)-----
T1_pheno_survival <- T1_pheno[, c("Overall.Survival.Status", "Overall.Survival..Months.", 
                                  "Relapse.Free.Status", "Relapse.Free.Status..Months.", 
                                  "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC",
                                  "X3.Gene.classifier.subtype", 
                                  "immune_clusters_i20_logReg"
                                  )]

#&&&&&&&&&&&&&&&&&&&&&&&&&&
CoxData_metabric_T1 <- data.frame(T1_pheno_survival)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## survival analysis for just T1-----
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

### OS------
Fit_sig_metabric_os_T1_i20 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ immune_clusters_i20, data = CoxData_metabric_T1)
Fit_sig_metabric_os_T1_i45 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ immune_clusters_i45, data = CoxData_metabric_T1)
Fit_sig_metabric_os_T1_i20_logReg <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ immune_clusters_i20_logReg, data = CoxData_metabric_T1)


### RFS-----
Fit_sig_metabric_RFS_T1_i20 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ immune_clusters_i20, data = CoxData_metabric_T1)
Fit_sig_metabric_RFS_T1_i45 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ immune_clusters_i45, data = CoxData_metabric_T1)
Fit_sig_metabric_RFS_T1_i20_logReg <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ immune_clusters_i20_logReg, data = CoxData_metabric_T1)


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## recombine T1 with the rest---------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

T1 <- data.frame(T1_immune_clusters_i20_logReg = T1_pheno$immune_clusters_i20_logReg,
                 `Sample.ID` = rownames(T1_pheno))

rownames(T1) <- rownames(T1_pheno)


# merge
Pheno_metabric$`THR clusters`[Pheno_metabric$`THR clusters` == 'PQNBC'] <- NA


Pheno_metabric2 <- merge(x = T1, y = Pheno_metabric, by = "Sample.ID", all.y = TRUE)

Pheno_metabric2 <- Pheno_metabric2 %>% 
  mutate(merged_THR_clusters_i20_logReg = coalesce(T1_immune_clusters_i20_logReg,`THR clusters`))

Pheno_metabric2$merged_THR_clusters_i20 <- droplevels(Pheno_metabric2$merged_THR_clusters_i20)
Pheno_metabric2$merged_THR_clusters_i20_logReg <- droplevels(Pheno_metabric2$merged_THR_clusters_i20_logReg)
Pheno_metabric2$merged_THR_clusters_i45 <- droplevels(Pheno_metabric2$merged_THR_clusters_i45)

table(Pheno_metabric2$merged_THR_clusters_i20)
table(Pheno_metabric2$merged_THR_clusters_i20_logReg)
table(Pheno_metabric2$merged_THR_clusters_i45)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## survival analysis---------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## Keep only the relevant information (Metastasis Event and Time)
survival_metabric <- Pheno_metabric2[, c("Overall.Survival.Status", "Overall.Survival..Months.", 
                                         "Relapse.Free.Status", "Relapse.Free.Status..Months.", 
                                         "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC",
                                         "X3.Gene.classifier.subtype", 
                                         "merged_THR_clusters_i20_logReg" 
)] 


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
### OS---------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# fix the levels
#levels(survival_metabric$merged_THR_clusters_i20) <- c("PQNBC_i-", "PQNBC_i+", "E2", 'E3', 'E1', 'HER2+')

cluster_colors <- c("#1960b2", "#DF271F", "#66A61E", "#CD4174", "#A253A2", "#9D696C", "#1B9E77")

Fit_metabric_os_i20 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ merged_THR_clusters_i20, data = survival_metabric)

png("./figures/THR70_metabric_clusters/THR70_metabric_os_5clusters_HER2_DrawnOut_10yrs_T1_i20_2.png", width = 2000, height = 2000, res = 300)
ggsurvplot(Fit_metabric_os_i20,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = levels(survival_metabric$merged_THR_clusters_i20),
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

# Cox
cox_metabric <- survival_metabric
cox_metabric$merged_THR_clusters_i20 <- factor(cox_metabric$merged_THR_clusters_i20, levels = c('E3', 'E2', 'E1', 'HER2+', 'PQNBC_i+', 'PQNBC_i-')) 
cox_Fit_metabric_os_THR70_with_i20 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ merged_THR_clusters_i20, data = cox_metabric)
summary(cox_Fit_metabric_os_THR70_with_i20)

##################
Fit_metabric_os_i20_logReg <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ merged_THR_clusters_i20_logReg, data = survival_metabric)

png("./figures/THR70_metabric_clusters/THR70_metabric_os_5clusters_HER2_DrawnOut_10yrs_T1_i20_logReg.png", width = 2000, height = 2000, res = 300)
ggsurvplot(Fit_metabric_os_i20_logReg,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = levels(survival_metabric$merged_THR_clusters_i20_logReg),
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

# Cox
cox_metabric$merged_THR_clusters_i20_logReg <- factor(cox_metabric$merged_THR_clusters_i20_logReg, levels = c('E3', 'E2', 'E1', 'HER2+', 'PQNBC_i+', 'PQNBC_i-')) 
cox_Fit_metabric_os_THR70_with_i20_logReg <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ merged_THR_clusters_i20_logReg, data = cox_metabric)
summary(cox_Fit_metabric_os_THR70_with_i20_logReg)



#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
### RFS---------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Fit_metabric_rfs_i20 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ merged_THR_clusters_i20, data = survival_metabric)

png("./figures/THR70_metabric_clusters/THR70_metabric_rfs_5clusters_HER2_DrawnOut_10yrs_T1_i20_2.png", width = 2000, height = 2000, res = 300)
ggsurvplot(Fit_metabric_rfs_i20,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = levels(survival_metabric$merged_THR_clusters_i20),
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

# Cox
cox_Fit_metabric_rfs_THR70_with_i20 <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ merged_THR_clusters_i20, data = cox_metabric)
summary(cox_Fit_metabric_rfs_THR70_with_i20)

#&&&&&&&&&&&&
# i20_logReg
Fit_metabric_rfs_i20_logReg <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ merged_THR_clusters_i20_logReg, data = survival_metabric)
png("./figures/THR70_metabric_clusters/THR70_metabric_rfs_5clusters_HER2_DrawnOut_10yrs_T1_i20_logReg.png", width = 2000, height = 2000, res = 300)
ggsurvplot(Fit_metabric_rfs_i20_logReg,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = c('PNBC.i-', 'PNBC.i+', 'E2', 'E1', 'E3', 'HER2+'),
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


# Cox
cox_metabric <- survival_metabric
cox_metabric$merged_THR_clusters_i20_logReg <- factor(cox_metabric$merged_THR_clusters_i20_logReg, levels = c('PQNBC_i+','PQNBC_i-', 'E3', 'E2', 'E1', 'HER2+')) 

cox_Fit_metabric_rfs_THR70_with_i20_logReg <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ merged_THR_clusters_i20_logReg, data = cox_metabric)
summary(cox_Fit_metabric_rfs_THR70_with_i20_logReg)


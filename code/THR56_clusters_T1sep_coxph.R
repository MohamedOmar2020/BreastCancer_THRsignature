

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

#################
THR_signature <- readxl::read_xlsx("./data/THR Signatures_sep23.xlsx")

# get the THR50 signature
THR_50 <- THR_signature$`THR-50.1`[!is.na(THR_signature$`THR-50.1`)]

THR_50 <- gsub('-', '', THR_50)

################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')

# fix gene names
rownames(Expr_metabric_refAll)[grep('^ZNF652', rownames(Expr_metabric_refAll))]

# filter the THR signatures to include only the genes present in the expr matrices
THR_50_fil <- THR_50[THR_50 %in% rownames(Expr_metabric_refAll)]

THR_56 <- c(THR_50_fil, 'LMNB2', 'CDC20', 'KIF2C', 'FAM64A', 'KIF4A', 'TPX2')

#############################################################################################
#############################################################################################
## heatmap (THR 50)
Expr_metabric_refAll_heatmap <- Expr_metabric_refAll[THR_56, ] 

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
                          silent = TRUE,
                          cex = 1,
                          cutree_cols = 5,
                          cutree_rows = 5,
                          breaks = seq(-1, 1, by = 0.1),
                          main = "")

clusters_metabric <- as.data.frame(cbind(t(Expr_metabric_refAll_heatmap), 
                                         'THR clusters' = cutree(heat_metabric$tree_col, 
                                                                 k = 5)))

table(clusters_metabric$`THR clusters`)

# add the cluster info to the phenotype table
all(rownames(clusters_metabric) == rownames(Pheno_metabric))
Pheno_metabric$`THR clusters` <- clusters_metabric$`THR clusters`

# add the cluster info to the Ann dataframe and re-plot the heatmap
all(rownames(clusters_metabric) == rownames(AnnAll_metabric))
AnnAll_metabric$`THR clusters` <- as.factor(paste0('c', clusters_metabric$`THR clusters`))
table(AnnAll_metabric$`THR clusters`)

# re-order the annotation dataframe then the expression matrix by cluster
#AnnAll_metabric <- AnnAll_metabric[order(AnnAll_metabric$cluster, decreasing = FALSE), ]
#Expr_metabric_refAll_heatmap <- Expr_metabric_refAll_heatmap[, rownames(AnnAll_metabric)]


ann_colors$`THR clusters` <- colorRampPalette(colors = rev(brewer.pal(5,"Dark2")))(5)
ann_colors$`THR clusters` <- c("#66A61E", "#7570B3", "#1B9E77" , "#D95F02", "#E7298A") 
levels(AnnAll_metabric$`THR clusters`) <- c('E1', 'T1', 'E3', 'E2', 'E2')
names(ann_colors$`THR clusters`) <- levels(AnnAll_metabric$`THR clusters`)

# fix the cluster names in the pheno table
table(Pheno_metabric$`THR clusters`)
Pheno_metabric$`THR clusters` <- as.factor(Pheno_metabric$`THR clusters`)
levels(Pheno_metabric$`THR clusters`) <- c('E1', 'T1', 'E3', 'E2', 'E2')
table(Pheno_metabric$`THR clusters`)


#############################################################################################################
# get cluster 3 with crossing curves
T1_pheno <- Pheno_metabric[Pheno_metabric$`THR clusters` == 'T1', ]
T1_expr <- Expr_metabric_refAll[, rownames(T1_pheno)]

all(rownames(T1_pheno) == colnames(T1_expr))

######################################
# load the THR50-derived i20 genes
THR50_i20 <- readxl::read_xlsx("./figures/c3_DE_THR50_RFS/THR50_c3_longVSshortSurv_DE.xlsx")$gene

# Subset expression data to include only THR50_i20 genes
T1_expr_i20 <- as.data.frame(t(T1_expr[THR50_i20, ]))

# Prepare survival data
surv_data <- Surv(T1_pheno$Relapse.Free.Status..Months., T1_pheno$Relapse.Free.Status)

# Fit Cox model using only THR50_i20 genes
cox_model <- coxph(surv_data ~ ., data = T1_expr_i20)

# Get summary
summary(cox_model)

# Predict risk scores
risk_scores <- predict(cox_model, newdata = T1_expr_i20)

cutoff_25th <- quantile(risk_scores, 0.5)

# Assign samples to groups based on risk scores
T1_pheno$group <- ifelse(risk_scores <= cutoff_25th, "T1a", "T1b")

# Estimate survival curves
surv_fit <- survfit(surv_data ~ group, data = T1_pheno)

# Plot survival curves
ggsurvplot(surv_fit)

###########################################################################
T1 <- data.frame(THR_clusters_i20model = T1_pheno$group,
                 `Sample.ID` = rownames(T1_pheno))

rownames(T1) <- rownames(T1_pheno)


# merge
Pheno_metabric$`THR clusters`[Pheno_metabric$`THR clusters` == 'T1'] <- NA


Pheno_metabric2 <- merge(x = T1, y = Pheno_metabric, by="Sample.ID", all.y = TRUE)

Pheno_metabric2 <- Pheno_metabric2 %>% 
  mutate(`THR clusters` = as.factor(`THR clusters`), THR_clusters_i20model = as.factor(THR_clusters_i20model)) %>%
  mutate(THR.clusters_i20model_Merged = coalesce(THR_clusters_i20model,`THR clusters`))

###########################################################################################
##########################################################################################
## survival analysis

## Keep only the relevant information (Metastasis Event and Time)
survival_metabric <- Pheno_metabric2[, c("Overall.Survival.Status", "Overall.Survival..Months.", 
                                         "Relapse.Free.Status", "Relapse.Free.Status..Months.", 
                                         "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC",
                                         "X3.Gene.classifier.subtype", 
                                         "THR.clusters_i20model_Merged"
)] 

survival_metabric$THR.clusters_i20model_Merged <- as.factor(survival_metabric$THR.clusters_i20model_Merged)
survival_metabric$THR.clusters_i20model_Merged <- droplevels(survival_metabric$THR.clusters_i20model_Merged)


cluster_colors <- as.vector(ann_colors$THR.clusters_i20model_Merged) # same order for the others

################################################################
## Survival curves: model 20
################################################################

# OS
Fit_metabric_os_i20model <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_i20model_Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_i20model <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_i20model_Merged, data = survival_metabric)

pdf("./figures/logreg/THR56_clusters/metabric_os_5clusters_model20Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_i20model,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_i20model_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and OS: THR56 + THR50-derived i20')
dev.off()

## RFS: 
pdf("./figures/logreg/THR56_clusters/metabric_rfs_5clusters_model20Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_i20model,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$THR.clusters_i20model_Merged),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR56 clusters and RFS: THR56 + THR50-derived i20')
dev.off()

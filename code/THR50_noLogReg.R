
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
library(dendextend)

#################
THR_signature <- readxl::read_xlsx("./data/THR Signatures_sep23.xlsx")

# get the THR50 signature
THR_50 <- THR_signature$`THR-50.1`[!is.na(THR_signature$`THR-50.1`)]

THR_50 <- gsub('-', '', THR_50)

################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')

# fix gene names
rownames(Expr_metabric)[grep('^ZNF652', rownames(Expr_metabric))]

# filter the THR signatures to include only the genes present in the expr matrices
THR_50_fil <- THR_50[THR_50 %in% rownames(Expr_metabric)]

#############################################################################################
#############################################################################################
## heatmap (THR 50)
Expr_metabric_heatmap <- Expr_metabric[THR_50_fil, ] 

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
Expr_metabric_heatmap <- Expr_metabric_heatmap[, rownames(AnnAll_metabric)]
Expr_metabric_heatmap_t <- t(Expr_metabric_heatmap)

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
heat_metabric <- pheatmap(Expr_metabric_heatmap, 
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

clusters_metabric <- as.data.frame(cbind(t(Expr_metabric_heatmap), 
                                         'THR clusters' = cutree(heat_metabric$tree_col, 
                                                                 k = 5)))
# plot dendrogram
#plot(heat_metabric$tree_col, labels = FALSE, hang = -1)

table(clusters_metabric$`THR clusters`)

#Re-order original data (samples) to match ordering in heatmap (left-to-right)
#samples_ordered <- as.data.frame(sort(cutree(heat_metabric$tree_col, k=5)))
#table(samples_ordered$`sort(cutree(heat_metabric$tree_col, k = 5))`)


# heat_metabric[[2]] %>%
#   as.dendrogram() %>%
#   set("branches_k_color", k=5) %>%
#   set("labels_colors", k=5) %>%
#   set("branches_lwd", 2) %>%  
#   plot(horiz=T, lwd=0.1, nodePar = list(cex = .007))

dend <- heat_metabric[[2]] %>% as.dendrogram() 
dend %>% head
#dend %>% set("branches_k_color", k = 5) %>% plot(main = 'original')

# h = 9.71
Cluster3 <- cut(dend, h = 9.71)
T1_a <- Cluster3$lower[[1]] %>% labels()
T1_b <- Cluster3$lower[[2]] %>% labels()


# re-add the new info to clusters_metabric
clusters_metabric[rownames(clusters_metabric) %in% T1_a, ]['THR clusters'] <- 'T1_a'
clusters_metabric[rownames(clusters_metabric) %in% T1_b, ]['THR clusters'] <- 'T1_b'

# add the cluster info to the phenotype table
all(rownames(clusters_metabric) == rownames(Pheno_metabric))
Pheno_metabric$`THR clusters` <- clusters_metabric$`THR clusters`

# add the cluster info to the Ann dataframe and re-plot the heatmap
all(rownames(clusters_metabric) == rownames(AnnAll_metabric))
AnnAll_metabric$`THR clusters` <- as.factor(clusters_metabric$`THR clusters`)
table(AnnAll_metabric$`THR clusters`)
levels(AnnAll_metabric$`THR clusters`) <- c('E3', 'E1', 'E4', 'E2', 'T1_a', 'T1_b')

# re-order the annotation dataframe then the expression matrix by cluster
#AnnAll_metabric <- AnnAll_metabric[order(AnnAll_metabric$cluster, decreasing = FALSE), ]
#Expr_metabric_heatmap <- Expr_metabric_heatmap[, rownames(AnnAll_metabric)]


ann_colors$`THR clusters` <- colorRampPalette(colors = rev(brewer.pal(6,"Dark2")))(6)
levels(AnnAll_metabric$`THR clusters`) <- c('E3', 'E1', 'E4', 'E2', 'T1_a', 'T1_b')
names(ann_colors$`THR clusters`) <- levels(AnnAll_metabric$`THR clusters`)

# heatmap with clinical annotation
tiff('./figures/THR50_c3Split_noLogReg/THR50_heatmap_metabric_clusters_C3split_noLogReg.tiff', width=3000, height=2000, res = 300)
pheatmap(Expr_metabric_heatmap, 
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
         fontsize_row = 8,
         cex = 1,
         cutree_cols = 5,
         cutree_rows = 5,
         breaks = seq(-1, 1, by = 0.1),
         main = "")
dev.off()


##############################################################################################################
## survival analysis:

## Keep only the relevant information (Metastasis Event and Time)
survival_metabric <- Pheno_metabric[, c("Overall.Survival.Status", "Overall.Survival..Months.", 
                                        "Relapse.Free.Status", "Relapse.Free.Status..Months.", 
                                        "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC",
                                        "X3.Gene.classifier.subtype", "THR clusters")] 

survival_metabric$`THR clusters` <- as.factor(survival_metabric$`THR clusters`)
levels(survival_metabric$`THR clusters`) <- c('E3', 'E1', 'E4', 'E2', 'T1_a', 'T1_b')

# OS
Fit_metabric_os <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ as.factor(`THR clusters`), data = survival_metabric)

# RFS
Fit_metabric_RFS <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ as.factor(`THR clusters`), data = survival_metabric)

###########################################################################
############################################################################
# plot OS
cluster_colors <- as.vector(ann_colors$`THR clusters`)

pdf("./figures/THR50_c3Split_noLogReg/metabric_os_5clusters_clustersFromHeatmap.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           #ylim = c(0.5, 1.00),
           legend.labs = c('E3', 'E1', 'E4', 'E2', 'T1_a', 'T1_b'),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR50 clusters and OS')
dev.off()


# plot RFS
pdf("./figures/THR50_c3Split_noLogReg/metabric_rfs_5clusters_clustersFromHeatmap.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('E3', 'E1', 'E4', 'E2', 'T1_a', 'T1_b'),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR50 clusters and RFS')
dev.off()




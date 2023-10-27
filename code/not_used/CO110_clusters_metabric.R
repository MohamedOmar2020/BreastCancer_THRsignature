
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

#############################
## load the ET signatures
All <- read_xlsx('./data/ET-9 Selection Steps.xlsx')
ET60 <- All$`ET-60`[!is.na(All$`ET-60`)]

#################
# CO110 signature
THR_signature <- readxl::read_xlsx("./data/THR Signatures_sep23.xlsx")

# get CO110 50.1
THR_50 <- THR_signature$`THR-50.1`[!is.na(THR_signature$`THR-50.1`)]

THR_50 <- gsub('-', '', THR_50)


# combine
CO110 <- c(ET60, THR_50)
CO110 <- CO110[!CO110 %in% 'Gene']


###############################################
# Load the expression and pheno data
##############################################
load('./objs/forKTSP.rda')

##############################################
# fix gene names
##############################################
setdiff(CO110, rownames(Expr_metabric_refAll))

grep('^SCNN1A', rownames(Expr_metabric_refAll), value = TRUE)
grep('^TNFAIP', rownames(Expr_metabric_refAll), value = TRUE)
grep('^TBC1D30', rownames(Expr_metabric_refAll), value = TRUE)
grep('^PSMB8AS1', rownames(Expr_metabric_refAll), value = TRUE)
grep('^C12orf73', rownames(Expr_metabric_refAll), value = TRUE)
grep('^PXNAS1', rownames(Expr_metabric_refAll), value = TRUE)
grep('^FLJ36848', rownames(Expr_metabric_refAll), value = TRUE)
grep('^DSCAMAS1', rownames(Expr_metabric_refAll), value = TRUE)
grep('^LOC100131541', rownames(Expr_metabric_refAll), value = TRUE)
grep('^LOC645513', rownames(Expr_metabric_refAll), value = TRUE)
grep('^USP3AS1', rownames(Expr_metabric_refAll), value = TRUE)
grep('^KCNK6', rownames(Expr_metabric_refAll), value = TRUE)
grep('^CCDC30', rownames(Expr_metabric_refAll), value = TRUE)


# ADGRG1 is the same as GPR56
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'GPR56'] <- 'ADGRG1'

# HSA011916 is the same as CTDNEP1
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'CTDNEP1'] <- 'HSA011916'

# DENND6B is the same as FAM116B
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'FAM116B'] <- 'DENND6B'

# TENT5B is the same as FAM46B
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'FAM46B'] <- 'TENT5B'

# SDHAF3 is the same as ACN9
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'ACN9'] <- 'SDHAF3'

# EVA1A is the same as FAM176A
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'FAM176A'] <- 'EVA1A'

# MINDY1 is the same as FAM63A
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'FAM63A'] <- 'MINDY1'

# P3H2 is the same as LEPREL1
rownames(Expr_metabric_refAll)[rownames(Expr_metabric_refAll) == 'LEPREL1'] <- 'P3H2'

# filter the CO110 signatures to include only the genes present in the expr matrices
CO110_fil <- CO110[CO110 %in% rownames(Expr_metabric_refAll)]

#############################################################################################
#############################################################################################
## heatmap (CO110 50)
Expr_metabric_refAll_heatmap <- Expr_metabric_refAll[CO110_fil, ] 

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

# heatmap with clinical annotation
tiff('./figures/CO110_metabric/clusters/CO110_heatmap_metabric.tiff', width=3200, height=2500, res = 300)
pheatmap(Expr_metabric_refAll_heatmap,
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
         fontsize_row = 6,
         cex = 1,
         cutree_cols = 5,
         cutree_rows = 5,
         breaks = seq(-1, 1, by = 0.1),
         main = "")
dev.off()

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
                                         'CO110 clusters' = cutree(heat_metabric$tree_col, 
                                                                 k = 5)))

table(clusters_metabric$`CO110 clusters`)

# add the cluster info to the phenotype table
all(rownames(clusters_metabric) == rownames(Pheno_metabric))
Pheno_metabric$`CO110 clusters` <- clusters_metabric$`CO110 clusters`

# add the cluster info to the Ann dataframe and re-plot the heatmap
all(rownames(clusters_metabric) == rownames(AnnAll_metabric))
AnnAll_metabric$`CO110 clusters` <- as.factor(paste0('c', clusters_metabric$`CO110 clusters`))
table(AnnAll_metabric$`CO110 clusters`)

# re-order the annotation dataframe then the expression matrix by cluster
#AnnAll_metabric <- AnnAll_metabric[order(AnnAll_metabric$cluster, decreasing = FALSE), ]
#Expr_metabric_refAll_heatmap <- Expr_metabric_refAll_heatmap[, rownames(AnnAll_metabric)]


#ann_colors$`CO110 clusters` <- colorRampPalette(colors = rev(brewer.pal(5,"Dark2")))(5)
levels(AnnAll_metabric$`CO110 clusters`) <- c('E1', 'E2', 'E3', 'T1', 'E4')
#names(ann_colors$`CO110 clusters`) <- levels(AnnAll_metabric$`CO110 clusters`)
ann_colors$`CO110 clusters` <- c("#66A61E", "#E7298A", "#1B9E77", "#7570B3", "#D95F02") 
names(ann_colors$`CO110 clusters`) <- levels(AnnAll_metabric$`CO110 clusters`)

# heatmap with clinical annotation
tiff('./figures/CO110_metabric/clusters/CO110_heatmap_metabric_clusters.tiff', width=3200, height=2500, res = 300)
pheatmap(Expr_metabric_refAll_heatmap, 
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
         fontsize_row = 6,
         cex = 1,
         cutree_cols = 5,
         cutree_rows = 5,
         breaks = seq(-1, 1, by = 0.1),
         main = "")
dev.off()

##############################################################################################################

## Keep only the relevant information (Metastasis Event and Time)
survival_metabric <- Pheno_metabric[, c("Overall.Survival.Status", "Overall.Survival..Months.", 
                                        "Relapse.Free.Status", "Relapse.Free.Status..Months.", 
                                        "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC",
                                        "X3.Gene.classifier.subtype", "CO110 clusters")] 

survival_metabric$`CO110 clusters` <- as.factor(survival_metabric$`CO110 clusters`)
levels(survival_metabric$`CO110 clusters`) <- c('E1', 'E2', 'E3', 'T1', 'E4')

# OS
Fit_metabric_os <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ as.factor(`CO110 clusters`), data = survival_metabric)

# RFS
Fit_metabric_RFS <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ as.factor(`CO110 clusters`), data = survival_metabric)

############################################################################
############################################################################
# plot OS
cluster_colors <- as.vector(ann_colors$`CO110 clusters`)

pdf("./figures/CO110_metabric/clusters/metabric_os_5clusters.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$`CO110 clusters`) ,
           legend.title	= 'CO110 clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'CO110 clusters and OS in METABRIC')
dev.off()


# plot RFS
pdf("./figures/CO110_metabric/clusters/metabric_rfs_5clusters.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = levels(survival_metabric$`CO110 clusters`) ,
           legend.title	= 'CO110 clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'CO110 clusters and RFS in METABRIC')
dev.off()

#######################################
# 10 years: 120 months
pdf("./figures/CO110_metabric/clusters/metabric_os_5clusters_10yrs.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           xlim = c(0,120),
           legend.labs = levels(survival_metabric$`CO110 clusters`) ,
           legend.title	= 'CO110 clusters',
           pval.size = 12,
           break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'CO110 clusters and 10 years OS in METABRIC')
dev.off()


# plot RFS
pdf("./figures/CO110_metabric/clusters/metabric_rfs_5clusters_10yrs.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           xlim = c(0,120),
           legend.labs = levels(survival_metabric$`CO110 clusters`) ,
           legend.title	= 'CO110 clusters',
           pval.size = 12,
           break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'CO110 clusters and 10 years RFS in METABRIC')
dev.off()

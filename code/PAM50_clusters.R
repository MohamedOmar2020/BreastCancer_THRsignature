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

#################
PAM50_signature <- readxl::read_xlsx("./data/THR_Signatures_Jan25_2023.xlsx")

# get the THR50 signature
PAM50 <- PAM50_signature$PAM50[!is.na(PAM50_signature$PAM50)]

PAM50 <- gsub('-', '', PAM50)

################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')

# fix gene names
rownames(Expr_metabric_refAll)[grep('^ZNF652', rownames(Expr_metabric_refAll))]

# filter the THR signatures to include only the genes present in the expr matrices
PAM50_fil <- PAM50[PAM50 %in% rownames(Expr_metabric_refAll)]

#############################################################################################
#############################################################################################
## heatmap (THR 70)
Expr_metabric_refAll_heatmap <- Expr_metabric_refAll[PAM50_fil, ] 

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
                HER2.Status = as.factor(HER2.Status), 
                PR.Status = as.factor(PR.Status),
                #Overall.Survival.Status = as.factor(Overall.Survival.Status),
                #Relapse.Free.Status = as.factor(Relapse.Free.Status), 
                #Tumor.Stage = as.factor(Tumor.Stage), 
                Neoplasm.Histologic.Grade = as.factor(Neoplasm.Histologic.Grade),
                Pam50...Claudin.low.subtype = as.factor(Pam50...Claudin.low.subtype)
                ) %>%
  arrange(Pam50...Claudin.low.subtype)


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
tiff('./figures/logreg/PAM50_clusters/PAM50_heatmap_metabric.tiff', width=3000, height=2000, res = 300)
pheatmap(Expr_metabric_refAll_heatmap,
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
                                         'PAM50 clusters' = cutree(heat_metabric$tree_col, 
                                                                 k = 5)))

table(clusters_metabric$`PAM50 clusters`)

# add the cluster info to the phenotype table
all(rownames(clusters_metabric) == rownames(Pheno_metabric))
Pheno_metabric$`PAM50 clusters` <- clusters_metabric$`PAM50 clusters`

# add the cluster info to the Ann dataframe and re-plot the heatmap
all(rownames(clusters_metabric) == rownames(AnnAll_metabric))
AnnAll_metabric$`PAM50 clusters` <- as.factor(paste0('c', clusters_metabric$`PAM50 clusters`))
table(AnnAll_metabric$`PAM50 clusters`)

# re-order the annotation dataframe then the expression matrix by cluster
#AnnAll_metabric <- AnnAll_metabric[order(AnnAll_metabric$cluster, decreasing = FALSE), ]
#Expr_metabric_refAll_heatmap <- Expr_metabric_refAll_heatmap[, rownames(AnnAll_metabric)]


ann_colors$`PAM50 clusters` <- colorRampPalette(colors = rev(brewer.pal(5,"Dark2")))(5)
#levels(AnnAll_metabric$`PAM50 clusters`) <- c('E3', 'E1', 'E2', 'E4', 'T1')
# fix the color pallete to match that of THR50
names(ann_colors$`PAM50 clusters`) <- levels(AnnAll_metabric$`PAM50 clusters`)

# heatmap with clinical annotation
tiff('./figures/logreg/PAM50_clusters/PAM50_heatmap_metabric_clusters.tiff', width=3000, height=2000, res = 300)
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
         fontsize_row = 5,
         cex = 1,
         cutree_cols = 5,
         cutree_rows = 5,
         breaks = seq(-1, 1, by = 0.1),
         main = "")
dev.off()

#############################################################################################################
##############################################################################################################

## Keep only the relevant information (Metastasis Event and Time)
survival_metabric <- Pheno_metabric[, c("Overall.Survival.Status", "Overall.Survival..Months.", 
                                        "Relapse.Free.Status", "Relapse.Free.Status..Months.", 
                                        "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC",
                                        "X3.Gene.classifier.subtype", "PAM50 clusters")] 

survival_metabric$`PAM50 clusters` <- as.factor(survival_metabric$`PAM50 clusters`)
#levels(survival_metabric$`THR clusters`) <- c('E3', 'E1', 'E2', 'E4', 'T1')

# OS
Fit_metabric_os <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ as.factor(`PAM50 clusters`), data = survival_metabric)

# RFS
Fit_metabric_RFS <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ as.factor(`PAM50 clusters`), data = survival_metabric)

# os with original PAM50
Fit_metabric_os_PAM50original <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Pam50...Claudin.low.subtype, data = survival_metabric)

# RFS with original PAM50
Fit_metabric_rfs_PAM50original <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Pam50...Claudin.low.subtype, data = survival_metabric)

############################################################################
############################################################################
# plot OS
cluster_colors <- as.vector(ann_colors$`PAM50 clusters`)

pdf("./figures/logreg/PAM50_clusters/metabric_os_PAM50Clusters.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('c1', 'c2', 'c3', 'c4', 'c5'),
           legend.title	= 'PAM50 clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'PAM50 clusters and OS')
dev.off()


# plot RFS
pdf("./figures/logreg/PAM50_clusters/metabric_rfs_PAM50Clusters.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('c1', 'c2', 'c3', 'c4', 'c5'),
           legend.title	= 'PAM50 clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'PAM50 clusters and RFS')
dev.off()

##############################
## original PAM50 variable

# OS
pdf("./figures/logreg/PAM50_clusters/metabric_os_PAM50original.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_PAM50original,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = c('Basal', 'Claudin-low', 'Her2+', 'Luminal A', 'Luminal B'),
           legend.title	= 'clusters',
           pval.size = 12,
           break.x.by = 40,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(17, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'Original PAM50 subtypes and OS')
dev.off()


# RFS
png("./figures/logreg/PAM50_clusters/metabric_rfs_PAM50original.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_metabric_rfs_PAM50original,
           risk.table = FALSE,
           pval = FALSE,
           #palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = c('Basal', 'Claudin-low', 'Her2+', 'Luminal A', 'Luminal B'),
           legend.title	= '',
           pval.size = 11,
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
           #title = 'Original PAM50 subtypes and RFS'
           ) + guides(
             colour = guide_legend(ncol = 2))
dev.off()

############################################################
## original X3
############################################################
# os with X3
Fit_metabric_os_X3 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ X3.Gene.classifier.subtype, data = survival_metabric)

# RFS with X3
Fit_metabric_rfs_X3 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ X3.Gene.classifier.subtype, data = survival_metabric)

# OS
tiff("./figures/logreg/PAM50_clusters/metabric_os_X3.tiff", width = 2000, height = 2000, res = 300)
ggsurvplot(Fit_metabric_os_X3,
           risk.table = FALSE,
           pval = FALSE,
           #palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = c('TNBC', 'ER+ high proliferation', 'ER+ low proliferation', 'HER2+'),
           legend.title	= '',
           pval.size = 12,
           break.x.by = 40,
           ggtheme = theme_survminer(
             base_size = 18,
             font.x = c(18, 'bold.italic', 'black'),
             font.y = c(18, 'bold.italic', 'black'),
             font.tickslab = c(18, 'plain', 'black'),
             font.legend = c(17, 'bold', 'black'),
             axis.line = element_line(colour = "black"),
             panel.grid.major = element_line(colour = "grey90"),
             panel.grid.minor = element_line(colour = "grey90"),
             panel.border = element_blank(),
             panel.background = element_blank()
             ),
           
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE
           ) + guides(
             colour = guide_legend(ncol = 2))
dev.off()


# RFS
tiff("./figures/logreg/PAM50_clusters/metabric_rfs_X3.tiff", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_metabric_rfs_X3,
           risk.table = FALSE,
           pval = FALSE,
           #palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = c('TNBC', 'ER+ HP', 'ER+ LP', 'HER2+'),
           legend.title	= '',
           pval.size = 11,
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
           #title = 'Original PAM50 subtypes and RFS'
) + guides(
  colour = guide_legend(ncol = 2))
dev.off()

########################################################################
## basal vs Claudin
########################################################################
# OS

table(survival_metabric$Pam50...Claudin.low.subtype)

# keep only basal and claudin subtypes
survival_metabric_basalClaudin <- survival_metabric[survival_metabric$Pam50...Claudin.low.subtype %in% c("Basal", "claudin-low"), ] 
survival_metabric_basalClaudin$Pam50...Claudin.low.subtype <- as.factor(survival_metabric_basalClaudin$Pam50...Claudin.low.subtype)


# RFS
Fit_sig_metabric_rfs_basalClaudin <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Pam50...Claudin.low.subtype, data = survival_metabric_basalClaudin)

# RFS COXPH
survival_metabric_basalClaudin2 <- survival_metabric_basalClaudin
survival_metabric_basalClaudin2$Pam50...Claudin.low.subtype <- factor(survival_metabric_basalClaudin2$Pam50...Claudin.low.subtype, levels = c('claudin-low', 'Basal'))
Fit_sig_metabric_rfs_basalClaudin_coxph <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Pam50...Claudin.low.subtype, data = survival_metabric_basalClaudin2)
summary(Fit_sig_metabric_rfs_basalClaudin_coxph)


# RFS
png("./figures/logreg/PAM50_clusters/metabric_rfs_Basal_vs_claudin.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_sig_metabric_rfs_basalClaudin,
           risk.table = FALSE,
           pval = FALSE,
           palette = 'jco',
           xlim = c(0,240),
           legend.labs = c('Basal', 'Claudin-low'),
           legend.title	= '',
           pval.size = 11,
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
           #title = 'Original PAM50 subtypes and RFS'
) + guides(
  colour = guide_legend(ncol = 2))
dev.off()


########################################################################
## luminal A vs B
########################################################################
# OS

table(survival_metabric$Pam50...Claudin.low.subtype)

# keep only basal and claudin subtypes
survival_metabric_luminal <- survival_metabric[survival_metabric$Pam50...Claudin.low.subtype %in% c("LumA", "LumB"), ] 
survival_metabric_luminal$Pam50...Claudin.low.subtype <- factor(survival_metabric_luminal$Pam50...Claudin.low.subtype, levels = c('LumB', 'LumA'))
table(survival_metabric_luminal$Pam50...Claudin.low.subtype)


# RFS
Fit_sig_metabric_rfs_luminal <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Pam50...Claudin.low.subtype, data = survival_metabric_luminal)

# RFS COXPH
survival_metabric_luminal2 <- survival_metabric_luminal
survival_metabric_luminal2$Pam50...Claudin.low.subtype <- factor(survival_metabric_luminal2$Pam50...Claudin.low.subtype, levels = c('LumA', 'LumB'))
Fit_sig_metabric_rfs_luminal_coxph <- coxph(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Pam50...Claudin.low.subtype, data = survival_metabric_luminal2)
summary(Fit_sig_metabric_rfs_luminal_coxph)

# RFS
png("./figures/logreg/PAM50_clusters/metabric_rfs_LuminalA_vs_B.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_sig_metabric_rfs_luminal,
           risk.table = FALSE,
           pval = FALSE,
           #palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = c('Luminal B', 'Luminal A'),
           legend.title	= '',
           pval.size = 11,
           break.x.by = 40,
           palette = 'jco',
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
           #title = 'Original PAM50 subtypes and RFS'
) + guides(
  colour = guide_legend(ncol = 2))
dev.off()



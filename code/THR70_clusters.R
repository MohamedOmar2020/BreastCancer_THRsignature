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
THR_signature <- readxl::read_xlsx("./data/THR_Signatures_Jan25_2023.xlsx")

# get the THR50 signature
THR_70 <- THR_signature$`THR-70`[!is.na(THR_signature$`THR-70`)]

THR_70 <- gsub('-', '', THR_70)

################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')

# fix gene names
rownames(Expr_metabric_refAll)[grep('^ZNF652', rownames(Expr_metabric_refAll))]

# filter the THR signatures to include only the genes present in the expr matrices
THR_70_fil <- THR_70[THR_70 %in% rownames(Expr_metabric_refAll)]

#############################################################################################
#############################################################################################
## heatmap (THR 70)
Expr_metabric_refAll_heatmap <- Expr_metabric_refAll[THR_70_fil, ] 

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
tiff('./figures/logreg/THR70_clusters/THR70_heatmap_metabric.tiff', width=3000, height=2000, res = 300)
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
levels(AnnAll_metabric$`THR clusters`) <- c('E2b', 'E2a', 'E1', 'E3', 'PQNBC')
# fix the color pallete to match that of THR50
ann_colors$`THR clusters` <- c("#66A61E", "#E7298A", "#1B9E77" , "#D95F02", "#7570B3") 
names(ann_colors$`THR clusters`) <- levels(AnnAll_metabric$`THR clusters`)
table(AnnAll_metabric$`THR clusters`)

# heatmap with clinical annotation
tiff('./figures/logreg/THR70_clusters/THR70_heatmap_metabric_clusters.tiff', width=3000, height=2000, res = 300)
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
                                        "X3.Gene.classifier.subtype", "THR clusters")] 

survival_metabric$`THR clusters` <- as.factor(survival_metabric$`THR clusters`)
levels(survival_metabric$`THR clusters`) <- c('E2b', 'E2a', 'E1', 'E3', 'PNBC')
survival_metabric$`THR clusters` <- factor(survival_metabric$`THR clusters`, levels = c('E1', 'E2a', 'E2b', 'E3', 'PNBC'))
# OS
Fit_metabric_os <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ as.factor(`THR clusters`), data = survival_metabric)

# RFS
Fit_metabric_RFS <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ as.factor(`THR clusters`), data = survival_metabric)

####################
# with E2 merged
survival_metabric2 <- survival_metabric
levels(survival_metabric2$`THR clusters`) <- c("E1", "E2", 'E2', 'E3', 'PQNBC') 

# OS
Fit_metabric_os_2 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~  as.factor(`THR clusters`), data = survival_metabric2)

# RFS
Fit_metabric_RFS_2 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~  as.factor(`THR clusters`), data = survival_metabric2)

############################################################################
############################################################################
# plot OS
#cluster_colors <- as.vector(ann_colors$`THR clusters`)
#cluster_colors <- levels(survival_metabric$`THR clusters`)

cluster_colors <- c("#1B9E77", "#E7298A", "#66A61E" , "#D95F02", "#7570B3") 
cluster_colors2 <- c('#1B9E77', "#c1b026", "#D95F02", '#7570B3')

#names(ann_colors$`THR clusters`) <- levels(survival_metabric$`THR clusters`)
#cluster_colors <- as.vector(ann_colors$`THR clusters`)


pdf("./figures/logreg/THR70_clusters/metabric_os_5clusters_20yrs.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = levels(survival_metabric$`THR clusters`),
           legend.title	= 'THR clusters',
           pval.size = 12,
           break.x.by = 40,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR70 clusters and OS'
           )
dev.off()


# plot RFS
png("./figures/logreg/THR70_clusters/metabric_rfs_5clusters_20yrs.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_metabric_RFS,
           risk.table = FALSE,
           pval = TRUE,
           palette = cluster_colors,
           xlim = c(0,240),
           legend.labs = c("E1", "E2a", "E2b", "E3", "PQNBC"),
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

## same with merged E2a and E2b
png("./figures/logreg/THR70_clusters/metabric_rfs_5clusters_20yrs_E2merged.png", width = 2000, height = 2000, res = 350)
ggsurvplot(Fit_metabric_RFS_2,
           risk.table = FALSE,
           pval = TRUE,
           palette =  cluster_colors2,
           xlim = c(0,240),
           legend.labs = c("E1", "E2","E3", "PQNBC"),
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
           #title = 'THR70 clusters and RFS: THR70 + I20'
) + guides(
  colour = guide_legend(ncol = 2))
dev.off()

###############################################################
# Venn diagram of PAM50 groups within THR-70 clusters
library(ggvenn)
library(UpSetR)

## prepare the data
# create a dataframe
df <- data.frame(
  'PAM50 groups' = survival_metabric2$Pam50...Claudin.low.subtype,
  'THR groups' = survival_metabric2$`THR clusters`,
  'sample_id' = rownames(survival_metabric2)
)

# Convert each unique group to a list
pam50_list <- split(df$sample_id, df$PAM50.groups)
thr_list <- split(df$sample_id, df$THR.groups)

# Combine the lists
group_lists <- c(pam50_list, thr_list)

ggvenn(
  group_lists,
  #fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  #stroke_size = 0.5, set_name_size = 4
)

# Convert the list to a data frame suitable for UpSet plot
upset_data <- fromList(group_lists)

# Generate the UpSet plot
png("./figures/logreg/THR70_clusters/THR70_PAM50_intersection.png", width = 2500, height = 2000, res = 300)
upset(upset_data, 
      nsets = 9, 
      sets = rev(c('LumA', 'LumB', 'Her2', 'Basal', 'claudin-low', 'E1', 'E2', 'E3', 'PQNBC')),
      mainbar.y.label = "N of samples in common",
      sets.x.label = "N of samples per group",
      keep.order = TRUE,
      sets.bar.color = c('Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Red', 'Red', 'Red', 'Red')
      )
dev.off()




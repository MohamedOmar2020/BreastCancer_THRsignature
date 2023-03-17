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

# get the THR70 signature
THR_70 <- THR_signature$`THR-70`[!is.na(THR_signature$`THR-70`)]

THR_70 <- gsub('-', '', THR_70)

################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')

# fix gene names
rownames(Expr_metabric)[grep('^ZNF652', rownames(Expr_metabric))]

# filter the THR signatures to include only the genes present in the expr matrices
THR_70_fil <- THR_70[THR_70 %in% rownames(Expr_metabric)]

#############################################################################################
#############################################################################################
## heatmap (THR 50)
Expr_metabric_heatmap <- Expr_metabric[THR_70_fil, ] 

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
#Expr_metabric_heatmap <- Expr_metabric_heatmap[, rownames(AnnAll_metabric)]


ann_colors$`THR clusters` <- colorRampPalette(colors = rev(brewer.pal(5,"Dark2")))(5)
levels(AnnAll_metabric$`THR clusters`) <- c('E3', 'E1', 'T1', 'E4', 'E2')
names(ann_colors$`THR clusters`) <- levels(AnnAll_metabric$`THR clusters`)


#############################################################################################################
# get cluster 3 with crossing curves
cluster3_pheno <- Pheno_metabric[Pheno_metabric$`THR clusters` == '3', ]
cluster3_expr <- Expr_metabric[, rownames(cluster3_pheno)]

all(rownames(cluster3_pheno) == colnames(cluster3_expr))

######################################
## divide based on survival (RFS)
summary(cluster3_pheno$Relapse.Free.Status..Months.)
####
# c3 summary
# sumtable(Pheno_metabric,
#          group = 'THR clusters',
#          file='metabric_clusters_summary',
#          out = 'browser',
#          title='METABRIC clusters Summary Statistics',
#          simple.kable=FALSE,
#          opts=list())

#################################################################################################
## training

### combine in 1 dataset: Training
RFS_c3 <- as.factor(cluster3_pheno$Relapse.Free.Status)
Data_c3 <- as.data.frame(cbind(t(cluster3_expr), RFS_c3))
Data_c3$RFS_c3 <- as.factor(Data_c3$RFS_c3)
levels(Data_c3$RFS_c3) <- c('0', '1')

#####################################
# the model
THR70_model_c3 <- glm(as.formula((paste("RFS_c3 ~", paste(THR_70_fil, collapse = "+")))), data = Data_c3, family = "binomial")
summary(THR70_model_c3)

############################################################################
# Make predictions

Train_prob_THR70_c3 <- THR70_model_c3 %>% predict(Data_c3 , type = "response")

### Threshold
thr_THR70_c3 <- coords(roc(RFS_c3, Train_prob_THR70_c3, levels = c("0", "1"), direction = "<"), "best")["threshold"]
thr_THR70_c3

### ROC Curve
ROCTrain_THR70_c3 <- roc(RFS_c3, Train_prob_THR70_c3, plot = F, print.thres=thr_THR70_c3$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("0", "1"), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR70_c3

### Get predictions based on best threshold from ROC curve
predClasses_THR70_c3 <- ifelse(Train_prob_THR70_c3 >= thr_THR70_c3$threshold, "1", "0")
table(predClasses_THR70_c3)
predClasses_THR70_c3 <- factor(predClasses_THR70_c3, levels = c('0', '1'))

##########################
## Keep only the relevant information (Metastasis Event and Time)
cluster3_pheno <- cbind(cluster3_pheno[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Relapse.Free.Status", "Relapse.Free.Status..Months.", "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC", "X3.Gene.classifier.subtype")], 
                        Train_prob_THR70_c3, predClasses_THR70_c3)


CoxData_metabric_c3 <- data.frame(cluster3_pheno)

##########################################################################################
##########################################################################################
##########################################################################################
## survival analysis for just c3

# OS
Fit_sig_metabric_os_THR70_c3 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR70_c3, data = CoxData_metabric_c3)

# RFS
Fit_sig_metabric_RFS_THR70_c3 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR70_c3, data = CoxData_metabric_c3)


# plot OS
tiff("./figures/c3_logReg_RFS_THR70/THR70_metabric_os_c3.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR70_c3,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC c3'
)
dev.off()

######################################
# plot RFS
tiff("./figures/c3_logReg_RFS_THR70/THR70_metabric_RFS_c3.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR70_c3,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC c3'
)
dev.off()

##########################################################################################
##########################################################################################
## recombine c3 with the rest 
cluster3_pheno$`THR clusters` <- as.factor(cluster3_pheno$predClasses_THR70_c3)
levels(cluster3_pheno$`THR clusters`) <- c('3_1', '3_2')

c3 <- data.frame(`THR clusters` = cluster3_pheno$`THR clusters`, `Sample.ID` = rownames(cluster3_pheno))
rownames(c3) <- rownames(cluster3_pheno)


# merge
Pheno_metabric$`THR clusters`[Pheno_metabric$`THR clusters` == '3'] <- NA
Pheno_metabric2 <- merge(x = c3, y = Pheno_metabric, by="Sample.ID", all.y = TRUE)

Pheno_metabric2 <- Pheno_metabric2 %>% 
  mutate(`THR clusters` = as.factor(`THR clusters`), THR.clusters = as.factor(THR.clusters)) %>%
  mutate(THR.clusters = coalesce(THR.clusters,`THR clusters`))

###########################################################################################
##########################################################################################
## survival analysis

## Keep only the relevant information (Metastasis Event and Time)
survival_metabric <- Pheno_metabric2[, c("Overall.Survival.Status", "Overall.Survival..Months.", 
                                         "Relapse.Free.Status", "Relapse.Free.Status..Months.", 
                                         "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC",
                                         "X3.Gene.classifier.subtype", "THR.clusters")] 

survival_metabric$THR.clusters <- as.factor(survival_metabric$THR.clusters)
levels(survival_metabric$THR.clusters) <- c('T1_a', 'T1_b', 'E3', 'E1', 'E4', 'E2')


cluster_colors <- as.vector(ann_colors$THR.clusters)

# OS
Fit_metabric_os <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters, data = survival_metabric)

# RFS
Fit_metabric_RFS <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters, data = survival_metabric)

pdf("./figures/c3_logReg_RFS_THR70/metabric_os_5clusters_newC3.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('T1_a', 'Ta_b', 'E3', 'E1', 'E4', 'E2'),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR70 clusters and OS')
dev.off()

## RFS: 
pdf("./figures/c3_logReg_RFS_THR70/metabric_rfs_5clusters_newC3_10yrs.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           xlim = c(0,120),
           legend.labs = c('T1_a', 'T1_b', 'E3', 'E1', 'E4', 'E2'),
           legend.title	= 'THR clusters',
           pval.size = 12,
           break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR70 clusters and RFS')
dev.off()


##############################################################
## heatmap

Expr_metabric_heatmap <- Expr_metabric[THR_70_fil, ] 

# Create annotation for columns/samples based on some clinical variables:
Pheno_metabric3 <- merge(x = c3, y = Pheno_metabric, by="Sample.ID", all.y = TRUE)

Pheno_metabric3 <- Pheno_metabric3 %>% 
  mutate(`THR clusters` = as.factor(`THR clusters`), THR.clusters = as.factor(THR.clusters)) %>%
  mutate(THR.clusters = coalesce(THR.clusters,`THR clusters`))

rownames(Pheno_metabric3) <- Pheno_metabric3$Sample.ID

Pheno_metabric_forHeatmap <- Pheno_metabric3
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

# add the cluster info to the Ann dataframe and re-plot the heatmap
all(rownames(Pheno_metabric3) == rownames(AnnAll_metabric))
levels(Pheno_metabric3$THR.clusters) <- c('T1_a', 'T1_b', 'E3', 'E1', 'E4', 'E2')
AnnAll_metabric$THR.clusters <- Pheno_metabric3$THR.clusters
table(AnnAll_metabric$THR.clusters)

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

ann_colors$THR.clusters <- colorRampPalette(colors = rev(brewer.pal(5,"Dark2")))(6)
#levels(AnnAll_metabric$`THR clusters`) <- c('E3', 'E1', 'T1', 'E4', 'E2')
names(ann_colors$THR.clusters) <- levels(AnnAll_metabric$THR.clusters)


breaksList = seq(-4, 4, by = 1)
ColPal <- colorRampPalette(colors = rev(brewer.pal(11,"RdYlBu")))(20)
ColPal2 <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(20))

# heatmap with clinical annotation
tiff('./figures/c3/THR70_heatmap_metabric_clusters_newC3.tiff', width=3000, height=2000, res = 300)
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





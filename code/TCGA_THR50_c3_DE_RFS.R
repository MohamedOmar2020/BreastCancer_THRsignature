

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
rownames(Expr_tcga)[grep('^ZNF652', rownames(Expr_tcga))]

# filter the THR signatures to include only the genes present in the expr matrices
THR_50_fil <- THR_50[THR_50 %in% rownames(Expr_tcga)]

#######################################
## read TCGA PAM50
Pam50_tcga <- read.delim('data/brca_tcga/data_clinical_patient.txt') 
Pam50_tcga <- Pam50_tcga[-c(1:4), ]
Pam50_tcga$Patient.ID <- gsub("\\-", "\\.", Pam50_tcga$X.Patient.Identifier)
Pam50_tcga <- Pam50_tcga[, c('Patient.ID', 'Subtype')]

# merge with TCGA pheno table
summary(Pheno_tcga$Patient.ID %in% Pam50_tcga$Patient.ID)

Pheno_tcga <- merge(Pheno_tcga, Pam50_tcga, by = "Patient.ID")

table(Pheno_tcga$Subtype)
table(Pam50_tcga$Subtype)


#############################################################################################
#############################################################################################
## heatmap (THR 50)
Expr_tcga_heatmap <- Expr_tcga[THR_50_fil, ] 

# Create annotation for columns/samples based on some clinical variables:
Pheno_tcga_forHeatmap <- Pheno_tcga
rownames(Pheno_tcga_forHeatmap) <- NULL

AnnAll_tcga <- Pheno_tcga_forHeatmap %>% 
  as.data.frame() %>%
  dplyr::select(Patient.ID, ER.Status.By.IHC, PR.status.by.ihc, HER2.fish.status, IHC.HER2, Subtype) %>%
  filter(Subtype %in% c('BRCA_Basal', 'BRCA_Her2', 'BRCA_LumA', 'BRCA_LumB')) %>%
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
Expr_tcga_heatmap <- Expr_tcga_heatmap[, rownames(AnnAll_tcga)]
Expr_tcga_heatmap_t <- t(Expr_tcga_heatmap)

# filter pheno (above we remove normal and NC)
Pheno_tcga <- Pheno_tcga[rownames(AnnAll_tcga), ]

# colors
ann_colors = list()
ann_colors$Pam50_subtypes <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(4)
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
heat_tcga <- pheatmap(Expr_tcga_heatmap, 
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

clusters_tcga <- as.data.frame(cbind(t(Expr_tcga_heatmap), 
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
#Expr_tcga_heatmap <- Expr_tcga_heatmap[, rownames(AnnAll_tcga)]


ann_colors$`THR clusters` <- colorRampPalette(colors = rev(brewer.pal(5,"Dark2")))(5)
names(ann_colors$`THR clusters`) <- levels(AnnAll_tcga$`THR clusters`)

# heatmap with cluster annotation
tiff('./figures/tcga/THR50_original_clusters/THR50_heatmap_tcga_clusters.tiff', width=3000, height=2000, res = 300)
pheatmap(Expr_tcga_heatmap, 
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
         fontsize_row = 8,
         cex = 1,
         cutree_cols = 5,
         cutree_rows = 5,
         breaks = seq(-1, 1, by = 0.1),
         main = "")
dev.off()

# another one with T and E annotation
levels(AnnAll_tcga$`THR clusters`) <- c('E3', 'E1', 'E4', 'E2', 'T1')
ann_colors$`THR clusters` <- c('#66a61e', "#E7298A", '#D95F02', "#1B9E77", '#7570b3')

names(ann_colors$`THR clusters`) <- levels(AnnAll_tcga$`THR clusters`)

tiff('./figures/tcga/THR50_original_clusters/THR50_heatmap_tcga_clusters_E.tiff', width=3000, height=2000, res = 300)
pheatmap(Expr_tcga_heatmap, 
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
         fontsize_row = 8,
         cex = 1,
         cutree_cols = 5,
         cutree_rows = 5,
         breaks = seq(-1, 1, by = 0.1),
         main = "")
dev.off()


#############################################################################################################
# get cluster 3 with crossing curves
cluster3_pheno <- Pheno_tcga[Pheno_tcga$`THR clusters` == '3', ]
cluster3_expr <- Expr_tcga[, rownames(cluster3_pheno)]

all(rownames(cluster3_pheno) == colnames(cluster3_expr))

######################################
## divide based on survival (RFS)
summary(cluster3_pheno$Relapse.Free.Status..Months.)

summary(cluster3_pheno$Relapse.Free.Status..Months. >= 50)

cluster3_pheno$c3_rfs_binary <- ifelse(cluster3_pheno$Relapse.Free.Status..Months. >= 50, 'longSurv', 'shortSurv')
table(cluster3_pheno$c3_rfs_binary)

######################################
# differential expression 
design <- model.matrix( ~ cluster3_pheno$c3_rfs_binary)
colnames(design)[2] <- "longVSshortSurvival"

fit <- lmFit(cluster3_expr, design)

fitted.ebayes <- eBayes(fit)

c3_top20 <- topTable(fitted.ebayes, number = 20)
c3_top20$gene <- rownames(c3_top20)

c3_top200 <- topTable(fitted.ebayes, number = 200)
c3_top200$gene <- rownames(c3_top200)

#write_csv(as.data.frame(c3_top20), file = './figures/c3_DE_THR50_RFS/c3_longVSshortSurv_DE.csv')
library("writexl")
write_xlsx(c3_top20,"./figures/c3_DE_THR50_RFS/c3_longVSshortSurv_DE.xlsx")

# save top200 DE genes
write_xlsx(c3_top200,"./figures/c3_DE_THR50_RFS/c3_longVSshortSurv_DE_top200.xlsx")





c3_gns <- rownames(topTable(fitted.ebayes, number = 20))

# genes in common with THR50
summary(c3_gns %in% THR_50)

#summary(decideTests(fitted.ebayes[,"longVSshortSurvival"],lfc=0))



####
# c3 summary
# sumtable(Pheno_tcga,
#          group = 'THR clusters',
#          file='tcga_clusters_summary',
#          out = 'browser',
#          title='tcga clusters Summary Statistics',
#          simple.kable=FALSE,
#          opts=list())

#################################################################################################
## training

### combine in 1 dataset: Training
RFS_c3 <- as.factor(cluster3_pheno$c3_rfs_binary)
Data_c3 <- as.data.frame(cbind(t(cluster3_expr), RFS_c3))
Data_c3$RFS_c3 <- as.factor(Data_c3$RFS_c3)
levels(Data_c3$RFS_c3) <- c('longSurv', 'shortSurv')
table(Data_c3$RFS_c3)

# the model
model20_c3 <- glm(as.formula((paste("RFS_c3 ~", paste(c3_gns, collapse = "+")))), data = Data_c3, family = "binomial")
summary(model20_c3)

save(model20_c3, file = './objs/THR50_20_model_c3.rda')

#####################################
# the model
############################################################################
# 
# pred_c3 <- as.matrix(t(cluster3_expr))
# fit <- cv.glmnet(x=pred_c3, y=RFS_c3, type.measure = "class", alpha = 0.5, family="binomial", nlambda=200)
# 
# plot(fit)
# print(fit)
# 
# 
# tmp_coeffs <- coef(fit, s=fit$lambda.1se)
# Predictor_genes <- data.frame(name=tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
# Predictor_genes_sorted <- Predictor_genes[order(Predictor_genes$coefficient, decreasing = TRUE),]
# save(Predictor_genes_sorted, file = "./Objs/Predictor_genes_ENR.rda")
# 
# 
# lambda <- cv.glmnet(pred_c3, RFS_c3, alpha = 1, family="binomial")$lambda.1se
# mod <- glmnet(pred_c3, RFS_c3, lambda = lambda, alpha = 1, family="binomial" )
# 
# 
# # use the broom package to get the row names you want:
# broom::tidy(mod) %>% 
#   slice(-1) %>%  # drop the intercept
#   arrange(desc(abs(estimate)))


############################################################################
# Make predictions

Train_prob_THR50_c3 <- model20_c3 %>% predict(Data_c3 , type = "response")

### Threshold
thr_THR50_c3 <- coords(roc(RFS_c3, Train_prob_THR50_c3, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR50_c3

### ROC Curve
ROCTrain_THR50_c3 <- roc(RFS_c3, Train_prob_THR50_c3, plot = F, print.thres=thr_THR50_c3$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR50_c3

### Get predictions based on best threshold from ROC curve
predClasses_THR50_c3 <- ifelse(Train_prob_THR50_c3 >= thr_THR50_c3$threshold, "longSurv", "shortSurv")
table(predClasses_THR50_c3)
predClasses_THR50_c3 <- factor(predClasses_THR50_c3, levels = c('longSurv', 'shortSurv'))

##########################
## Keep only the relevant information (Metastasis Event and Time)
cluster3_pheno <- cbind(cluster3_pheno[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Relapse.Free.Status", "Relapse.Free.Status..Months.", "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC", "X3.Gene.classifier.subtype")], 
                        Train_prob_THR50_c3, predClasses_THR50_c3)


CoxData_tcga_c3 <- data.frame(cluster3_pheno)

##########################################################################################
##########################################################################################
##########################################################################################
## survival analysis for just c3

# OS
Fit_sig_tcga_os_THR50_c3 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR50_c3, data = CoxData_tcga_c3)

# RFS
Fit_sig_tcga_RFS_THR50_c3 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR50_c3, data = CoxData_tcga_c3)


# plot OS
tiff("./figures/c3_DE_THR50_RFS/THR50_tcga_os_c3.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_tcga_os_THR50_c3,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in tcga in T1 class derived from THR50'
)
dev.off()

######################################
# plot RFS
tiff("./figures/c3_DE_THR50_RFS/THR50_tcga_RFS_c3.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_tcga_RFS_THR50_c3,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in tcga in T1 class derived from THR50'
)
dev.off()

##########################################################################################
##########################################################################################
## recombine c3 with the rest 
cluster3_pheno$`THR clusters` <- as.factor(cluster3_pheno$predClasses_THR50_c3)
levels(cluster3_pheno$`THR clusters`) <- c('3_1', '3_2')

c3 <- data.frame(`THR clusters` = cluster3_pheno$`THR clusters`, `Sample.ID` = rownames(cluster3_pheno))
rownames(c3) <- rownames(cluster3_pheno)


# merge
Pheno_tcga$`THR clusters`[Pheno_tcga$`THR clusters` == '3'] <- NA
Pheno_tcga2 <- merge(x = c3, y = Pheno_tcga, by="Sample.ID", all.y = TRUE)

Pheno_tcga2 <- Pheno_tcga2 %>% 
  mutate(`THR clusters` = as.factor(`THR clusters`), THR.clusters = as.factor(THR.clusters)) %>%
  mutate(THR.clusters = coalesce(THR.clusters,`THR clusters`))

###########################################################################################
##########################################################################################
## survival analysis

## Keep only the relevant information (Metastasis Event and Time)
survival_tcga <- Pheno_tcga2[, c("Overall.Survival.Status", "Overall.Survival..Months.", 
                                         "Relapse.Free.Status", "Relapse.Free.Status..Months.", 
                                         "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC",
                                         "X3.Gene.classifier.subtype", "THR.clusters")] 

survival_tcga$THR.clusters <- as.factor(survival_tcga$THR.clusters)
levels(survival_tcga$THR.clusters) <- c('T1_a', 'T1_b', 'E3', 'E1', 'E4', 'E2')


cluster_colors <- as.vector(ann_colors$THR.clusters)

# OS
Fit_tcga_os <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters, data = survival_tcga)

# RFS
Fit_tcga_RFS <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters, data = survival_tcga)

pdf("./figures/c3_DE_THR50_RFS/tcga_os_5clusters_newC3.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_tcga_os,
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
           risk.table.y.text = FALSE, title = 'THR50 clusters and OS')
dev.off()

## RFS: 
pdf("./figures/c3_DE_THR50_RFS/tcga_rfs_5clusters_newC3.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_tcga_RFS,
           risk.table = FALSE,
           pval = TRUE,
           #palette = cluster_colors,
           #xlim = c(0,120),
           legend.labs = c('T1_a', 'T1_b', 'E3', 'E1', 'E4', 'E2'),
           legend.title	= 'THR clusters',
           pval.size = 12,
           #break.x.by = 20,
           ggtheme = theme_survminer(base_size = 18, font.x = c(18, 'bold.italic', 'black'), font.y = c(18, 'bold.italic', 'black'), font.tickslab = c(18, 'plain', 'black'), font.legend = c(18, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR50 clusters and RFS')
dev.off()


##############################################################
## heatmap

Expr_tcga_heatmap <- Expr_tcga[THR_50_fil, ] 

# Create annotation for columns/samples based on some clinical variables:
Pheno_tcga3 <- merge(x = c3, y = Pheno_tcga, by="Sample.ID", all.y = TRUE)

Pheno_tcga3 <- Pheno_tcga3 %>% 
  mutate(`THR clusters` = as.factor(`THR clusters`), THR.clusters = as.factor(THR.clusters)) %>%
  mutate(THR.clusters = coalesce(THR.clusters,`THR clusters`))

rownames(Pheno_tcga3) <- Pheno_tcga3$Sample.ID

Pheno_tcga_forHeatmap <- Pheno_tcga3
rownames(Pheno_tcga_forHeatmap) <- NULL

AnnAll_tcga <- Pheno_tcga_forHeatmap %>% 
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
all(rownames(Pheno_tcga3) == rownames(AnnAll_tcga))
levels(Pheno_tcga3$THR.clusters) <- c('T1_a', 'T1_b', 'E3', 'E1', 'E4', 'E2')
AnnAll_tcga$THR.clusters <- Pheno_tcga3$THR.clusters
table(AnnAll_tcga$THR.clusters)

# filter and transpose the expression matrix
Expr_tcga_heatmap <- Expr_tcga_heatmap[, rownames(AnnAll_tcga)]
Expr_tcga_heatmap_t <- t(Expr_tcga_heatmap)

# filter pheno (above we remove normal and NC)
Pheno_tcga <- Pheno_tcga[rownames(AnnAll_tcga), ]

# colors
ann_colors = list()
ann_colors$Pam50...Claudin.low.subtype <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(5)
names(ann_colors$Pam50...Claudin.low.subtype) <- levels(AnnAll_tcga$Pam50...Claudin.low.subtype)

ann_colors$ER.status.measured.by.IHC <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(2)
names(ann_colors$ER.status.measured.by.IHC) <- levels(AnnAll_tcga$ER.status.measured.by.IHC)

ann_colors$X3.Gene.classifier.subtype <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(4)
names(ann_colors$X3.Gene.classifier.subtype) <- levels(AnnAll_tcga$X3.Gene.classifier.subtype)

ann_colors$Neoplasm.Histologic.Grade <- colorRampPalette(colors = rev(brewer.pal(8,"RdYlBu")))(3)
names(ann_colors$Neoplasm.Histologic.Grade) <- levels(AnnAll_tcga$Neoplasm.Histologic.Grade)

ann_colors$THR.clusters <- colorRampPalette(colors = rev(brewer.pal(5,"Dark2")))(6)
#levels(AnnAll_tcga$`THR clusters`) <- c('E3', 'E1', 'T1', 'E4', 'E2')
names(ann_colors$THR.clusters) <- levels(AnnAll_tcga$THR.clusters)


breaksList = seq(-4, 4, by = 1)
ColPal <- colorRampPalette(colors = rev(brewer.pal(11,"RdYlBu")))(20)
ColPal2 <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(20))

# heatmap with clinical annotation
tiff('./figures/c3/THR50_heatmap_tcga_clusters_newC3.tiff', width=3000, height=2000, res = 300)
pheatmap(Expr_tcga_heatmap, 
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
         fontsize_row = 8,
         cex = 1,
         cutree_cols = 5,
         cutree_rows = 5,
         breaks = seq(-1, 1, by = 0.1),
         main = "")
dev.off()



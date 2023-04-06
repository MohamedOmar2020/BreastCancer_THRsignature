

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

#############################################################################################
#############################################################################################
## heatmap (THR 50)
Expr_metabric_refAll_heatmap <- Expr_metabric_refAll[THR_50_fil, ] 

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
levels(AnnAll_metabric$`THR clusters`) <- c('E3', 'E1', 'T1', 'E4', 'E2')
names(ann_colors$`THR clusters`) <- levels(AnnAll_metabric$`THR clusters`)


#############################################################################################################
# get cluster 3 with crossing curves
cluster3_pheno <- Pheno_metabric[Pheno_metabric$`THR clusters` == '3', ]
cluster3_expr <- Expr_metabric_refAll[, rownames(cluster3_pheno)]

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
write_xlsx(c3_top20,"./figures/c3_DE_THR50_RFS/THR50_c3_longVSshortSurv_DE.xlsx")

# save top200 DE genes
write_xlsx(c3_top200,"./figures/c3_DE_THR50_RFS/THR50_c3_longVSshortSurv_DE_top200.xlsx")





c3_gns <- rownames(topTable(fitted.ebayes, number = 20))

# genes in common with THR50
summary(c3_gns %in% THR_50)

#summary(decideTests(fitted.ebayes[,"longVSshortSurvival"],lfc=0))



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
RFS_c3 <- as.factor(cluster3_pheno$c3_rfs_binary)
Data_c3 <- as.data.frame(cbind(t(cluster3_expr), RFS_c3))
Data_c3$RFS_c3 <- as.factor(Data_c3$RFS_c3)
levels(Data_c3$RFS_c3) <- c('longSurv', 'shortSurv')
table(Data_c3$RFS_c3)

# the model
model20_c3 <- glm(as.formula((paste("RFS_c3 ~", paste(c3_gns, collapse = "+")))), data = Data_c3, family = "binomial")
summary(model20_c3)
save(model20_c3, file = './objs/THR50_20_model_c3.rda')


CTLA4_model_c3 <- glm(RFS_c3 ~ CTLA4, data = Data_c3, family = "binomial")
summary(CTLA4_model_c3)

CTLA4_ICOS_CXCL13_model_c3 <- glm(RFS_c3 ~ CTLA4+ICOS+CXCL13, data = Data_c3, family = "binomial")
summary(CTLA4_ICOS_CXCL13_model_c3)

KLF7_model_c3 <- glm(RFS_c3 ~ KLF7, data = Data_c3, family = "binomial")
summary(KLF7_model_c3)

ZNF627_model_c3 <- glm(RFS_c3 ~ ZNF627, data = Data_c3, family = "binomial")
summary(ZNF627_model_c3)

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

Train_prob_THR50_c3_model20 <- model20_c3 %>% predict(Data_c3 , type = "response")
Train_prob_THR50_c3_CTLA4 <- CTLA4_model_c3 %>% predict(Data_c3 , type = "response")
Train_prob_THR50_c3_CTLA4_ICOS_CXCL13 <- CTLA4_ICOS_CXCL13_model_c3 %>% predict(Data_c3 , type = "response")
Train_prob_THR50_c3_KLF7 <- KLF7_model_c3 %>% predict(Data_c3 , type = "response")
Train_prob_THR50_c3_ZNF627 <- ZNF627_model_c3 %>% predict(Data_c3 , type = "response")

### Threshold
thr_THR50_c3_model20 <- coords(roc(RFS_c3, Train_prob_THR50_c3_model20, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR50_c3_CTLA4 <- coords(roc(RFS_c3, Train_prob_THR50_c3_CTLA4, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR50_c3_CTLA4_ICOS_CXCL13 <- coords(roc(RFS_c3, Train_prob_THR50_c3_CTLA4_ICOS_CXCL13, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR50_c3_KLF7 <- coords(roc(RFS_c3, Train_prob_THR50_c3_KLF7, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]
thr_THR50_c3_ZNF627 <- coords(roc(RFS_c3, Train_prob_THR50_c3_ZNF627, levels = c('longSurv', 'shortSurv'), direction = "<"), "best")["threshold"]


### ROC Curve
ROCTrain_THR50_c3_model20 <- roc(RFS_c3, Train_prob_THR50_c3_model20, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR50_c3_model20

ROCTrain_THR50_c3_CTLA4 <- roc(RFS_c3, Train_prob_THR50_c3_CTLA4, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR50_c3_CTLA4

ROCTrain_THR50_c3_CTLA4_ICOS_CXCL13 <- roc(RFS_c3, Train_prob_THR50_c3_CTLA4_ICOS_CXCL13, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR50_c3_CTLA4_ICOS_CXCL13

ROCTrain_THR50_c3_KLF7 <- roc(RFS_c3, Train_prob_THR50_c3_KLF7, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR50_c3_KLF7

ROCTrain_THR50_c3_ZNF627 <- roc(RFS_c3, Train_prob_THR50_c3_ZNF627, plot = F, print.auc=TRUE, print.auc.col="black", ci = T, levels = c('longSurv', 'shortSurv'), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR50_c3_ZNF627

### Get predictions based on best threshold from ROC curve
predClasses_THR50_c3_model20 <- ifelse(Train_prob_THR50_c3_model20 >= thr_THR50_c3_model20$threshold, "longSurv", "shortSurv")
table(predClasses_THR50_c3_model20)
predClasses_THR50_c3_model20 <- factor(predClasses_THR50_c3_model20, levels = c('longSurv', 'shortSurv'))

predClasses_THR50_c3_CTLA4 <- ifelse(Train_prob_THR50_c3_CTLA4 >= thr_THR50_c3_CTLA4$threshold, "longSurv", "shortSurv")
table(predClasses_THR50_c3_CTLA4)
predClasses_THR50_c3_CTLA4 <- factor(predClasses_THR50_c3_CTLA4, levels = c('longSurv', 'shortSurv'))

predClasses_THR50_c3_CTLA4_ICOS_CXCL13 <- ifelse(Train_prob_THR50_c3_CTLA4_ICOS_CXCL13 >= thr_THR50_c3_CTLA4_ICOS_CXCL13$threshold, "longSurv", "shortSurv")
table(predClasses_THR50_c3_CTLA4_ICOS_CXCL13)
predClasses_THR50_c3_CTLA4_ICOS_CXCL13 <- factor(predClasses_THR50_c3_CTLA4_ICOS_CXCL13, levels = c('longSurv', 'shortSurv'))

predClasses_THR50_c3_KLF7 <- ifelse(Train_prob_THR50_c3_KLF7 >= thr_THR50_c3_KLF7$threshold, "longSurv", "shortSurv")
table(predClasses_THR50_c3_KLF7)
predClasses_THR50_c3_KLF7 <- factor(predClasses_THR50_c3_KLF7, levels = c('longSurv', 'shortSurv'))

predClasses_THR50_c3_ZNF627 <- ifelse(Train_prob_THR50_c3_ZNF627 >= thr_THR50_c3_ZNF627$threshold, "longSurv", "shortSurv")
table(predClasses_THR50_c3_ZNF627)
predClasses_THR50_c3_ZNF627 <- factor(predClasses_THR50_c3_ZNF627, levels = c('longSurv', 'shortSurv'))

##########################
## Keep only the relevant information (Metastasis Event and Time)
cluster3_pheno <- cbind(cluster3_pheno[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Relapse.Free.Status", "Relapse.Free.Status..Months.", "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC", "X3.Gene.classifier.subtype")], 
                        Train_prob_THR50_c3_model20, predClasses_THR50_c3_model20,
                        Train_prob_THR50_c3_CTLA4, predClasses_THR50_c3_CTLA4,
                        Train_prob_THR50_c3_CTLA4_ICOS_CXCL13, predClasses_THR50_c3_CTLA4_ICOS_CXCL13,
                        Train_prob_THR50_c3_KLF7, predClasses_THR50_c3_KLF7,
                        Train_prob_THR50_c3_ZNF627, predClasses_THR50_c3_ZNF627
                        )


CoxData_metabric_c3 <- data.frame(cluster3_pheno)

##########################################################################################
##########################################################################################
##########################################################################################
## survival analysis for just c3

#########################################
# Model 20
#########################################

# OS
Fit_sig_metabric_os_THR50_c3_model20 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR50_c3_model20, data = CoxData_metabric_c3)

# RFS
Fit_sig_metabric_RFS_THR50_c3_model20 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR50_c3_model20, data = CoxData_metabric_c3)


# plot OS
tiff("./figures/c3_DE_THR50_RFS/THR50_metabric_os_c3_model20.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR50_c3_model20,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR50: 20 genes model'
)
dev.off()

######################################
# plot RFS
tiff("./figures/c3_DE_THR50_RFS/THR50_metabric_RFS_c3_model20.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR50_c3_model20,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR50: 20 genes model'
)
dev.off()

#########################################
# CTLA4
#########################################

# OS
Fit_sig_metabric_os_THR50_c3_CTLA4 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR50_c3_CTLA4, data = CoxData_metabric_c3)

# RFS
Fit_sig_metabric_RFS_THR50_c3_CTLA4 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR50_c3_CTLA4, data = CoxData_metabric_c3)


# plot OS
tiff("./figures/c3_DE_THR50_RFS/THR50_metabric_os_c3_CTLA4.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR50_c3_CTLA4,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR50: CTLA4'
)
dev.off()

######################################
# plot RFS
tiff("./figures/c3_DE_THR50_RFS/THR50_metabric_RFS_c3_CTLA4.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR50_c3_CTLA4,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR50: CTLA4'
)
dev.off()


#########################################
# CTLA4_ICOS_CXCL13
#########################################

# OS
Fit_sig_metabric_os_THR50_c3_CTLA4_ICOS_CXCL13 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR50_c3_CTLA4_ICOS_CXCL13, data = CoxData_metabric_c3)

# RFS
Fit_sig_metabric_RFS_THR50_c3_CTLA4_ICOS_CXCL13 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR50_c3_CTLA4_ICOS_CXCL13, data = CoxData_metabric_c3)


# plot OS
tiff("./figures/c3_DE_THR50_RFS/THR50_metabric_os_c3_CTLA4_ICOS_CXCL13.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR50_c3_CTLA4_ICOS_CXCL13,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR50: CTLA4 + ICOS + CXCL13'
)
dev.off()

######################################
# plot RFS
tiff("./figures/c3_DE_THR50_RFS/THR50_metabric_RFS_c3_CTLA4_ICOS_CXCL13.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR50_c3_CTLA4_ICOS_CXCL13,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR50: CTLA4 + ICOS + CXCL13'
)
dev.off()

#########################################
# KLF7
#########################################

# OS
Fit_sig_metabric_os_THR50_c3_KLF7 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR50_c3_KLF7, data = CoxData_metabric_c3)

# RFS
Fit_sig_metabric_RFS_THR50_c3_KLF7 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR50_c3_KLF7, data = CoxData_metabric_c3)


# plot OS
tiff("./figures/c3_DE_THR50_RFS/THR50_metabric_os_c3_KLF7.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR50_c3_KLF7,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR50: KLF7'
)
dev.off()

######################################
# plot RFS
tiff("./figures/c3_DE_THR50_RFS/THR50_metabric_RFS_c3_KLF7.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR50_c3_KLF7,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR50: KLF7'
)
dev.off()

#########################################
# KLF7
#########################################

# OS
Fit_sig_metabric_os_THR50_c3_ZNF627 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ predClasses_THR50_c3_ZNF627, data = CoxData_metabric_c3)

# RFS
Fit_sig_metabric_RFS_THR50_c3_ZNF627 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ predClasses_THR50_c3_ZNF627, data = CoxData_metabric_c3)


# plot OS
tiff("./figures/c3_DE_THR50_RFS/THR50_metabric_os_c3_ZNF627.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_THR50_c3_ZNF627,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           title = 'OS in METABRIC in T1 class derived from THR50: ZNF627'
)
dev.off()

######################################
# plot RFS
tiff("./figures/c3_DE_THR50_RFS/THR50_metabric_RFS_c3_ZNF627.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_THR50_c3_ZNF627,
           risk.table = FALSE,
           pval = TRUE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           title = 'RFS in METABRIC in T1 class derived from THR50: ZNF627'
)
dev.off()
##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################
## recombine c3 with the rest 
cluster3_pheno$THR_clusters_model20 <- as.factor(cluster3_pheno$predClasses_THR50_c3_model20)
levels(cluster3_pheno$THR_clusters_model20) <- c('3_1', '3_2')

cluster3_pheno$THR_clusters_CTLA4 <- as.factor(cluster3_pheno$predClasses_THR50_c3_CTLA4)
levels(cluster3_pheno$THR_clusters_CTLA4) <- c('3_1', '3_2')

cluster3_pheno$THR_clusters_CTLA4_ICOS_CXCL13 <- as.factor(cluster3_pheno$predClasses_THR50_c3_CTLA4_ICOS_CXCL13)
levels(cluster3_pheno$THR_clusters_CTLA4_ICOS_CXCL13) <- c('3_1', '3_2')

cluster3_pheno$THR_clusters_KLF7 <- as.factor(cluster3_pheno$predClasses_THR50_c3_KLF7)
levels(cluster3_pheno$THR_clusters_KLF7) <- c('3_1', '3_2')

cluster3_pheno$THR_clusters_ZNF627 <- as.factor(cluster3_pheno$predClasses_THR50_c3_ZNF627)
levels(cluster3_pheno$THR_clusters_ZNF627) <- c('3_1', '3_2')


c3 <- data.frame(THR_clusters_model20 = cluster3_pheno$THR_clusters_model20,
                 THR_clusters_CTLA4 = cluster3_pheno$THR_clusters_CTLA4,
                 THR_clusters_CTLA4_ICOS_CXCL13 = cluster3_pheno$THR_clusters_CTLA4_ICOS_CXCL13,
                 THR_clusters_KLF7 = cluster3_pheno$THR_clusters_KLF7,
                 THR_clusters_ZNF627 = cluster3_pheno$THR_clusters_ZNF627,
                 `Sample.ID` = rownames(cluster3_pheno))

rownames(c3) <- rownames(cluster3_pheno)


# merge
Pheno_metabric$`THR clusters`[Pheno_metabric$`THR clusters` == '3'] <- NA


Pheno_metabric2 <- merge(x = c3, y = Pheno_metabric, by="Sample.ID", all.y = TRUE)

Pheno_metabric2 <- Pheno_metabric2 %>% 
  mutate(`THR clusters` = as.factor(`THR clusters`), THR_clusters_model20 = as.factor(THR_clusters_model20)) %>%
  mutate(THR.clusters_model20Merged = coalesce(THR_clusters_model20,`THR clusters`)) %>%
  mutate(THR.clusters_CTLA4Merged = coalesce(THR_clusters_CTLA4,`THR clusters`)) %>%
  mutate(THR.clusters_CTLA4_ICOS_CXCL13Merged = coalesce(THR_clusters_CTLA4_ICOS_CXCL13,`THR clusters`)) %>%
  mutate(THR.clusters_KLF7Merged = coalesce(THR_clusters_KLF7,`THR clusters`)) %>%
  mutate(THR.clusters_ZNF627Merged = coalesce(THR_clusters_ZNF627,`THR clusters`)) 

###########################################################################################
##########################################################################################
## survival analysis

## Keep only the relevant information (Metastasis Event and Time)
survival_metabric <- Pheno_metabric2[, c("Overall.Survival.Status", "Overall.Survival..Months.", 
                                         "Relapse.Free.Status", "Relapse.Free.Status..Months.", 
                                         "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC",
                                         "X3.Gene.classifier.subtype", 
                                         "THR.clusters_model20Merged", 
                                         'THR.clusters_CTLA4Merged',
                                         'THR.clusters_CTLA4_ICOS_CXCL13Merged',
                                         'THR.clusters_KLF7Merged',
                                         'THR.clusters_ZNF627Merged'
                                         )] 

survival_metabric$THR.clusters_model20Merged <- as.factor(survival_metabric$THR.clusters_model20Merged)
survival_metabric$THR.clusters_CTLA4Merged <- as.factor(survival_metabric$THR.clusters_CTLA4Merged)
survival_metabric$THR.clusters_CTLA4_ICOS_CXCL13Merged <- as.factor(survival_metabric$THR.clusters_CTLA4_ICOS_CXCL13Merged)
survival_metabric$THR.clusters_KLF7Merged <- as.factor(survival_metabric$THR.clusters_KLF7Merged)
survival_metabric$THR.clusters_ZNF627Merged <- as.factor(survival_metabric$THR.clusters_ZNF627Merged)

levels(survival_metabric$THR.clusters_model20Merged) <- c('T1_a', 'T1_b', 'E3', 'E1', 'E4', 'E2')
levels(survival_metabric$THR.clusters_CTLA4Merged) <- c('T1_a', 'T1_b', 'E3', 'E1', 'E4', 'E2')
levels(survival_metabric$THR.clusters_CTLA4_ICOS_CXCL13Merged) <- c('T1_a', 'T1_b', 'E3', 'E1', 'E4', 'E2')
levels(survival_metabric$THR.clusters_KLF7Merged) <- c('T1_a', 'T1_b', 'E3', 'E1', 'E4', 'E2')
levels(survival_metabric$THR.clusters_ZNF627Merged) <- c('T1_a', 'T1_b', 'E3', 'E1', 'E4', 'E2')


cluster_colors <- as.vector(ann_colors$THR.clusters_model20Merged) # same order for the others



################################################################
## Survival curves: model 20
################################################################

# OS
Fit_metabric_os_model20 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_model20Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_model20 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_model20Merged, data = survival_metabric)

pdf("./figures/c3_DE_THR50_RFS/metabric_os_5clusters_model20Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_model20,
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
           risk.table.y.text = FALSE, title = 'THR50 clusters and OS: THR50 + I20')
dev.off()

## RFS: 
pdf("./figures/c3_DE_THR50_RFS/metabric_rfs_5clusters_model20Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_model20,
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
           risk.table.y.text = FALSE, title = 'THR50 clusters and RFS: THR50 + I20')
dev.off()

################################################################
## Survival curves: CTLA4
################################################################

# OS
Fit_metabric_os_CTLA4 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_CTLA4Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_CTLA4 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_CTLA4Merged, data = survival_metabric)

pdf("./figures/c3_DE_THR50_RFS/metabric_os_5clusters_CTLA4Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_CTLA4,
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
           risk.table.y.text = FALSE, title = 'THR50 clusters and OS: THR50 + CTLA4')
dev.off()

## RFS: 
pdf("./figures/c3_DE_THR50_RFS/metabric_rfs_5clusters_CTLA4Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_CTLA4,
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
           risk.table.y.text = FALSE, title = 'THR50 clusters and RFS: THR50 + CTLA4')
dev.off()


################################################################
## Survival curves: CTLA4_ICOS_CXCL13
################################################################

# OS
Fit_metabric_os_CTLA4_ICOS_CXCL13 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_CTLA4_ICOS_CXCL13Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_CTLA4_ICOS_CXCL13 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_CTLA4_ICOS_CXCL13Merged, data = survival_metabric)

pdf("./figures/c3_DE_THR50_RFS/metabric_os_5clusters_CTLA4_ICOS_CXCL13Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_CTLA4_ICOS_CXCL13,
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
           risk.table.y.text = FALSE, title = 'THR50 clusters and OS: THR50 + CTLA4 + ICOS + CXCL13')
dev.off()

## RFS: 
pdf("./figures/c3_DE_THR50_RFS/metabric_rfs_5clusters_CTLA4_ICOS_CXCL13Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_CTLA4_ICOS_CXCL13,
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
           risk.table.y.text = FALSE, title = 'THR50 clusters and RFS: THR50 + CTLA4 + ICOS + CXCL13')
dev.off()



################################################################
## Survival curves: KLF7
################################################################

# OS
Fit_metabric_os_KLF7 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_KLF7Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_KLF7 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_KLF7Merged, data = survival_metabric)

pdf("./figures/c3_DE_THR50_RFS/metabric_os_5clusters_KLF7Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_KLF7,
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
           risk.table.y.text = FALSE, title = 'THR50 clusters and OS: THR50 + KLF7')
dev.off()

## RFS: 
pdf("./figures/c3_DE_THR50_RFS/metabric_rfs_5clusters_KLF7Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_KLF7,
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
           risk.table.y.text = FALSE, title = 'THR50 clusters and RFS: THR50 + KLF7')
dev.off()


################################################################
## Survival curves: ZNF627
################################################################

# OS
Fit_metabric_os_ZNF627 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ THR.clusters_ZNF627Merged, data = survival_metabric)

# RFS
Fit_metabric_RFS_ZNF627 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ THR.clusters_ZNF627Merged, data = survival_metabric)

pdf("./figures/c3_DE_THR50_RFS/metabric_os_5clusters_ZNF627Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_os_ZNF627,
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
           risk.table.y.text = FALSE, title = 'THR50 clusters and OS: THR50 + ZNF627')
dev.off()

## RFS: 
pdf("./figures/c3_DE_THR50_RFS/metabric_rfs_5clusters_ZNF627Merged.pdf", width = 10, height = 8, onefile = F)
ggsurvplot(Fit_metabric_RFS_ZNF627,
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
           risk.table.y.text = FALSE, title = 'THR50 clusters and RFS: THR50 + ZNF627')
dev.off()











##############################################################
## heatmap

Expr_metabric_refAll_heatmap <- Expr_metabric_refAll[THR_50_fil, ] 

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

ann_colors$THR.clusters <- colorRampPalette(colors = rev(brewer.pal(5,"Dark2")))(6)
#levels(AnnAll_metabric$`THR clusters`) <- c('E3', 'E1', 'T1', 'E4', 'E2')
names(ann_colors$THR.clusters) <- levels(AnnAll_metabric$THR.clusters)


breaksList = seq(-4, 4, by = 1)
ColPal <- colorRampPalette(colors = rev(brewer.pal(11,"RdYlBu")))(20)
ColPal2 <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(20))

# heatmap with clinical annotation
tiff('./figures/c3/THR50_heatmap_metabric_clusters_newC3.tiff', width=3000, height=2000, res = 300)
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





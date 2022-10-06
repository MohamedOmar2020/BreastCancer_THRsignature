
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

#################
THR_signature <- readxl::read_xlsx("./data/THR Signatures_sep23.xlsx")

# get THR 25 and 50
THR_25 <- THR_signature$`THR-25`[!is.na(THR_signature$`THR-25`)]
THR_50_1 <- THR_signature$`THR-50.1`[!is.na(THR_signature$`THR-50.1`)]
THR_50_2 <- THR_signature$`THR-50.2`[!is.na(THR_signature$`THR-50.2`)]

THR_50_1 <- gsub('-', '', THR_50_1)
THR_50_2 <- gsub('-', '', THR_50_2)

################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')

### combine in 1 dataset: Training
Data_metabric <- as.data.frame(cbind(t(Expr_metabric), group_metabric))
Data_metabric$group_metabric <- as.factor(Data_metabric$group_metabric)
levels(Data_metabric$group_metabric) <- c('0', '1')
colnames(Data_metabric)[colnames(Data_metabric) %in% c('group_metabric')] <- c('os')

##########
### tcga
Data_tcga <- as.data.frame(cbind(t(Expr_tcga), group_tcga))
Data_tcga$group_tcga <- as.factor(Data_tcga$group_tcga)
levels(Data_tcga$group_tcga) <- c('0', '1')
colnames(Data_tcga)[colnames(Data_tcga) %in% c('group_tcga')] <- c('os')

###########################################################################
### TRAINING using logistic regression
###########################################################################

setdiff(THR_50_1, THR_50_1_fil)
setdiff(THR_50_2, THR_50_2_fil)

rownames(Expr_metabric)[grep('^ZNF652', rownames(Expr_metabric))]
rownames(Expr_tcga)[grep('^ZNF652', rownames(Expr_tcga))]


# filter the THR signatures to include only the genes present in the expr matrices
THR_25_fil <- THR_25[THR_25 %in% rownames(Expr_metabric) & THR_25 %in% rownames(Expr_tcga)]
THR_50_1_fil <- THR_50_1[THR_50_1 %in% rownames(Expr_metabric) & THR_50_1 %in% rownames(Expr_tcga)]
THR_50_2_fil <- THR_50_2[THR_50_2 %in% rownames(Expr_metabric) & THR_50_2 %in% rownames(Expr_tcga)]


THR25_model <- glm(as.formula((paste("os ~", paste(THR_25_fil, collapse = "+")))), data = Data_metabric, family = "binomial")
summary(THR25_model)

THR50_1_model <- glm(as.formula((paste("os ~", paste(THR_50_1_fil, collapse = "+")))), data = Data_metabric, family = "binomial")
summary(THR50_1_model)

THR50_2_model <- glm(as.formula((paste("os ~", paste(THR_50_2_fil, collapse = "+")))), data = Data_metabric, family = "binomial")
summary(THR50_2_model)

save(THR25_model, THR50_1_model, THR50_2_model, file = "./objs/THR_models_logreg.rda")

###########################################################################
############################################################################
### predict in the training dataset
# Make predictions

Train_prob_THR25 <- THR25_model %>% predict(Data_metabric , type = "response")
Train_prob_THR50_1 <- THR50_1_model %>% predict(Data_metabric , type = "response")
Train_prob_THR50_2 <- THR50_2_model %>% predict(Data_metabric , type = "response")

### Threshold
thr_THR25 <- coords(roc(group_metabric, Train_prob_THR25, levels = c("0", "1"), direction = "<"), "best")["threshold"]
thr_THR25

thr_THR50_1 <- coords(roc(group_metabric, Train_prob_THR50_1, levels = c("0", "1"), direction = "<"), "best")["threshold"]
thr_THR50_1

thr_THR50_2 <- coords(roc(group_metabric, Train_prob_THR50_2, levels = c("0", "1"), direction = "<"), "best")["threshold"]
thr_THR50_2

### ROC Curve
ROCTrain_THR25 <- roc(group_metabric, Train_prob_THR25, plot = F, print.thres=thr_THR25$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("0", "1"), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR25

ROCTrain_THR50_1 <- roc(group_metabric, Train_prob_THR50_1, plot = F, print.thres=thr_THR50_1$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("0", "1"), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR50_1

ROCTrain_THR50_2 <- roc(group_metabric, Train_prob_THR50_2, plot = F, print.thres=thr_THR50_2$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("0", "1"), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_THR50_2

### Get predictions based on best threshold from ROC curve
Train_predClasses_THR25 <- ifelse(Train_prob_THR25 >= thr_THR25$threshold, "1", "0")
table(Train_predClasses_THR25)
Train_predClasses_THR25 <- factor(Train_predClasses_THR25, levels = c('0', '1'))

Train_predClasses_THR50_1 <- ifelse(Train_prob_THR50_1 >= thr_THR50_1$threshold, "1", "0")
table(Train_predClasses_THR50_1)
Train_predClasses_THR50_1 <- factor(Train_predClasses_THR50_1, levels = c('0', '1'))

Train_predClasses_THR50_2 <- ifelse(Train_prob_THR50_2 >= thr_THR50_2$threshold, "1", "0")
table(Train_predClasses_THR50_2)
Train_predClasses_THR50_2 <- factor(Train_predClasses_THR50_2, levels = c('0', '1'))


### Resubstitution performance in the TRAINING set
ConfusionTrain_THR25 <- confusionMatrix(Train_predClasses_THR25, group_metabric, positive = "1", mode = "everything")
ConfusionTrain_THR25

ConfusionTrain_THR50_1 <- confusionMatrix(Train_predClasses_THR50_1, group_metabric, positive = "1", mode = "everything")
ConfusionTrain_THR50_1

ConfusionTrain_THR50_2 <- confusionMatrix(Train_predClasses_THR50_2, group_metabric, positive = "1", mode = "everything")
ConfusionTrain_THR50_2

## MCC
MCC_Train_THR25 <- mltools::mcc(pred = Train_predClasses_THR25, actuals = group_metabric)
MCC_Train_THR25

MCC_Train_THR50_1 <- mltools::mcc(pred = Train_predClasses_THR50_1, actuals = group_metabric)
MCC_Train_THR50_1

MCC_Train_THR50_2 <- mltools::mcc(pred = Train_predClasses_THR50_2, actuals = group_metabric)
MCC_Train_THR50_2

#########################################################################
#########################################################################
### Testing

tcga_prob_THR25 <- THR25_model %>% predict(Data_tcga , type = "response")
tcga_prob_THR50_1 <-  THR50_1_model %>% predict(Data_tcga , type = "response")
tcga_prob_THR50_2 <- THR50_2_model %>% predict(Data_tcga , type = "response")


### ROC
ROC_tcga_THR25 <- roc(group_tcga, tcga_prob_THR25, plot = F, print.thres=thr_THR25$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("0", "1"), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_tcga_THR25

ROC_tcga_THR50_1 <- roc(group_tcga, tcga_prob_THR50_1, plot = F, print.thres=thr_THR50_1$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("0", "1"), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_tcga_THR50_1

ROC_tcga_THR50_2 <- roc(group_tcga, tcga_prob_THR50_2, plot = F, print.thres=thr_THR50_2$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("0", "1"), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_tcga_THR50_2

############################
### Get predictions based on best threshold from ROC curve
tcga_predClasses_THR25 <- ifelse(tcga_prob_THR25 >= thr_THR25$threshold, "1", "0")
table(tcga_predClasses_THR25)
tcga_predClasses_THR25 <- factor(tcga_predClasses_THR25, levels = c('0', '1'))

tcga_predClasses_THR50_1 <- ifelse(tcga_prob_THR50_1 >= thr_THR50_1$threshold, "1", "0")
table(tcga_predClasses_THR50_1)
tcga_predClasses_THR50_1 <- factor(tcga_predClasses_THR50_1, levels = c('0', '1'))

tcga_predClasses_THR50_2 <- ifelse(tcga_prob_THR50_2 >= thr_THR50_2$threshold, "1", "0")
table(tcga_predClasses_THR50_2)
tcga_predClasses_THR50_2 <- factor(tcga_predClasses_THR50_2, levels = c('0', '1'))

##################################
### CI  in testing 1
Confusion_tcga_THR25 <- confusionMatrix(tcga_predClasses_THR25, group_tcga, positive = "1", mode = "everything")
Confusion_tcga_THR25

Confusion_tcga_THR50_1 <- confusionMatrix(tcga_predClasses_THR50_1, group_tcga, positive = "1", mode = "everything")
Confusion_tcga_THR50_1

Confusion_tcga_THR50_2 <- confusionMatrix(tcga_predClasses_THR50_2, group_tcga, positive = "1", mode = "everything")
Confusion_tcga_THR50_2


################
## MCC
MCC_tcga_THR25 <- mltools::mcc(pred = tcga_predClasses_THR25, actuals = group_tcga)
MCC_tcga_THR25

MCC_tcga_THR50_1 <- mltools::mcc(pred = tcga_predClasses_THR50_1, actuals = group_tcga)
MCC_tcga_THR50_1

MCC_tcga_THR50_2 <- mltools::mcc(pred = tcga_predClasses_THR50_2, actuals = group_tcga)
MCC_tcga_THR50_2


##########################
## Keep only the relevant information (Metastasis Event and Time)
Phenotype_metabric <- cbind(Pheno_metabric[, c("Overall.Survival.Status", "Overall.Survival..Months.")], 
                                  Train_prob_THR25, Train_prob_THR50_1, Train_prob_THR50_2, Train_predClasses_THR25, Train_predClasses_THR50_1, Train_predClasses_THR50_2)

Phenotype_tcga <- cbind(Pheno_tcga[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Progression.Free.Status", "Progress.Free.Survival..Months.")], 
                              tcga_prob_THR25, tcga_prob_THR50_1, tcga_prob_THR50_2, tcga_predClasses_THR25, tcga_predClasses_THR50_1, tcga_predClasses_THR50_2)

#Expr_metabric <- Expr_metabric[ClassifierGenes, ]
#Expr_tcga <- Expr_tcga[ClassifierGenes, ]


# create a merged pdata and Z-scores object
CoxData_metabric <- data.frame(Phenotype_metabric)
CoxData_tcga <- data.frame(Phenotype_tcga)

# divide the probabilities into quartiles
CoxData_metabric <- CoxData_metabric %>%
  mutate(metabric_prob_THR25_quartiles = ntile(Train_prob_THR25, 4), 
         metabric_prob_THR50_1_quartiles = ntile(Train_prob_THR50_1, 4), 
         metabric_prob_THR50_2_quartiles = ntile(Train_prob_THR50_2, 4))

CoxData_tcga <- CoxData_tcga %>%
  mutate(tcga_prob_THR25_quartiles = ntile(tcga_prob_THR25, 4), 
         tcga_prob_THR50_1_quartiles = ntile(tcga_prob_THR50_1, 4), 
         tcga_prob_THR50_2_quartiles = ntile(tcga_prob_THR50_2, 4))

#CutPoint <- surv_cutpoint(data = CoxData, time = "Time", event = "Event", variables = "ResidualDisease_Score")
#CutPoint

#SurvData <- surv_categorize(CutPoint)

#SurvData$ResidualDisease_Score <- factor(SurvData$ResidualDisease_Score, levels = c("low", "high"))


########################################################################  
## Fit survival curves

# Metabric: 

# init a list for classifier pairs
# THR_25_fil_list <- as.list(unlist(THR_25_fil))
# THR_50_1_fil_list <- as.list(unlist(THR_50_1_fil))
# THR_50_2_fil_list <- as.list(unlist(THR_50_2_fil))
# 
# 
# names(THR_25_fil_list) <- THR_25_fil
# names(THR_50_1_fil_list) <- THR_50_1_fil
# names(THR_50_2_fil_list) <- THR_50_2_fil
# 
# surv_func_metabric_os <- function(x){
#   f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
#   return(surv_fit(f, data = CoxData_metabric))
# }
# 
# fit_list_metabric_os <- lapply(THR_25_fil_list, surv_func_metabric_os)
# 
# # calculate the pvalue
# Pval_list_metabric_os <- surv_pvalue(fit_list_metabric_os)
# 
# Pval_df_metabric_os <- do.call(rbind.data.frame, Pval_list_metabric_os)

############
## Plot survival curves

#plot_list_metabric_os <- ggsurvplot_list(fit_list_metabric_os, CoxData_metabric, legend.title = names(fit_list_metabric_os), pval = TRUE)


#Splot_metabric_os <- arrange_ggsurvplots(plot_list_metabric_os, title = "Survival plots using the 10 pairs individually (THR25)", ncol = 2, nrow = 5)
#ggsave("./figures/sep28/THR25_Metabric_os.pdf", Splot_metabric_os, width = 40, height = 40, units = "cm")


## metabric all genes
Fit_sig_metabric_THR25 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_THR25, data = CoxData_metabric)
Fit_sig_metabric_THR50_1 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_THR50_1, data = CoxData_metabric)
Fit_sig_metabric_THR50_2 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_THR50_2, data = CoxData_metabric)

## by quartiles
Fit_sig_metabric_THR25_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR25_quartiles, data = CoxData_metabric)
Fit_sig_metabric_THR50_1_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR50_1_quartiles, data = CoxData_metabric)
Fit_sig_metabric_THR50_2_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_THR50_2_quartiles, data = CoxData_metabric)

######
surv_pvalue(Fit_sig_metabric_THR25)
surv_pvalue(Fit_sig_metabric_THR50_1)
surv_pvalue(Fit_sig_metabric_THR50_2)

pdf("./figures/logreg/THR25_metabric_allpairs.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_THR25,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 25 (logistic regression) and METABRIC OS')
dev.off()

pdf("./figures/logreg/THR50_1_metabric_allpairs.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_THR50_1,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and METABRIC OS')
dev.off()

pdf("./figures/logreg/THR50_2_metabric_allpairs.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_THR50_2,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_2 (logistic regression) and METABRIC OS')
dev.off()

########
# by quartiles
pdf("./figures/logreg/THR25_metabric_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_THR25_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 25 (logistic regression) and METABRIC OS: quartiles')
dev.off()

pdf("./figures/logreg/THR50_1_metabric_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_THR50_1_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles')
dev.off()

pdf("./figures/logreg/THR50_2_metabric_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_THR50_2_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_2 (logistic regression) and METABRIC OS: quartiles')
dev.off()

#############
## fit coxph model:

# by gene
# surv_func_metabric_os_coxph <- function(x){
#   f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
#   return(coxph(f, data = CoxData_metabric))
# }
# 
# fit_list_metabric_os_coxph_THR25 <- lapply(pairs_list_THR25, surv_func_metabric_os_coxph_THR25)
# names(fit_list_metabric_os_coxph_THR25) <- pairs_list_THR25
# 
# summary_list_metabric_os_coxph_THR25 <- lapply(fit_list_metabric_os_coxph_THR25, summary)
# 
# # get the HR
# HR_list_metabric_os_coxph_THR25 <- lapply(summary_list_metabric_os_coxph_THR25, function(x){
#   HR <- x$conf.int[, 'exp(coef)']
#   Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
#   Pvalue_logrank_test <- x$sctest['pvalue']
#   Pvalue_wald_test <- x$waldtest['pvalue']
#   data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
#              Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
# })


# HR_df_metabric_os_coxph_THR25 <- as.data.frame(do.call(rbind, HR_list_metabric_os_coxph_THR25))
# HR_df_metabric_os_coxph_THR25$variable <- rownames(HR_df_metabric_os_coxph_THR25)
# HR_df_metabric_os_coxph_THR25 <- HR_df_metabric_os_coxph_THR25[order(HR_df_metabric_os_coxph_THR25$HR, decreasing = T), ]

# save the results
#write.csv(HR_df_metabric_os_coxph_THR25, 'objs/sep28/THR25_HR_df_metabric_os_coxph.csv')


########
# by predictions
Fit_sig_metabric_coxph_THR25 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_prob_THR25, data = CoxData_metabric)
summary(Fit_sig_metabric_coxph_THR25)

Fit_sig_metabric_coxph_THR50_1 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_prob_THR50_1, data = CoxData_metabric)
summary(Fit_sig_metabric_coxph_THR50_1)

Fit_sig_metabric_coxph_THR50_2 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_prob_THR50_2, data = CoxData_metabric)
summary(Fit_sig_metabric_coxph_THR50_2)

png('./figures/logreg/THR25_HR_metabric_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_metabric_coxph_THR25, fontsize = 0.5)
dev.off()

png('./figures/logreg/THR50_1_HR_metabric_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_metabric_coxph_THR50_1, fontsize = 0.5)
dev.off()

png('./figures/logreg/THR50_2_HR_metabric_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_metabric_coxph_THR50_2, fontsize = 0.5)
dev.off()

########################################################################  
########################################################################  
## Fit survival curves: TCGA: 

# fit surv curves: os

# surv_func_TCGA_os_THR25 <- function(x){
#   f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
#   return(surv_fit(f, data = CoxData_tcga_THR25))
# }
# 
# fit_list_tcga_os_THR25 <- lapply(pairs_list_THR25, surv_func_TCGA_os_THR25)
# 
# # calculate the pvalue
# Pval_list_tcga_os_THR25 <- surv_pvalue(fit_list_tcga_os_THR25)
# 
# Pval_df_tcga_os_THR25 <- do.call(rbind.data.frame, Pval_list_tcga_os_THR25)

#Pval_df_fil <- Pval_df[Pval_df$pval < 0.05, ] 

############
## Plot survival curves

# plot_list_tcga_os_THR25 <- ggsurvplot_list(fit_list_tcga_os_THR25, CoxData_tcga_THR25, legend.title = names(fit_list_tcga_os_THR25), pval = TRUE)
# 
# 
# Splot_tcga_os_THR25 <- arrange_ggsurvplots(plot_list_tcga_os_THR25, title = "Overall survival in the TCGA using the 10 pairs individually (THR25)", ncol = 2, nrow = 5)
# ggsave("./figures/sep28/THR25_10TSPs_tcga_os.pdf", Splot_tcga_os_THR25, width = 40, height = 40, units = "cm")


## TCGA all pairs
Fit_sig_TCGA_os_THR25 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_predClasses_THR25, data = CoxData_tcga)
Fit_sig_TCGA_os_THR50_1 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_predClasses_THR50_1, data = CoxData_tcga)
Fit_sig_TCGA_os_THR50_2 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_predClasses_THR50_2, data = CoxData_tcga)

## by quartiles
Fit_sig_TCGA_os_THR25_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_THR25_quartiles, data = CoxData_tcga)
Fit_sig_TCGA_os_THR50_1_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_THR50_1_quartiles, data = CoxData_tcga)
Fit_sig_TCGA_os_THR50_2_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_THR50_2_quartiles, data = CoxData_tcga)

pdf("./figures/logreg/THR25_Allpairs_TCGA_os.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_THR25,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 25 (logistic regression) and TCGA OS')
dev.off()


pdf("./figures/logreg/THR50_1_Allpairs_TCGA_os.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_THR50_1,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and TCGA OS')
dev.off()


pdf("./figures/logreg/THR50_2_Allpairs_TCGA_os.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_THR50_1,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_2 (logistic regression) and TCGA OS')
dev.off()

########
# by quartiles
pdf("./figures/logreg/THR25_tcga_os_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_THR25_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 25 (logistic regression) and TCGA OS: quartiles')
dev.off()

pdf("./figures/logreg/THR50_1_tcga_os_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_THR50_1_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and TCGA OS: quartiles')
dev.off()

pdf("./figures/logreg/THR50_2_tcga_os_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_THR50_2_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_2 (logistic regression) and TCGA OS: quartiles')
dev.off()


#############
## fit coxph model:

# by pair
# surv_func_TCGA_os_coxph_THR25 <- function(x){
#   f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
#   return(coxph(f, data = CoxData_tcga_THR25))
# }
# 
# fit_list_TCGA_os_coxph_THR25 <- lapply(pairs_list_THR25, surv_func_TCGA_os_coxph_THR25)
# names(fit_list_TCGA_os_coxph_THR25) <- pairs_list_THR25
# 
# summary_list_TCGA_os_coxph_THR25 <- lapply(fit_list_TCGA_os_coxph_THR25, summary)
# 
# # get the HR
# HR_list_TCGA_os_coxph_THR25 <- lapply(summary_list_TCGA_os_coxph_THR25, function(x){
#   HR <- x$conf.int[, 'exp(coef)']
#   Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
#   Pvalue_logrank_test <- x$sctest['pvalue']
#   Pvalue_wald_test <- x$waldtest['pvalue']
#   data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
#              Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
# })
# 
# 
# HR_df_TCGA_os_coxph_THR25 <- as.data.frame(do.call(rbind, HR_list_TCGA_os_coxph_THR25))
# HR_df_TCGA_os_coxph_THR25$variable <- rownames(HR_df_TCGA_os_coxph_THR25)
# HR_df_TCGA_os_coxph_THR25 <- HR_df_TCGA_os_coxph_THR25[order(HR_df_TCGA_os_coxph_THR25$HR, decreasing = T), ]
# 
# # save the results
# write.csv(HR_df_TCGA_os_coxph_THR25, 'objs/sep28/THR25_HR_df_TCGA_os_coxph.csv')


########
# by predictions
Fit_sig_TCGA_coxph_THR25 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_THR25, data = CoxData_tcga)
summary(Fit_sig_TCGA_coxph_THR25)

Fit_sig_TCGA_coxph_THR50_1 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_THR50_1, data = CoxData_tcga)
summary(Fit_sig_TCGA_coxph_THR50_1)

Fit_sig_TCGA_coxph_THR50_2 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ tcga_prob_THR50_2, data = CoxData_tcga)
summary(Fit_sig_TCGA_coxph_THR50_2)


png('./figures/logreg/THR25_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_coxph_THR25, fontsize = 0.5)
dev.off()

png('./figures/logreg/THR50_1_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_coxph_THR50_1, fontsize = 0.5)
dev.off()

png('./figures/logreg/THR50_2_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_coxph_THR50_2, fontsize = 0.5)
dev.off()
########################################################################  
########################################################################  
## Fit survival curves: TCGA: 

# fit surv curves: PFS

# surv_func_TCGA_pfs_THR25 <- function(x){
#   f <- as.formula(paste("Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~", x))
#   return(surv_fit(f, data = CoxData_tcga_THR25))
# }
# 
# fit_list_tcga_pfs_THR25 <- lapply(pairs_list_THR25, surv_func_TCGA_pfs_THR25)
# 
# # calculate the pvalue
# Pval_list_tcga_pfs_THR25 <- surv_pvalue(fit_list_tcga_pfs_THR25)
# 
# Pval_df_tcga_pfs_THR25 <- do.call(rbind.data.frame, Pval_list_tcga_pfs_THR25)

#Pval_df_fil <- Pval_df[Pval_df$pval < 0.05, ] 

############
## Plot survival curves

# plot_list_tcga_pfs_THR25 <- ggsurvplot_list(fit_list_tcga_pfs_THR25, CoxData_tcga_THR25, legend.title = names(fit_list_tcga_pfs_THR25), pval = TRUE)
# 
# 
# Splot_tcga_pfs_THR25 <- arrange_ggsurvplots(plot_list_tcga_pfs_THR25, title = "Progression free survival in the TCGA using the 10 pairs individually (THR25)", ncol = 2, nrow = 5)
# ggsave("./figures/sep28/10TSPs_tcga_pfs_THR25.pdf", Splot_tcga_pfs_THR25, width = 40, height = 40, units = "cm")


## TCGA all genes
Fit_sig_TCGA_pfs_THR25 <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_predClasses_THR25, data = CoxData_tcga)
Fit_sig_TCGA_pfs_THR50_1 <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_predClasses_THR50_1, data = CoxData_tcga)
Fit_sig_TCGA_pfs_THR50_2 <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_predClasses_THR50_2, data = CoxData_tcga)

## by quartiles
Fit_sig_TCGA_pfs_THR25_quartiles <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_THR25_quartiles, data = CoxData_tcga)
Fit_sig_TCGA_pfs_THR50_1_quartiles <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_THR50_1_quartiles, data = CoxData_tcga)
Fit_sig_TCGA_pfs_THR50_2_quartiles <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_THR50_2_quartiles, data = CoxData_tcga)

pdf("./figures/logreg/THR25_Allpairs_TCGA_PFS.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_THR25,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 25 (logistic regression) and TCGA PFS')
dev.off()

pdf("./figures/logreg/THR50_1_Allpairs_TCGA_PFS.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_THR50_1,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and TCGA PFS')
dev.off()

pdf("./figures/logreg/THR50_2_Allpairs_TCGA_PFS.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_THR50_2,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_2 (logistic regression) and TCGA PFS')
dev.off()

########
# by quartiles
pdf("./figures/logreg/THR25_tcga_pfs_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_THR25_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 25 (logistic regression) and TCGA PFS: quartiles')
dev.off()

pdf("./figures/logreg/THR50_1_tcga_pfs_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_THR50_1_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_1 (logistic regression) and TCGA PFS: quartiles')
dev.off()

pdf("./figures/logreg/THR50_2_tcga_pfs_quartiles.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_THR50_2_quartiles,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'THR 50_2 (logistic regression) and TCGA PFS: quartiles')
dev.off()
#############
## fit coxph model:

# by gene
# surv_func_TCGA_pfs_coxph_THR25 <- function(x){
#   f <- as.formula(paste("Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~", x))
#   return(coxph(f, data = CoxData_tcga_THR25))
# }
# 
# fit_list_TCGA_pfs_coxph_THR25 <- lapply(pairs_list_THR25, surv_func_TCGA_pfs_coxph_THR25)
# names(fit_list_TCGA_pfs_coxph_THR25) <- pairs_list_THR25
# 
# summary_list_TCGA_pfs_coxph_THR25 <- lapply(fit_list_TCGA_pfs_coxph_THR25, summary)
# 
# # get the HR
# HR_list_TCGA_pfs_coxph_THR25 <- lapply(summary_list_TCGA_pfs_coxph_THR25, function(x){
#   HR <- x$conf.int[, 'exp(coef)']
#   Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
#   Pvalue_logrank_test <- x$sctest['pvalue']
#   Pvalue_wald_test <- x$waldtest['pvalue']
#   data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
#              Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
# })
# 
# 
# HR_df_TCGA_pfs_coxph_THR25 <- as.data.frame(do.call(rbind, HR_list_TCGA_pfs_coxph_THR25))
# HR_df_TCGA_pfs_coxph_THR25$variable <- rownames(HR_df_TCGA_pfs_coxph_THR25)
# HR_df_TCGA_pfs_coxph_THR25 <- HR_df_TCGA_pfs_coxph_THR25[order(HR_df_TCGA_pfs_coxph_THR25$HR, decreasing = T), ]
# 
# # save the results
# write.csv(HR_df_TCGA_pfs_coxph_THR25, 'objs/sep28/HR_df_TCGA_pfs_coxph_THR25.csv')


###########
# by predictions
Fit_sig_TCGA_pfs_coxph_THR25 <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_THR25, data = CoxData_tcga)
summary(Fit_sig_TCGA_pfs_coxph_THR25)

Fit_sig_TCGA_pfs_coxph_THR50_1 <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_THR50_1, data = CoxData_tcga)
summary(Fit_sig_TCGA_pfs_coxph_THR50_1)

Fit_sig_TCGA_pfs_coxph_THR50_2 <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ tcga_prob_THR50_2, data = CoxData_tcga)
summary(Fit_sig_TCGA_pfs_coxph_THR50_2)

png('./figures/logreg/THR25_HR_tcga_pfs.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_pfs_coxph_THR25, fontsize = 0.5)
dev.off()

png('./figures/logreg/THR50_1_HR_tcga_pfs.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_pfs_coxph_THR50_1, fontsize = 0.5)
dev.off()

png('./figures/logreg/THR50_2_HR_tcga_pfs.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_pfs_coxph_THR50_2, fontsize = 0.5)
dev.off()




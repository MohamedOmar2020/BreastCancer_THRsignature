###### 
# Clean Work space
rm(list = ls())

############################################################################
### Load library
require(switchBox)
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
library(readxl)
library(survival)
library(survminer)

#################
THR_signature <- readxl::read_xlsx("./data/THR Signatures_sep23.xlsx")

# get THR 25 and 50
THR_25 <- THR_signature$`THR-25`[!is.na(THR_signature$`THR-25`)]
THR_50_1 <- THR_signature$`THR-50.1`[!is.na(THR_signature$`THR-50.1`)]
THR_50_2 <- THR_signature$`THR-50.2`[!is.na(THR_signature$`THR-50.2`)]


# make TSPs
#myTSPs <- t(combn(THR_signature$Gene, 2))
myTSPs_25 <- t(combn(THR_25, 2))
colnames(myTSPs_25) <- c("gene1", "gene2")

myTSPs_50_1 <- t(combn(THR_50_1, 2))
colnames(myTSPs_50_1) <- c("gene1", "gene2")

myTSPs_50_2 <- t(combn(THR_50_2, 2))
colnames(myTSPs_50_2) <- c("gene1", "gene2")

################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')

###################################

# Missing genes: TBC1D30 (absent in both) and TNFAIP2 (absent in metabric)
rownames(Expr_metabric)[grep('^TBC1D30', rownames(Expr_metabric))]
rownames(Expr_tcga)[grep('^TBC1D30', rownames(Expr_tcga))]

rownames(Expr_metabric)[grep('^TNFAIP2', rownames(Expr_metabric))]
rownames(Expr_tcga)[grep('^TNFAIP2', rownames(Expr_tcga))]

### Common genes
keepGns_datasets <- intersect(rownames(Expr_tcga), rownames(Expr_metabric))
keepGns_tsp_THR25 <- intersect(as.vector(myTSPs_25), keepGns_datasets)
keepGns_tsp_THR50_1 <- intersect(as.vector(myTSPs_50_1), keepGns_datasets)
keepGns_tsp_THR50_2 <- intersect(as.vector(myTSPs_50_2), keepGns_datasets)

#Expr_metabric <- Expr_metabric[keepGns_tsp, ]
#Expr_tcga <- Expr_tcga[keepGns_tsp, ]

### For the TSP
myTSPs_25 <- myTSPs_25[myTSPs_25[,1] %in% keepGns_tsp_THR25 & myTSPs_25[,2] %in% keepGns_tsp_THR25 , ]
myTSPs_50_1 <- myTSPs_50_1[myTSPs_50_1[,1] %in% keepGns_tsp_THR50_1 & myTSPs_50_1[,2] %in% keepGns_tsp_THR50_1 , ]
myTSPs_50_2 <- myTSPs_50_2[myTSPs_50_2[,1] %in% keepGns_tsp_THR50_2 & myTSPs_50_2[,2] %in% keepGns_tsp_THR50_2 , ]

###########################################################################
### TRAINING using restricted pairs
###########################################################################
### Set Feature number and max k
ktsp <- c(3:25) #7
featNo <- nrow(Expr_metabric)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333) # 9/7 # 25

ktsp_metabric_THR25 <- SWAP.Train.KTSP(
  Expr_metabric, group_metabric, krange=ktsp, disjoint = T, 
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs_25)

ktsp_metabric_THR50_1 <- SWAP.Train.KTSP(
  Expr_metabric, group_metabric, krange=ktsp, disjoint = T, 
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs_50_1)

ktsp_metabric_THR50_2 <- SWAP.Train.KTSP(
  Expr_metabric, group_metabric, krange=ktsp, disjoint = T, 
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs_50_2)

ktsp_metabric_THR25
ktsp_metabric_THR50_1
ktsp_metabric_THR50_2

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStats_metabric_THR25 <- SWAP.KTSP.Statistics(inputMat = Expr_metabric, classifier = ktsp_metabric_THR25, CombineFunc = sum)
ktspStats_metabric_THR50_1 <- SWAP.KTSP.Statistics(inputMat = Expr_metabric, classifier = ktsp_metabric_THR50_1, CombineFunc = sum)
ktspStats_metabric_THR50_2 <- SWAP.KTSP.Statistics(inputMat = Expr_metabric, classifier = ktsp_metabric_THR50_2, CombineFunc = sum)

### Threshold
thr_THR25 <- coords(roc(group_metabric, ktspStats_metabric_THR25$statistics, levels = c(0, 1), direction = "<" ), "best")["threshold"]
thr_THR50_1 <- coords(roc(group_metabric, ktspStats_metabric_THR50_1$statistics, levels = c(0, 1), direction = "<" ), "best")["threshold"]
thr_THR50_2 <- coords(roc(group_metabric, ktspStats_metabric_THR50_2$statistics, levels = c(0, 1), direction = "<" ), "best")["threshold"]


### Print ROC curve local maximas
coords(roc(group_metabric, ktspStats_metabric_THR25$statistics, levels = c(0, 1), direction = "<" ), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROC_metabric_THR25 <- roc(group_metabric, ktspStats_metabric_THR25$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_metabric_THR25

ROC_metabric_THR50_1 <- roc(group_metabric, ktspStats_metabric_THR50_1$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_metabric_THR50_1

ROC_metabric_THR50_2 <- roc(group_metabric, ktspStats_metabric_THR50_2$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_metabric_THR50_2

### Get predictions based on best threshold from ROC curve
prediction_metabric_THR25 <- SWAP.KTSP.Classify(Expr_metabric, ktsp_metabric_THR25, DecisionFunc = function(x) sum(x) > thr_THR25)
prediction_metabric_THR50_1 <- SWAP.KTSP.Classify(Expr_metabric, ktsp_metabric_THR50_1, DecisionFunc = function(x) sum(x) > thr_THR50_1)
prediction_metabric_THR50_2 <- SWAP.KTSP.Classify(Expr_metabric, ktsp_metabric_THR50_2, DecisionFunc = function(x) sum(x) > thr_THR50_2)

### Resubstitution performance in the TRAINING set
confusionMatrix(prediction_metabric_THR25, group_metabric, positive = '1', mode = "everything")
confusionMatrix(prediction_metabric_THR50_1, group_metabric, positive = '1', mode = "everything")
confusionMatrix(prediction_metabric_THR50_2, group_metabric, positive = '1', mode = "everything")

mltools::mcc(pred = prediction_metabric_THR25, actuals = group_metabric)
mltools::mcc(pred = prediction_metabric_THR50_1, actuals = group_metabric)
mltools::mcc(pred = prediction_metabric_THR50_2, actuals = group_metabric)

########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStats_tcga_THR25 <- SWAP.KTSP.Statistics(inputMat = Expr_tcga, classifier = ktsp_metabric_THR25, CombineFunc = sum)
ktspStats_tcga_THR50_1 <- SWAP.KTSP.Statistics(inputMat = Expr_tcga, classifier = ktsp_metabric_THR50_1, CombineFunc = sum)
ktspStats_tcga_THR50_2 <- SWAP.KTSP.Statistics(inputMat = Expr_tcga, classifier = ktsp_metabric_THR50_2, CombineFunc = sum)

## Plot curve
ROC_tcga_THR25 <- roc(group_tcga, ktspStats_tcga_THR25$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")
ROC_tcga_THR25

ROC_tcga_THR50_1 <- roc(group_tcga, ktspStats_tcga_THR50_1$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")
ROC_tcga_THR50_1

ROC_tcga_THR50_2 <- roc(group_tcga, ktspStats_tcga_THR50_2$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")
ROC_tcga_THR50_2

### Get predictions based on best threshold from ROC curve
prediction_tcga_THR25 <- SWAP.KTSP.Classify(Expr_tcga, ktsp_metabric_THR25, DecisionFunc = function(x) sum(x) > thr_THR25)
prediction_tcga_THR50_1 <- SWAP.KTSP.Classify(Expr_tcga, ktsp_metabric_THR50_1, DecisionFunc = function(x) sum(x) > thr_THR50_1)
prediction_tcga_THR50_2 <- SWAP.KTSP.Classify(Expr_tcga, ktsp_metabric_THR50_2, DecisionFunc = function(x) sum(x) > thr_THR50_2)

### Resubstitution performance in the Test set
confusionMatrix(prediction_tcga_THR25, group_tcga, positive = "1", mode = "everything")
confusionMatrix(prediction_tcga_THR50_1, group_tcga, positive = "1", mode = "everything")
confusionMatrix(prediction_tcga_THR50_2, group_tcga, positive = "1", mode = "everything")


mltools::mcc(pred = prediction_tcga_THR25, actuals = group_tcga)
mltools::mcc(pred = prediction_tcga_THR50_1, actuals = group_tcga)
mltools::mcc(pred = prediction_tcga_THR50_2, actuals = group_tcga)


############################################################
############################################################
############################################################
# test the ktsp pairs with survival
ClassifierGenes_THR25 <- as.vector(ktsp_metabric_THR25$TSPs)


##########################
## Keep only the relevant information (Metastasis Event and Time)
Phenotype_metabric <- cbind(Pheno_metabric[, c("Overall.Survival.Status", "Overall.Survival..Months.")], 
                            ktspStats_metabric_THR25$comparisons, ktspStats_metabric_THR50_1$comparisons, ktspStats_metabric_THR50_2$comparisons,  
                            prediction_metabric_THR25, prediction_metabric_THR50_1, prediction_metabric_THR50_2)

Phenotype_tcga <- cbind(Pheno_tcga[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Progression.Free.Status", "Progress.Free.Survival..Months.")], 
                        ktspStats_tcga_THR25$comparisons, ktspStats_tcga_THR50_1$comparisons, ktspStats_tcga_THR50_2$comparisons, 
                        prediction_tcga_THR25, prediction_tcga_THR50_1, prediction_tcga_THR50_2)

# create a merged pdata and Z-scores object
CoxData_metabric <- data.frame(Phenotype_metabric)
CoxData_tcga <- data.frame(Phenotype_tcga)

#CutPoint <- surv_cutpoint(data = CoxData, time = "Time", event = "Event", variables = "ResidualDisease_Score")
#CutPoint

#SurvData <- surv_categorize(CutPoint)

#SurvData$ResidualDisease_Score <- factor(SurvData$ResidualDisease_Score, levels = c("low", "high"))


########################################################################  
## Fit survival curves

# Metabric: 

# init a list for classifier pairs
pairs_list_THR25 <- list()
pairs_list_THR50_1 <- list()
pairs_list_THR50_2 <- list()

for(i in seq(1,nrow(ktsp_metabric_THR25$TSPs))){
  pairs_list_THR25[i] <- paste0(ktsp_metabric_THR25$TSPs[i,1], '.', ktsp_metabric_THR25$TSPs[i,2])
}
names(pairs_list_THR25) <- pairs_list_THR25

for(i in seq(1,nrow(ktsp_metabric_THR50_1$TSPs))){
  pairs_list_THR50_1[i] <- paste0(ktsp_metabric_THR50_1$TSPs[i,1], '.', ktsp_metabric_THR50_1$TSPs[i,2])
}
names(pairs_list_THR50_1) <- pairs_list_THR50_1

for(i in seq(1,nrow(ktsp_metabric_THR50_2$TSPs))){
  pairs_list_THR50_2[i] <- paste0(ktsp_metabric_THR50_2$TSPs[i,1], '.', ktsp_metabric_THR50_2$TSPs[i,2])
}
names(pairs_list_THR50_2) <- pairs_list_THR50_2

######
# fit surv curves

surv_func_metabric_os <- function(x){
  f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
  return(surv_fit(f, data = CoxData_metabric))
}

fit_list_metabric_os_THR25 <- lapply(pairs_list_THR25, surv_func_metabric_os)
fit_list_metabric_os_THR50_1 <- lapply(pairs_list_THR50_1, surv_func_metabric_os)
fit_list_metabric_os_THR50_2 <- lapply(pairs_list_THR50_2, surv_func_metabric_os)

# calculate the pvalue
Pval_list_metabric_os_THR25 <- surv_pvalue(fit_list_metabric_os_THR25)
Pval_df_metabric_os_THR25 <- do.call(rbind.data.frame, Pval_list_metabric_os_THR25)

Pval_list_metabric_os_THR50_1 <- surv_pvalue(fit_list_metabric_os_THR50_1)
Pval_df_metabric_os_THR50_1 <- do.call(rbind.data.frame, Pval_list_metabric_os_THR50_1)

Pval_list_metabric_os_THR50_2 <- surv_pvalue(fit_list_metabric_os_THR50_2)
Pval_df_metabric_os_THR50_2 <- do.call(rbind.data.frame, Pval_list_metabric_os_THR50_2)

############
## Plot survival curves
plot_list_metabric_os_THR25 <- ggsurvplot_list(fit_list_metabric_os_THR25, CoxData_metabric, legend.title = names(fit_list_metabric_os_THR25), pval = TRUE)
plot_list_metabric_os_THR50_1 <- ggsurvplot_list(fit_list_metabric_os_THR50_1, CoxData_metabric, legend.title = names(fit_list_metabric_os_THR50_1), pval = TRUE)
plot_list_metabric_os_THR50_2 <- ggsurvplot_list(fit_list_metabric_os_THR50_2, CoxData_metabric, legend.title = names(fit_list_metabric_os_THR50_2), pval = TRUE)


Splot_metabric_os_THR25 <- arrange_ggsurvplots(plot_list_metabric_os_THR25, title = "Survival plots using the pairs individually (THR25)", ncol = 2, nrow = 6)
ggsave("./figures/sep28/THR25_Metabric_os.pdf", Splot_metabric_os_THR25, width = 40, height = 40, units = "cm")

Splot_metabric_os_THR50_1 <- arrange_ggsurvplots(plot_list_metabric_os_THR50_1, title = "Survival plots using the pairs individually (THR50_1)", ncol = 2, nrow = 5)
ggsave("./figures/sep28/THR50_1_Metabric_os.pdf", Splot_metabric_os_THR50_1, width = 40, height = 40, units = "cm")

Splot_metabric_os_THR50_2 <- arrange_ggsurvplots(plot_list_metabric_os_THR50_2, title = "Survival plots using the pairs individually (THR50_2)", ncol = 3, nrow = 5)
ggsave("./figures/sep28/THR50_2_Metabric_os.pdf", Splot_metabric_os_THR50_2, width = 40, height = 40, units = "cm")

## metabric all pairs
Fit_sig_metabric_THR25 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_metabric_THR25, data = CoxData_metabric)
surv_pvalue(Fit_sig_metabric_THR25)

Fit_sig_metabric_THR50_1 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_metabric_THR50_1, data = CoxData_metabric)
surv_pvalue(Fit_sig_metabric_THR50_1)

Fit_sig_metabric_THR50_2 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_metabric_THR50_2, data = CoxData_metabric)
surv_pvalue(Fit_sig_metabric_THR50_2)

#####
pdf("./figures/sep28/THR25_metabric_allpairs.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_THR25,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'k-TSPs and METABRIC OS (THR25)')
dev.off()

pdf("./figures/sep28/THR50_1_metabric_allpairs.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_THR50_1,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'k-TSPs and METABRIC OS (THR50_1)')
dev.off()

pdf("./figures/sep28/THR50_2_metabric_allpairs.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric_THR50_2,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'k-TSPs and METABRIC OS (THR50_2)')
dev.off()

#############
## fit coxph model:

# by pair
surv_func_metabric_os_coxph <- function(x){
  f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
  return(coxph(f, data = CoxData_metabric))
}

fit_list_metabric_os_coxph_THR25 <- lapply(pairs_list_THR25, surv_func_metabric_os_coxph)
names(fit_list_metabric_os_coxph_THR25) <- pairs_list_THR25
summary_list_metabric_os_coxph_THR25 <- lapply(fit_list_metabric_os_coxph_THR25, summary)

fit_list_metabric_os_coxph_THR50_1 <- lapply(pairs_list_THR50_1, surv_func_metabric_os_coxph)
names(fit_list_metabric_os_coxph_THR50_1) <- pairs_list_THR50_1
summary_list_metabric_os_coxph_THR50_1 <- lapply(fit_list_metabric_os_coxph_THR50_1, summary)

fit_list_metabric_os_coxph_THR50_2 <- lapply(pairs_list_THR50_2, surv_func_metabric_os_coxph)
names(fit_list_metabric_os_coxph_THR50_2) <- pairs_list_THR50_2
summary_list_metabric_os_coxph_THR50_2 <- lapply(fit_list_metabric_os_coxph_THR50_2, summary)

# get the HR
HR_list_metabric_os_coxph_THR25 <- lapply(summary_list_metabric_os_coxph_THR25, function(x){
  HR <- x$conf.int[, 'exp(coef)']
  Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
  Pvalue_logrank_test <- x$sctest['pvalue']
  Pvalue_wald_test <- x$waldtest['pvalue']
  data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
             Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
})


HR_df_metabric_os_coxph_THR25 <- as.data.frame(do.call(rbind, HR_list_metabric_os_coxph_THR25))
HR_df_metabric_os_coxph_THR25$variable <- rownames(HR_df_metabric_os_coxph_THR25)
HR_df_metabric_os_coxph_THR25 <- HR_df_metabric_os_coxph_THR25[order(HR_df_metabric_os_coxph_THR25$HR, decreasing = T), ]

#######
HR_list_metabric_os_coxph_THR50_1 <- lapply(summary_list_metabric_os_coxph_THR50_1, function(x){
  HR <- x$conf.int[, 'exp(coef)']
  Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
  Pvalue_logrank_test <- x$sctest['pvalue']
  Pvalue_wald_test <- x$waldtest['pvalue']
  data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
             Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
})


HR_df_metabric_os_coxph_THR50_1 <- as.data.frame(do.call(rbind, HR_list_metabric_os_coxph_THR50_1))
HR_df_metabric_os_coxph_THR50_1$variable <- rownames(HR_df_metabric_os_coxph_THR50_1)
HR_df_metabric_os_coxph_THR50_1 <- HR_df_metabric_os_coxph_THR50_1[order(HR_df_metabric_os_coxph_THR50_1$HR, decreasing = T), ]

######
HR_list_metabric_os_coxph_THR50_2 <- lapply(summary_list_metabric_os_coxph_THR50_2, function(x){
  HR <- x$conf.int[, 'exp(coef)']
  Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
  Pvalue_logrank_test <- x$sctest['pvalue']
  Pvalue_wald_test <- x$waldtest['pvalue']
  data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
             Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
})


HR_df_metabric_os_coxph_THR50_2 <- as.data.frame(do.call(rbind, HR_list_metabric_os_coxph_THR50_2))
HR_df_metabric_os_coxph_THR50_2$variable <- rownames(HR_df_metabric_os_coxph_THR50_2)
HR_df_metabric_os_coxph_THR50_2 <- HR_df_metabric_os_coxph_THR50_2[order(HR_df_metabric_os_coxph_THR50_2$HR, decreasing = T), ]

######
# save the results
write.csv(HR_df_metabric_os_coxph_THR25, 'objs/sep28/THR25_HR_df_metabric_os_coxph.csv')
write.csv(HR_df_metabric_os_coxph_THR50_1, 'objs/sep28/THR50_1_HR_df_metabric_os_coxph.csv')
write.csv(HR_df_metabric_os_coxph_THR50_2, 'objs/sep28/THR50_2_HR_df_metabric_os_coxph.csv')


########
# by predictions
Fit_sig_metabric_coxph_THR25 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_metabric_THR25, data = CoxData_metabric)
summary(Fit_sig_metabric_coxph_THR25)

Fit_sig_metabric_coxph_THR50_1 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_metabric_THR50_1, data = CoxData_metabric)
summary(Fit_sig_metabric_coxph_THR50_1)

Fit_sig_metabric_coxph_THR50_2 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_metabric_THR50_2, data = CoxData_metabric)
summary(Fit_sig_metabric_coxph_THR50_2)

png('./figures/sep28/THR25_HR_metabric_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_metabric_coxph_THR25, fontsize = 0.5)
dev.off()

png('./figures/sep28/THR50_1_HR_metabric_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_metabric_coxph_THR50_1, fontsize = 0.5)
dev.off()

png('./figures/sep28/THR50_2_HR_metabric_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_metabric_coxph_THR50_2, fontsize = 0.5)
dev.off()

########################################################################  
########################################################################  
## Fit survival curves: TCGA: 

# fit surv curves: os

surv_func_TCGA_os <- function(x){
  f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
  return(surv_fit(f, data = CoxData_tcga))
}

fit_list_tcga_os_THR25 <- lapply(pairs_list_THR25, surv_func_TCGA_os)
fit_list_tcga_os_THR50_1 <- lapply(pairs_list_THR50_1, surv_func_TCGA_os)
fit_list_tcga_os_THR50_2 <- lapply(pairs_list_THR50_2, surv_func_TCGA_os)

# calculate the pvalue
Pval_list_tcga_os_THR25 <- surv_pvalue(fit_list_tcga_os_THR25)
Pval_df_tcga_os_THR25 <- do.call(rbind.data.frame, Pval_list_tcga_os_THR25)

Pval_list_tcga_os_THR50_1 <- surv_pvalue(fit_list_tcga_os_THR50_1)
Pval_df_tcga_os_THR50_1 <- do.call(rbind.data.frame, Pval_list_tcga_os_THR50_1)

Pval_list_tcga_os_THR50_2 <- surv_pvalue(fit_list_tcga_os_THR50_2)
Pval_df_tcga_os_THR50_2 <- do.call(rbind.data.frame, Pval_list_tcga_os_THR50_2)

############
## Plot survival curves

plot_list_tcga_os_THR25 <- ggsurvplot_list(fit_list_tcga_os_THR25, CoxData_tcga, legend.title = names(fit_list_tcga_os_THR25), pval = TRUE)
plot_list_tcga_os_THR50_1 <- ggsurvplot_list(fit_list_tcga_os_THR50_1, CoxData_tcga, legend.title = names(fit_list_tcga_os_THR50_1), pval = TRUE)
plot_list_tcga_os_THR50_2 <- ggsurvplot_list(fit_list_tcga_os_THR50_2, CoxData_tcga, legend.title = names(fit_list_tcga_os_THR50_2), pval = TRUE)


Splot_tcga_os_THR25 <- arrange_ggsurvplots(plot_list_tcga_os_THR25, title = "Overall survival in the TCGA using the pairs individually (THR25)", ncol = 2, nrow = 6)
ggsave("./figures/sep28/THR25_10TSPs_tcga_os.pdf", Splot_tcga_os_THR25, width = 40, height = 40, units = "cm")

Splot_tcga_os_THR50_1 <- arrange_ggsurvplots(plot_list_tcga_os_THR50_1, title = "Overall survival in the TCGA using the pairs individually (THR50_1)", ncol = 2, nrow = 5)
ggsave("./figures/sep28/THR50_1_10TSPs_tcga_os.pdf", Splot_tcga_os_THR50_1, width = 40, height = 40, units = "cm")

Splot_tcga_os_THR50_2 <- arrange_ggsurvplots(plot_list_tcga_os_THR50_2, title = "Overall survival in the TCGA using the pairs individually (THR50_2)", ncol = 3, nrow = 5)
ggsave("./figures/sep28/THR50_2_10TSPs_tcga_os.pdf", Splot_tcga_os_THR50_2, width = 40, height = 40, units = "cm")

## TCGA all pairs
Fit_sig_TCGA_os_THR25 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_tcga_THR25, data = CoxData_tcga)
Fit_sig_TCGA_os_THR50_1 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_tcga_THR50_1, data = CoxData_tcga)
Fit_sig_TCGA_os_THR50_2 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_tcga_THR50_2, data = CoxData_tcga)

pdf("./figures/sep28/THR25_Allpairs_TCGA_os.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_THR25,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'k-TSPs and TCGA OS (THR25)')
dev.off()

########
pdf("./figures/sep28/THR50_1_Allpairs_TCGA_os.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_THR50_1,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'k-TSPs and TCGA OS (THR50_1)')
dev.off()

########
pdf("./figures/sep28/THR50_2_Allpairs_TCGA_os.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os_THR50_2,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'k-TSPs and TCGA OS (THR50_2)')
dev.off()

#############
## fit coxph model:

# by pair
surv_func_TCGA_os_coxph <- function(x){
  f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
  return(coxph(f, data = CoxData_tcga))
}

fit_list_TCGA_os_coxph_THR25 <- lapply(pairs_list_THR25, surv_func_TCGA_os_coxph)
names(fit_list_TCGA_os_coxph_THR25) <- pairs_list_THR25
summary_list_TCGA_os_coxph_THR25 <- lapply(fit_list_TCGA_os_coxph_THR25, summary)

fit_list_TCGA_os_coxph_THR50_1 <- lapply(pairs_list_THR50_1, surv_func_TCGA_os_coxph)
names(fit_list_TCGA_os_coxph_THR50_1) <- pairs_list_THR50_1
summary_list_TCGA_os_coxph_THR50_1 <- lapply(fit_list_TCGA_os_coxph_THR50_1, summary)

fit_list_TCGA_os_coxph_THR50_2 <- lapply(pairs_list_THR50_2, surv_func_TCGA_os_coxph)
names(fit_list_TCGA_os_coxph_THR50_2) <- pairs_list_THR50_2
summary_list_TCGA_os_coxph_THR50_2 <- lapply(fit_list_TCGA_os_coxph_THR50_2, summary)

##########
# get the HR
HR_list_TCGA_os_coxph_THR25 <- lapply(summary_list_TCGA_os_coxph_THR25, function(x){
  HR <- x$conf.int[, 'exp(coef)']
  Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
  Pvalue_logrank_test <- x$sctest['pvalue']
  Pvalue_wald_test <- x$waldtest['pvalue']
  data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
             Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
})


HR_df_TCGA_os_coxph_THR25 <- as.data.frame(do.call(rbind, HR_list_TCGA_os_coxph_THR25))
HR_df_TCGA_os_coxph_THR25$variable <- rownames(HR_df_TCGA_os_coxph_THR25)
HR_df_TCGA_os_coxph_THR25 <- HR_df_TCGA_os_coxph_THR25[order(HR_df_TCGA_os_coxph_THR25$HR, decreasing = T), ]

########
HR_list_TCGA_os_coxph_THR50_1 <- lapply(summary_list_TCGA_os_coxph_THR50_1, function(x){
  HR <- x$conf.int[, 'exp(coef)']
  Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
  Pvalue_logrank_test <- x$sctest['pvalue']
  Pvalue_wald_test <- x$waldtest['pvalue']
  data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
             Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
})


HR_df_TCGA_os_coxph_THR50_1 <- as.data.frame(do.call(rbind, HR_list_TCGA_os_coxph_THR50_1))
HR_df_TCGA_os_coxph_THR50_1$variable <- rownames(HR_df_TCGA_os_coxph_THR50_1)
HR_df_TCGA_os_coxph_THR50_1 <- HR_df_TCGA_os_coxph_THR50_1[order(HR_df_TCGA_os_coxph_THR50_1$HR, decreasing = T), ]

###########
HR_list_TCGA_os_coxph_THR50_2 <- lapply(summary_list_TCGA_os_coxph_THR50_2, function(x){
  HR <- x$conf.int[, 'exp(coef)']
  Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
  Pvalue_logrank_test <- x$sctest['pvalue']
  Pvalue_wald_test <- x$waldtest['pvalue']
  data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
             Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
})


HR_df_TCGA_os_coxph_THR50_2 <- as.data.frame(do.call(rbind, HR_list_TCGA_os_coxph_THR50_2))
HR_df_TCGA_os_coxph_THR50_2$variable <- rownames(HR_df_TCGA_os_coxph_THR50_2)
HR_df_TCGA_os_coxph_THR50_2 <- HR_df_TCGA_os_coxph_THR50_2[order(HR_df_TCGA_os_coxph_THR50_2$HR, decreasing = T), ]

##########################
# save the results
write.csv(HR_df_TCGA_os_coxph_THR25, 'objs/sep28/THR25_HR_df_TCGA_os_coxph.csv')
write.csv(HR_df_TCGA_os_coxph_THR50_1, 'objs/sep28/THR50_1_HR_df_TCGA_os_coxph.csv')
write.csv(HR_df_TCGA_os_coxph_THR50_2, 'objs/sep28/THR50_2_HR_df_TCGA_os_coxph.csv')


########
# by predictions
Fit_sig_TCGA_coxph_THR25 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_tcga_THR25, data = CoxData_tcga)
summary(Fit_sig_TCGA_coxph_THR25)

Fit_sig_TCGA_coxph_THR50_1 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_tcga_THR50_1, data = CoxData_tcga)
summary(Fit_sig_TCGA_coxph_THR50_1)

Fit_sig_TCGA_coxph_THR50_2 <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_tcga_THR50_2, data = CoxData_tcga)
summary(Fit_sig_TCGA_coxph_THR50_2)

png('./figures/sep28/THR25_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_coxph_THR25, fontsize = 0.5)
dev.off()

png('./figures/sep28/THR50_1_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_coxph_THR50_1, fontsize = 0.5)
dev.off()

png('./figures/sep28/THR50_2_HR_tcga_os.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_coxph_THR50_2, fontsize = 0.5)
dev.off()

########################################################################  
########################################################################  
## Fit survival curves: TCGA: 

# fit surv curves: PFS

surv_func_TCGA_pfs <- function(x){
  f <- as.formula(paste("Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~", x))
  return(surv_fit(f, data = CoxData_tcga))
}

fit_list_tcga_pfs_THR25 <- lapply(pairs_list_THR25, surv_func_TCGA_pfs)
fit_list_tcga_pfs_THR50_1 <- lapply(pairs_list_THR50_1, surv_func_TCGA_pfs)
fit_list_tcga_pfs_THR50_2 <- lapply(pairs_list_THR50_2, surv_func_TCGA_pfs)

# calculate the pvalue
Pval_list_tcga_pfs_THR25 <- surv_pvalue(fit_list_tcga_pfs_THR25)
Pval_df_tcga_pfs_THR25 <- do.call(rbind.data.frame, Pval_list_tcga_pfs_THR25)

Pval_list_tcga_pfs_THR50_1 <- surv_pvalue(fit_list_tcga_pfs_THR50_1)
Pval_df_tcga_pfs_THR50_1 <- do.call(rbind.data.frame, Pval_list_tcga_pfs_THR50_1)

Pval_list_tcga_pfs_THR50_2 <- surv_pvalue(fit_list_tcga_pfs_THR50_2)
Pval_list_tcga_pfs_THR50_2 <- do.call(rbind.data.frame, Pval_list_tcga_pfs_THR50_2)

############
## Plot survival curves
plot_list_tcga_pfs_THR25 <- ggsurvplot_list(fit_list_tcga_pfs_THR25, CoxData_tcga, legend.title = names(fit_list_tcga_pfs_THR25), pval = TRUE)
plot_list_tcga_pfs_THR50_1 <- ggsurvplot_list(fit_list_tcga_pfs_THR50_1, CoxData_tcga, legend.title = names(fit_list_tcga_pfs_THR50_1), pval = TRUE)
plot_list_tcga_pfs_THR50_2 <- ggsurvplot_list(fit_list_tcga_pfs_THR50_2, CoxData_tcga, legend.title = names(fit_list_tcga_pfs_THR50_2), pval = TRUE)


Splot_tcga_pfs_THR25 <- arrange_ggsurvplots(plot_list_tcga_pfs_THR25, title = "Progression free survival in the TCGA using the pairs individually (THR25)", ncol = 2, nrow = 6)
ggsave("./figures/sep28/10TSPs_tcga_pfs_THR25.pdf", Splot_tcga_pfs_THR25, width = 40, height = 40, units = "cm")

Splot_tcga_pfs_THR50_1 <- arrange_ggsurvplots(plot_list_tcga_pfs_THR50_1, title = "Progression free survival in the TCGA using the pairs individually (THR50_1)", ncol = 2, nrow = 5)
ggsave("./figures/sep28/10TSPs_tcga_pfs_THR50_1.pdf", Splot_tcga_pfs_THR50_1, width = 40, height = 40, units = "cm")

Splot_tcga_pfs_THR50_2 <- arrange_ggsurvplots(plot_list_tcga_pfs_THR50_2, title = "Progression free survival in the TCGA using the pairs individually (THR50_2)", ncol = 3, nrow = 5)
ggsave("./figures/sep28/10TSPs_tcga_pfs_THR50_2.pdf", Splot_tcga_pfs_THR50_2, width = 40, height = 40, units = "cm")

## TCGA all pairs
Fit_sig_TCGA_pfs_THR25 <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ prediction_tcga_THR25, data = CoxData_tcga)
Fit_sig_TCGA_pfs_THR50_1 <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ prediction_tcga_THR50_1, data = CoxData_tcga)
Fit_sig_TCGA_pfs_THR50_2 <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ prediction_tcga_THR50_2, data = CoxData_tcga)

pdf("./figures/sep28/THR25_Allpairs_TCGA_PFS.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_THR25,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'k-TSPs and TCGA PFS (THR25)')
dev.off()

pdf("./figures/sep28/THR50_1_Allpairs_TCGA_PFS.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_THR50_1,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'k-TSPs and TCGA PFS (THR50_1)')
dev.off()

pdf("./figures/sep28/THR50_2_Allpairs_TCGA_PFS.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs_THR50_2,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = 'k-TSPs and TCGA PFS (THR50_2)')
dev.off()

#############
## fit coxph model:

# by pair
surv_func_TCGA_pfs_coxph <- function(x){
  f <- as.formula(paste("Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~", x))
  return(coxph(f, data = CoxData_tcga))
}

fit_list_TCGA_pfs_coxph_THR25 <- lapply(pairs_list_THR25, surv_func_TCGA_pfs_coxph)
names(fit_list_TCGA_pfs_coxph_THR25) <- pairs_list_THR25
summary_list_TCGA_pfs_coxph_THR25 <- lapply(fit_list_TCGA_pfs_coxph_THR25, summary)

fit_list_TCGA_pfs_coxph_THR50_1 <- lapply(pairs_list_THR50_1, surv_func_TCGA_pfs_coxph)
names(fit_list_TCGA_pfs_coxph_THR50_1) <- pairs_list_THR50_1
summary_list_TCGA_pfs_coxph_THR50_1 <- lapply(fit_list_TCGA_pfs_coxph_THR50_1, summary)

fit_list_TCGA_pfs_coxph_THR50_2 <- lapply(pairs_list_THR50_2, surv_func_TCGA_pfs_coxph)
names(fit_list_TCGA_pfs_coxph_THR50_2) <- pairs_list_THR50_2
summary_list_TCGA_pfs_coxph_THR50_2 <- lapply(fit_list_TCGA_pfs_coxph_THR50_2, summary)

# get the HR
HR_list_TCGA_pfs_coxph_THR25 <- lapply(summary_list_TCGA_pfs_coxph_THR25, function(x){
  HR <- x$conf.int[, 'exp(coef)']
  Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
  Pvalue_logrank_test <- x$sctest['pvalue']
  Pvalue_wald_test <- x$waldtest['pvalue']
  data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
             Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
})

HR_df_TCGA_pfs_coxph_THR25 <- as.data.frame(do.call(rbind, HR_list_TCGA_pfs_coxph_THR25))
HR_df_TCGA_pfs_coxph_THR25$variable <- rownames(HR_df_TCGA_pfs_coxph_THR25)
HR_df_TCGA_pfs_coxph_THR25 <- HR_df_TCGA_pfs_coxph_THR25[order(HR_df_TCGA_pfs_coxph_THR25$HR, decreasing = T), ]

#############
# get the HR
HR_list_TCGA_pfs_coxph_THR50_1 <- lapply(summary_list_TCGA_pfs_coxph_THR50_1, function(x){
  HR <- x$conf.int[, 'exp(coef)']
  Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
  Pvalue_logrank_test <- x$sctest['pvalue']
  Pvalue_wald_test <- x$waldtest['pvalue']
  data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
             Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
})

HR_df_TCGA_pfs_coxph_THR50_1 <- as.data.frame(do.call(rbind, HR_list_TCGA_pfs_coxph_THR50_1))
HR_df_TCGA_pfs_coxph_THR50_1$variable <- rownames(HR_df_TCGA_pfs_coxph_THR50_1)
HR_df_TCGA_pfs_coxph_THR50_1 <- HR_df_TCGA_pfs_coxph_THR50_1[order(HR_df_TCGA_pfs_coxph_THR50_1$HR, decreasing = T), ]

#############
# get the HR
HR_list_TCGA_pfs_coxph_THR50_2 <- lapply(summary_list_TCGA_pfs_coxph_THR50_2, function(x){
  HR <- x$conf.int[, 'exp(coef)']
  Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
  Pvalue_logrank_test <- x$sctest['pvalue']
  Pvalue_wald_test <- x$waldtest['pvalue']
  data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
             Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
})

HR_df_TCGA_pfs_coxph_THR50_2 <- as.data.frame(do.call(rbind, HR_list_TCGA_pfs_coxph_THR50_2))
HR_df_TCGA_pfs_coxph_THR50_2$variable <- rownames(HR_df_TCGA_pfs_coxph_THR50_2)
HR_df_TCGA_pfs_coxph_THR50_2 <- HR_df_TCGA_pfs_coxph_THR50_1[order(HR_df_TCGA_pfs_coxph_THR50_2$HR, decreasing = T), ]

# save the results
write.csv(HR_df_TCGA_pfs_coxph_THR25, 'objs/sep28/HR_df_TCGA_pfs_coxph_THR25.csv')
write.csv(HR_df_TCGA_pfs_coxph_THR50_1, 'objs/sep28/HR_df_TCGA_pfs_coxph_THR50_1.csv')
write.csv(HR_df_TCGA_pfs_coxph_THR50_2, 'objs/sep28/HR_df_TCGA_pfs_coxph_THR50_2.csv')


########
# by predictions
Fit_sig_TCGA_pfs_coxph_THR25 <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ prediction_tcga_THR25, data = CoxData_tcga)
summary(Fit_sig_TCGA_pfs_coxph_THR25)

Fit_sig_TCGA_pfs_coxph_THR50_1 <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ prediction_tcga_THR50_1, data = CoxData_tcga)
summary(Fit_sig_TCGA_pfs_coxph_THR50_1)

Fit_sig_TCGA_pfs_coxph_THR50_2 <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ prediction_tcga_THR50_2, data = CoxData_tcga)
summary(Fit_sig_TCGA_pfs_coxph_THR50_2)

#####################
png('./figures/sep28/THR25_HR_tcga_pfs.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_pfs_coxph_THR25)
dev.off()

png('./figures/sep28/THR50_1_HR_tcga_pfs.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_pfs_coxph_THR50_1)
dev.off()

png('./figures/sep28/THR50_2_HR_tcga_pfs.png', width = 2000, height = 2000, res = 300)
ggforest(Fit_sig_TCGA_pfs_coxph_THR50_2)
dev.off()



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
THR_signature <- read.delim("./data/THR_CellExpress_new.csv", sep = ',')

THR_signature_top350 <- readxl::read_xlsx("./data/THR_top_350.xlsx")

# sort by p-value in ascending order
THR_signature <- THR_signature[order(THR_signature$p.value, decreasing = FALSE), ]
THR_signature_top350 <- THR_signature_top350[order(THR_signature_top350$`p-value`, decreasing = FALSE), ]

# remove NAs
THR_signature <- THR_signature[!is.na(THR_signature$Gene), ]
THR_signature_top350 <- THR_signature_top350[!is.na(THR_signature_top350$Gene), ]

# get the top genes
#THR_signature <- THR_signature[c(1:100), ]

# make TSPs
#myTSPs <- t(combn(THR_signature$Gene, 2))
myTSPs <- t(combn(THR_signature_top350$Gene, 2))
colnames(myTSPs) <- c("gene1", "gene2")

################
# Load the  expression and pheno data
load('./objs/forKTSP.rda')

###################################

### Common genes
keepGns_datasets <- intersect(rownames(Expr_tcga), rownames(Expr_metabric))
keepGns_tsp <- intersect(as.vector(myTSPs), keepGns_datasets)

Expr_metabric <- Expr_metabric[keepGns_tsp, ]
Expr_tcga <- Expr_tcga[keepGns_tsp, ]

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns_tsp & myTSPs[,2] %in% keepGns_tsp , ]

###########################################################################
### TRAINING using restricted pairs
###########################################################################
### Set Feature number and max k
ktsp <- c(3:25) #7
featNo <- nrow(Expr_metabric)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333) # 9/7 # 25

ktsp_metabric <- SWAP.Train.KTSP(
  Expr_metabric, group_metabric, krange=10, disjoint = T, 
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=25, RestrictedPairs = myTSPs)

ktsp_metabric

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStats_metabric <- SWAP.KTSP.Statistics(inputMat = Expr_metabric, classifier = ktsp_metabric, CombineFunc = sum)
summary(ktspStats_metabric$statistics)

### Threshold
thr <- coords(roc(group_metabric, ktspStats_metabric$statistics, levels = c(0, 1), direction = "<" ), "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(group_metabric, ktspStats_metabric$statistics, levels = c(0, 1), direction = "<" ), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROC_metabric <- roc(group_metabric, ktspStats_metabric$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_metabric

### Get predictions based on best threshold from ROC curve
prediction_metabric <- SWAP.KTSP.Classify(Expr_metabric, ktsp_metabric, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(prediction_metabric, group_metabric, positive = '1', mode = "everything")

MCC_metabric <- mltools::mcc(pred = prediction_metabric, actuals = group_metabric)
MCC_metabric

########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStats_tcga <- SWAP.KTSP.Statistics(inputMat = Expr_tcga, classifier = ktsp_metabric, CombineFunc = sum)
summary(ktspStats_tcga$statistics)

## Plot curve
ROC_tcga <- roc(group_tcga, ktspStats_tcga$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")
ROC_tcga

### Get predictions based on best threshold from ROC curve
prediction_tcga <- SWAP.KTSP.Classify(Expr_tcga, ktsp_metabric, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusion_tcga <- confusionMatrix(prediction_tcga, group_tcga, positive = "1", mode = "everything")
confusion_tcga

MCC_tcga <- mltools::mcc(pred = prediction_tcga, actuals = group_tcga)
MCC_tcga


############################################################
############################################################
############################################################
# test the ktsp pairs with survival
ClassifierGenes <- as.vector(ktsp_metabric$TSPs)


##########################
## Keep only the relevant information (Metastasis Event and Time)
Phenotype_metabric <- cbind(Pheno_metabric[, c("Overall.Survival.Status", "Overall.Survival..Months.")], 
                            ktspStats_metabric$comparisons, prediction_metabric)

Phenotype_tcga <- cbind(Pheno_tcga[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Progression.Free.Status", "Progress.Free.Survival..Months.")], 
                        ktspStats_tcga$comparisons, prediction_tcga)

#Expr_metabric <- Expr_metabric[ClassifierGenes, ]
#Expr_tcga <- Expr_tcga[ClassifierGenes, ]


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
pairs_list <- list()

for(i in seq(1,nrow(ktsp_metabric$TSPs))){
  pairs_list[i] <- paste0(ktsp_metabric$TSPs[i,1], '.', ktsp_metabric$TSPs[i,2])
}

names(pairs_list) <- pairs_list

# fit surv curves
#fit_list <- list()

surv_func_metabric_os <- function(x){
  f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
  return(surv_fit(f, data = CoxData_metabric))
}

fit_list_metabric_os <- lapply(pairs_list, surv_func_metabric_os)

# calculate the pvalue
Pval_list_metabric_os <- surv_pvalue(fit_list_metabric_os)

Pval_df_metabric_os <- do.call(rbind.data.frame, Pval_list_metabric_os)

#Pval_df_fil <- Pval_df[Pval_df$pval < 0.05, ] 

############
## Plot survival curves

plot_list_metabric_os <- ggsurvplot_list(fit_list_metabric_os, CoxData_metabric, legend.title = names(fit_list_metabric_os), pval = TRUE)


Splot_metabric_os <- arrange_ggsurvplots(plot_list_metabric_os, title = "Survival plots using the 10 pairs individually", ncol = 2, nrow = 5)
ggsave("./figures/10TSPs_Metabric_os.pdf", Splot_metabric_os, width = 40, height = 40, units = "cm")


## metabric all pairs
Fit_sig_metabric <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_metabric, data = CoxData_metabric)
surv_pvalue(Fit_sig_metabric)

pdf("./figures/THRtop350_Allpairs_metabric.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = '10-TSPs and METABRIC OS')
dev.off()

#############
## fit coxph model:

# by pair
surv_func_metabric_os_coxph <- function(x){
  f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
  return(coxph(f, data = CoxData_metabric))
}

fit_list_metabric_os_coxph <- lapply(pairs_list, surv_func_metabric_os_coxph)
names(fit_list_metabric_os_coxph) <- pairs_list

summary_list_metabric_os_coxph <- lapply(fit_list_metabric_os_coxph, summary)

# get the HR
HR_list_metabric_os_coxph <- lapply(summary_list_metabric_os_coxph, function(x){
  HR <- x$conf.int[, 'exp(coef)']
  Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
  Pvalue_logrank_test <- x$sctest['pvalue']
  Pvalue_wald_test <- x$waldtest['pvalue']
  data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
             Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
})


HR_df_metabric_os_coxph <- as.data.frame(do.call(rbind, HR_list_metabric_os_coxph))
HR_df_metabric_os_coxph$variable <- rownames(HR_df_metabric_os_coxph)
HR_df_metabric_os_coxph <- HR_df_metabric_os_coxph[order(HR_df_metabric_os_coxph$HR, decreasing = T), ]

# save the results
write.csv(HR_df_metabric_os_coxph, 'objs/HR_df_metabric_os_coxph.csv')


########
# by predictions
Fit_sig_metabric_coxph <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_metabric, data = CoxData_metabric)
summary(Fit_sig_metabric_coxph)

ggforest(Fit_sig_metabric_coxph)

########################################################################  
########################################################################  
## Fit survival curves: TCGA: 

# fit surv curves: os

surv_func_TCGA_os <- function(x){
  f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
  return(surv_fit(f, data = CoxData_tcga))
}

fit_list_tcga_os <- lapply(pairs_list, surv_func_TCGA_os)

# calculate the pvalue
Pval_list_tcga_os <- surv_pvalue(fit_list_tcga_os)

Pval_df_tcga_os <- do.call(rbind.data.frame, Pval_list_tcga_os)

#Pval_df_fil <- Pval_df[Pval_df$pval < 0.05, ] 

############
## Plot survival curves

plot_list_tcga_os <- ggsurvplot_list(fit_list_tcga_os, CoxData_tcga, legend.title = names(fit_list_tcga_os), pval = TRUE)


Splot_tcga_os <- arrange_ggsurvplots(plot_list_tcga_os, title = "Overall survival in the TCGA using the 10 pairs individually", ncol = 2, nrow = 5)
ggsave("./figures/10TSPs_tcga_os.pdf", Splot_tcga_os, width = 40, height = 40, units = "cm")


## TCGA all pairs
Fit_sig_TCGA_os <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_tcga, data = CoxData_tcga)

pdf("./figures/THRtop350_Allpairs_TCGA_os.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_os,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = '10-TSPs and TCGA OS')
dev.off()

#############
## fit coxph model:

# by pair
surv_func_TCGA_os_coxph <- function(x){
  f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
  return(coxph(f, data = CoxData_tcga))
}

fit_list_TCGA_os_coxph <- lapply(pairs_list, surv_func_TCGA_os_coxph)
names(fit_list_TCGA_os_coxph) <- pairs_list

summary_list_TCGA_os_coxph <- lapply(fit_list_TCGA_os_coxph, summary)

# get the HR
HR_list_TCGA_os_coxph <- lapply(summary_list_TCGA_os_coxph, function(x){
  HR <- x$conf.int[, 'exp(coef)']
  Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
  Pvalue_logrank_test <- x$sctest['pvalue']
  Pvalue_wald_test <- x$waldtest['pvalue']
  data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
             Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
})


HR_df_TCGA_os_coxph <- as.data.frame(do.call(rbind, HR_list_TCGA_os_coxph))
HR_df_TCGA_os_coxph$variable <- rownames(HR_df_TCGA_os_coxph)
HR_df_TCGA_os_coxph <- HR_df_TCGA_os_coxph[order(HR_df_TCGA_os_coxph$HR, decreasing = T), ]

# save the results
write.csv(HR_df_TCGA_os_coxph, 'objs/HR_df_TCGA_os_coxph.csv')


########
# by predictions
Fit_sig_TCGA_coxph <- coxph(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_tcga, data = CoxData_tcga)
summary(Fit_sig_TCGA_coxph)

ggforest(Fit_sig_TCGA_coxph)

########################################################################  
########################################################################  
## Fit survival curves: TCGA: 

# fit surv curves: PFS

surv_func_TCGA_pfs <- function(x){
  f <- as.formula(paste("Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~", x))
  return(surv_fit(f, data = CoxData_tcga))
}

fit_list_tcga_pfs <- lapply(pairs_list, surv_func_TCGA_pfs)

# calculate the pvalue
Pval_list_tcga_pfs <- surv_pvalue(fit_list_tcga_pfs)

Pval_df_tcga_pfs <- do.call(rbind.data.frame, Pval_list_tcga_pfs)

#Pval_df_fil <- Pval_df[Pval_df$pval < 0.05, ] 

############
## Plot survival curves

plot_list_tcga_pfs <- ggsurvplot_list(fit_list_tcga_pfs, CoxData_tcga, legend.title = names(fit_list_tcga_pfs), pval = TRUE)


Splot_tcga_pfs <- arrange_ggsurvplots(plot_list_tcga_pfs, title = "Progression free survival in the TCGA using the 10 pairs individually", ncol = 2, nrow = 5)
ggsave("./figures/10TSPs_tcga_pfs.pdf", Splot_tcga_pfs, width = 40, height = 40, units = "cm")


## TCGA all pairs
Fit_sig_TCGA_pfs <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ prediction_tcga, data = CoxData_tcga)

pdf("./figures/THRtop350_Allpairs_TCGA_PFS.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_TCGA_pfs,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = '10-TSPs and TCGA PFS')
dev.off()

#############
## fit coxph model:

# by pair
surv_func_TCGA_pfs_coxph <- function(x){
  f <- as.formula(paste("Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~", x))
  return(coxph(f, data = CoxData_tcga))
}

fit_list_TCGA_pfs_coxph <- lapply(pairs_list, surv_func_TCGA_pfs_coxph)
names(fit_list_TCGA_pfs_coxph) <- pairs_list

summary_list_TCGA_pfs_coxph <- lapply(fit_list_TCGA_pfs_coxph, summary)

# get the HR
HR_list_TCGA_pfs_coxph <- lapply(summary_list_TCGA_pfs_coxph, function(x){
  HR <- x$conf.int[, 'exp(coef)']
  Pvalue_Likelihood_ratio_test <- x$logtest['pvalue']
  Pvalue_logrank_test <- x$sctest['pvalue']
  Pvalue_wald_test <- x$waldtest['pvalue']
  data.frame(HR = HR, Pvalue_Likelihood_ratio_test = Pvalue_Likelihood_ratio_test, 
             Pvalue_logrank_test = Pvalue_logrank_test, Pvalue_wald_test = Pvalue_wald_test)
})


HR_df_TCGA_pfs_coxph <- as.data.frame(do.call(rbind, HR_list_TCGA_pfs_coxph))
HR_df_TCGA_pfs_coxph$variable <- rownames(HR_df_TCGA_pfs_coxph)
HR_df_TCGA_pfs_coxph <- HR_df_TCGA_pfs_coxph[order(HR_df_TCGA_pfs_coxph$HR, decreasing = T), ]

# save the results
write.csv(HR_df_TCGA_pfs_coxph, 'objs/HR_df_TCGA_pfs_coxph.csv')


########
# by predictions
Fit_sig_TCGA_pfs_coxph <- coxph(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ prediction_tcga, data = CoxData_tcga)
summary(Fit_sig_TCGA_pfs_coxph)

ggforest(Fit_sig_TCGA_pfs_coxph)


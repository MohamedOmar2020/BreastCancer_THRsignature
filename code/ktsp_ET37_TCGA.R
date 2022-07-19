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
library(survminer)
library(survival)

#################
load("./objs/ET_metabric.rda")

## get the ET125 genes
All <- read_xlsx('misc/ET-9 Selection Steps.xlsx')
# ET60 <- All$`ET-60`
# 
# ET60[ET60=='GPR56'] <- 'ADGRG1'
# ET60[ET60=='FAM116B'] <- 'DENND6B'
# #ET125[ET125=='FAM46B'] <- 'TENT5B'
# ET60[ET60=='HSA011916'] <- 'CTDNEP1'

# ET60 <- ET60[!is.na(ET60)]

myTSPs <- t(combn(All$`ET-125`, 2))
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
ktsp <- c(3:10) #7
featNo <- nrow(Expr_metabric)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktsp_tcga <- SWAP.Train.KTSP(
  Expr_tcga, group_tcga, krange=10, disjoint = T, 
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)

ktsp_tcga

############################################################
# Get the genes in common with ET9
ClassifierGenes <- as.vector(ktsp_tcga$TSPs)

# intersection with ET9
TSPs_ET9_int <- intersect(All$`ET-9`, ClassifierGenes)
TSPs_ET9_int

# ET4 intersection with ET9
TSPs_ET4_int <- intersect(ClassifierGenes, ET4)
TSPs_ET4_int

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStats_tcga <- SWAP.KTSP.Statistics(inputMat = Expr_tcga, classifier = ktsp_tcga, CombineFunc = sum)
summary(ktspStats_tcga$statistics)

### Threshold
thr <- coords(roc(group_tcga, ktspStats_tcga$statistics, levels = c(0, 1), direction = "<" ), "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(group_tcga, ktspStats_tcga$statistics, levels = c(0, 1), direction = "<" ), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROC_tcga <- roc(group_tcga, ktspStats_tcga$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE)
ROC_tcga

### Get predictions based on best threshold from ROC curve
prediction_tcga <- SWAP.KTSP.Classify(Expr_tcga, ktsp_tcga, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(prediction_tcga, group_tcga, positive = '1', mode = "everything")

MCC_tcga <- mltools::mcc(pred = prediction_tcga, actuals = group_tcga)
MCC_tcga

########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStats_metabric <- SWAP.KTSP.Statistics(inputMat = Expr_metabric, classifier = ktsp_tcga, CombineFunc = sum)
summary(ktspStats_metabric$statistics)

## Plot curve
ROC_metabric <- roc(group_metabric, ktspStats_metabric$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c(0, 1), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")
ROC_metabric

### Get predictions based on best threshold from ROC curve
prediction_metabric <- SWAP.KTSP.Classify(Expr_metabric, ktsp_tcga, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusion_metabric <- confusionMatrix(prediction_metabric, group_metabric, positive = "1", mode = "everything")
confusion_metabric

MCC_metabric <- mltools::mcc(pred = prediction_metabric, actuals = group_metabric)
MCC_metabric


############################################################
## Keep only the relevant information (Metastasis Event and Time)
Phenotype_metabric <- cbind(Pheno_metabric[, c("Overall.Survival.Status", "Overall.Survival..Months.")], 
                            ktspStats_metabric$comparisons, prediction_metabric)

Phenotype_tcga <- cbind(Pheno_tcga[, c("Progression.Free.Status", "Progress.Free.Survival..Months.")], 
                        ktspStats_tcga$comparisons, prediction_tcga)

# create a merged pdata and Z-scores object
CoxData_metabric <- data.frame(Phenotype_metabric)
CoxData_tcga <- data.frame(Phenotype_tcga)

########################################################################  
## Fit survival curves

pairs_list <- as.list(colnames(CoxData_metabric)[3:(ncol(CoxData_metabric)-1)])
names(pairs_list) <- colnames(CoxData_metabric)[3:(ncol(CoxData_metabric)-1)]

# Metabric: 
surv_func_metabric <- function(x){
  f <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~", x))
  return(surv_fit(f, data = CoxData_metabric))
}
metabric_fit_list <- lapply(pairs_list, surv_func_metabric)
fit_sig_metabric <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_metabric, data = CoxData_metabric)

# TCGA: 
surv_func_tcga <- function(x){
  f <- as.formula(paste("Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~", x))
  return(surv_fit(f, data = CoxData_tcga))
}

tcga_fit_list <- lapply(pairs_list, surv_func_tcga)
fit_sig_tcga <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ prediction_tcga, data = CoxData_tcga)


##################################################################
## Plot survival curves

## Metabric:
# all pairs combined
pdf("./figs/10TSPs_Allpairs_metabric.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(fit_sig_metabric,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = '10-TSPs combined (Metabric OS)')
dev.off()

# individual pairs
metabric_plot_list <- ggsurvplot_list(metabric_fit_list,
                                      data = CoxData_metabric,
                                      title =gsub("\\.", "\\>", names(metabric_fit_list)), 
                                      pval = TRUE)


Splot_metabric <- arrange_ggsurvplots(metabric_plot_list, title = "Metabric OS using the 10 TSPs pairs individually", ncol = 5, nrow = 2)
ggsave("./figs/10TSPs_indvPairs_Metabric.pdf", Splot_metabric, width = 45, height = 25, units = "cm")

##########
## TCGA:
# all pairs combined
pdf("./figs/10TSPs_Allpairs_tcga.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(fit_sig_tcga,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = '10-TSPs combined (TCGA PFS)')
dev.off()

# individual pairs
tcga_plot_list <- ggsurvplot_list(tcga_fit_list,
                                      data = CoxData_tcga,
                                      title =gsub("\\.", "\\>", names(tcga_fit_list)), 
                                      pval = TRUE)


Splot_tcga <- arrange_ggsurvplots(tcga_plot_list, title = "TCGA PFS using the 10 TSPs pairs individually", ncol = 5, nrow = 2)
ggsave("./figs/10TSPs_indvPairs_TCGA.pdf", Splot_tcga, width = 45, height = 25, units = "cm")




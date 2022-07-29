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

#################
THR_signature <- read.delim("./data/THR_CellExpress.csv", sep = ',')

# sort by p-value in ascending order
THR_signature <- THR_signature[order(THR_signature$p.value, decreasing = FALSE), ]

# remove NAs
THR_signature <- THR_signature[!is.na(THR_signature$Gene), ]

# get the top 50
#THR_signature <- THR_signature[c(1:50), ]

# make TSPs
myTSPs <- t(combn(THR_signature$Gene, 2))
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
set.seed(333)

ktsp_tcga <- SWAP.Train.KTSP(
  Expr_tcga, group_tcga, krange=ktsp, disjoint = T, 
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=50, RestrictedPairs = myTSPs)

ktsp_tcga

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
############################################################
############################################################
# test the ktsp pairs with survival
ClassifierGenes <- as.vector(ktsp_tcga$TSPs)

# intersection with ET9
All <- read_xlsx('misc/ET-9 Selection Steps.xlsx')

TSPs_ET9_int <- intersect(All$`ET-9`, ClassifierGenes)
TSPs_ET9_int

# ET4 intersection with ET9
ET9_ET4_int <- intersect(All$`ET-9`, ET4)
ET9_ET4_int
##########################
## Keep only the relevant information (Metastasis Event and Time)
Phenotype_metabric <- cbind(Pheno_metabric[, c("Overall.Survival.Status", "Overall.Survival..Months.")], 
                            ktspStats_metabric$comparisons, prediction_metabric)

Phenotype_tcga <- cbind(Pheno_tcga[, c("Progression.Free.Status", "Progress.Free.Survival..Months.")], 
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
Fit_IGFBP5_ID3_metabric <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ IGFBP5.ID3 , data = CoxData_metabric)
Fit_GDPD5_S100A6_metabric <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ GDPD5.S100A6, data = CoxData_metabric)
Fit_CCDC69_IER5_metabric <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ IER5L.CCDC69, data = CoxData_metabric)
Fit_CUX1_CX3CL1_metabric <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ CUX1.CX3CL1, data = CoxData_metabric)

Fit_sig_metabric <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ prediction_metabric, data = CoxData_metabric)


##################################################################
##################################################################
## Plot survival curves

plot_IGFBP5_ID3_metabric <- ggsurvplot(Fit_IGFBP5_ID3_metabric,
                                       risk.table = FALSE,
                                       pval = TRUE,
                                       ggtheme = theme_minimal(),
                                       risk.table.y.text.col = FALSE,
                                       risk.table.y.text = FALSE, title = "IGFBP5 > ID3")

plot_GDPD5_S100A6_metabric <- ggsurvplot(Fit_GDPD5_S100A6_metabric,
                                         risk.table = FALSE,
                                         pval = TRUE,
                                         ggtheme = theme_minimal(),
                                         risk.table.y.text.col = FALSE,
                                         risk.table.y.text = FALSE, title = "GDPD5 > S100A6")

plot_CCDC69_IER5_metabric <- ggsurvplot(Fit_CCDC69_IER5_metabric,
                                        risk.table = FALSE,
                                        pval = TRUE,
                                        ggtheme = theme_minimal(),
                                        risk.table.y.text.col = FALSE,
                                        risk.table.y.text = FALSE, title = "IER5L > CCDC69")

plot_CUX1_CX3CL1_metabric <- ggsurvplot(Fit_CUX1_CX3CL1_metabric,
                                        risk.table = FALSE,
                                        pval = TRUE,
                                        ggtheme = theme_minimal(),
                                        risk.table.y.text.col = FALSE,
                                        risk.table.y.text = FALSE, title = "CUX1 > CX3CL1")


pdf("./figs/ET4_Allpairs_metabric.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_metabric,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = '4-TSPs and METABRIC OS')
dev.off()

PlotList <- list(plot_IGFBP5_ID3_metabric, plot_GDPD5_S100A6_metabric, plot_CCDC69_IER5_metabric, plot_CUX1_CX3CL1_metabric)
names(PlotList) <- rownames(ktsp_metabric$TSPs)

Splot <- arrange_ggsurvplots(PlotList, title = "Survival plots using the ET4 pairs individually", ncol = 2, nrow = 2)
ggsave("./figs/ET4_Pairs_Metabric.pdf", Splot, width = 40, height = 30, units = "cm")


########################################################################  
## Fit survival curves: TCGA: 
Fit_IGFBP5_ID3_tcga <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ IGFBP5.ID3 , data = CoxData_tcga)
Fit_GDPD5_S100A6_tcga <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ GDPD5.S100A6, data = CoxData_tcga)
Fit_CCDC69_IER5_tcga <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ IER5L.CCDC69, data = CoxData_tcga)
Fit_CUX1_CX3CL1_tcga <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ CUX1.CX3CL1, data = CoxData_tcga)

Fit_sig_tcga <- survfit(Surv(Progress.Free.Survival..Months., Progression.Free.Status) ~ prediction_tcga, data = CoxData_tcga)


#########
## Plot survival curves

plot_IGFBP5_ID3_tcga <- ggsurvplot(Fit_IGFBP5_ID3_tcga,
                                   risk.table = FALSE,
                                   pval = TRUE,
                                   ggtheme = theme_minimal(),
                                   risk.table.y.text.col = FALSE,
                                   risk.table.y.text = FALSE, title = "IGFBP5 > ID3")

plot_GDPD5_S100A6_tcga <- ggsurvplot(Fit_GDPD5_S100A6_tcga,
                                     risk.table = FALSE,
                                     pval = TRUE,
                                     ggtheme = theme_minimal(),
                                     risk.table.y.text.col = FALSE,
                                     risk.table.y.text = FALSE, title = "GDPD5 > S100A6")

plot_CCDC69_IER5_tcga <- ggsurvplot(Fit_CCDC69_IER5_tcga,
                                    risk.table = FALSE,
                                    pval = TRUE,
                                    ggtheme = theme_minimal(),
                                    risk.table.y.text.col = FALSE,
                                    risk.table.y.text = FALSE, title = "IER5L > CCDC69")

plot_CUX1_CX3CL1_tcga <- ggsurvplot(Fit_CUX1_CX3CL1_tcga,
                                    risk.table = FALSE,
                                    pval = TRUE,
                                    ggtheme = theme_minimal(),
                                    risk.table.y.text.col = FALSE,
                                    risk.table.y.text = FALSE, title = "CUX1 > CX3CL1")


pdf("./figs/ET4_Allpairs_tcga.pdf", width = 8, height = 8, onefile = F)
ggsurvplot(Fit_sig_tcga,
           risk.table = FALSE,
           pval = TRUE,
           ggtheme = theme_minimal(),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, title = '4-TSPs and TCGA PFS')
dev.off()

PlotList <- list(plot_IGFBP5_ID3_tcga, plot_GDPD5_S100A6_tcga, plot_CCDC69_IER5_tcga, plot_CUX1_CX3CL1_tcga)
names(PlotList) <- rownames(ktsp_metabric$TSPs)

Splot <- arrange_ggsurvplots(PlotList, title = "Survival plots using the ET4 pairs individually (TCGA PFS)", ncol = 2, nrow = 2)
ggsave("./figs/ET4_Pairs_tcga.pdf", Splot, width = 40, height = 30, units = "cm")

###############################################################
###############################################################
################################################################
## Plot Cox Proportional Hazard Model
CoxData$T_Stage <- as.factor(Pheno$Tumor.Stage)
levels(CoxData$T_Stage)

CoxData$Radiotherapy <- as.factor(Pheno$Radio.Therapy)
levels(CoxData$Radiotherapy)

CoxData$ChemoTherapy <- as.factor(Pheno$Chemotherapy)
levels(CoxData$ChemoTherapy)

CoxData$Age <- Pheno$Age.at.Diagnosis
#ifelse(Pheno$Age.at.Diagnosis < 50, "<50", ">=50")


CoxData$Surgery <- as.factor(Pheno$Type.of.Breast.Surgery)
levels(CoxData$Surgery) 

CoxData$Menopausal_status <- as.factor(Pheno$Inferred.Menopausal.State)
levels(CoxData$Menopausal_status) 

CoxData$Cellularity <- as.factor(Pheno$Cellularity)
levels(CoxData$Cellularity) 

CoxData$Laterality <- as.factor(Pheno$Primary.Tumor.Laterality)
levels(CoxData$Laterality) 

# CoxData$ResidualDisease_Score <- ktspStatsTestRes$statistics
# CoxData$ResidualDisease_Score <- ifelse(CoxData$ResidualDisease_Score >= 2.5, "high", "low")
# CoxData$ResidualDisease_Score <- factor(CoxData$ResidualDisease_Score, levels = c("low", "high"))
# table(CoxData$ResidualDisease_Score)

# fit multivariable Cox model 
CoxModel_indPairs <- coxph(Surv(Time, Event) ~ GARS.PDCD10+CCND1.FNDC3A+CMA1.PTHLH+PCDH7.TOX+T_Stage+Radiotherapy+Age+Surgery+Menopausal_status, data = CoxData)
CoxModel_indPairs

CoxModel_AllPairs <- coxph(Surv(Time, Event) ~ ktspPredTest+T_Stage+Radiotherapy+Age+Surgery+Menopausal_status, data = CoxData)
CoxModel_AllPairs

pdf("./Figs/survival/CoxModel_METABRIC_IndPairs.pdf", width = 10, height = 10, onefile = F)
ggforest(CoxModel_indPairs, data = CoxData)
dev.off()

pdf("./Figs/survival/CoxModel_METABRIC_Allpairs.pdf", width = 10, height = 10, onefile = F)
ggforest(CoxModel_AllPairs, data = CoxData)
dev.off()

# ## Univariate Cox model
covariates <- c("ResidualDisease_Score", "T_Stage", "Radiotherapy", "Age", "Surgery", "ChemoTherapy")
# 
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(Time, Event)~', x)))
# 
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = SurvData)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (",
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })

res <- t(as.data.frame(univ_results, check.names = F))
as.data.frame(res)


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
rownames(Expr_metabric)[grep('^ZNF652', rownames(Expr_metabric))]

# filter the THR signatures to include only the genes present in the expr matrices
PAM50_fil <- PAM50[PAM50 %in% rownames(Expr_metabric)]

setdiff(PAM50, PAM50_fil)

#############################################################################################
### combine in 1 dataset: Training
Data_metabric <- as.data.frame(cbind(t(Expr_metabric), group_metabric))
Data_metabric$group_metabric <- as.factor(Data_metabric$group_metabric)
levels(Data_metabric$group_metabric) <- c('0', '1')
colnames(Data_metabric)[colnames(Data_metabric) %in% c('group_metabric')] <- c('os')

###########################################################################
### TRAINING using logistic regression
###########################################################################
#############################################################################################################
##############################################################################################################
# the model

PAM50_model <- glm(as.formula((paste("os ~", paste(PAM50_fil, collapse = "+")))), data = Data_metabric, family = "binomial")
summary(PAM50_model)

save(PAM50_fil, file = "./objs/PAM50_model_logreg.rda")



###########################################################################
############################################################################
### predict in the metabric
# Make predictions

Train_prob_PAM50 <- PAM50_model %>% predict(Data_metabric , type = "response")

### Threshold
thr_PAM50 <- coords(roc(group_metabric, Train_prob_PAM50, levels = c("0", "1"), direction = "<"), "best")["threshold"]
thr_PAM50

### ROC Curve
ROCTrain_PAM50 <- roc(group_metabric, Train_prob_PAM50, plot = F, print.thres=thr_PAM50$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("0", "1"), direction = "<", col="blue", lwd=2, grid=TRUE)
ROCTrain_PAM50

### Get predictions based on best threshold from ROC curve
Train_predClasses_PAM50 <- ifelse(Train_prob_PAM50 >= thr_PAM50$threshold, "1", "0")
table(Train_predClasses_PAM50)
Train_predClasses_PAM50 <- factor(Train_predClasses_PAM50, levels = c('0', '1'))


### Resubstitution performance in the TRAINING set
ConfusionTrain_PAM50 <- confusionMatrix(Train_predClasses_PAM50, group_metabric, positive = "1", mode = "everything")
ConfusionTrain_PAM50

## MCC
MCC_Train_PAM50 <- mltools::mcc(pred = Train_predClasses_PAM50, actuals = group_metabric)
MCC_Train_PAM50

##########################
## Keep only the relevant information (Metastasis Event and Time)
Phenotype_metabric <- cbind(Pheno_metabric[, c("Overall.Survival.Status", "Overall.Survival..Months.", "Relapse.Free.Status", "Relapse.Free.Status..Months.", "Pam50...Claudin.low.subtype", "ER.status.measured.by.IHC", "X3.Gene.classifier.subtype")], 
                            Train_prob_PAM50, Train_predClasses_PAM50)



# create a merged pdata and Z-scores object
CoxData_metabric <- data.frame(Phenotype_metabric)

# divide the probabilities into quartiles
CoxData_metabric <- CoxData_metabric %>%
  mutate(metabric_prob_PAM50_quartiles = ntile(Train_prob_PAM50, 4), 
         metabric_prob_PAM50_quintiles = ntile(Train_prob_PAM50, 5))


# OS
## metabric all genes
Fit_sig_metabric_os_PAM50 <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ Train_predClasses_PAM50, data = CoxData_metabric)


## by quartiles
Fit_sig_metabric_os_PAM50_quartiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_PAM50_quartiles, data = CoxData_metabric)

## by quintiles
Fit_sig_metabric_os_PAM50_quintiles <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ metabric_prob_PAM50_quintiles, data = CoxData_metabric)


# RFS
## metabric all genes
Fit_sig_metabric_RFS_PAM50 <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ Train_predClasses_PAM50, data = CoxData_metabric)

## by quartiles
Fit_sig_metabric_RFS_PAM50_quartiles <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_PAM50_quartiles, data = CoxData_metabric)

## by quintiles
Fit_sig_metabric_RFS_PAM50_quintiles <- survfit(Surv(Relapse.Free.Status..Months., Relapse.Free.Status) ~ metabric_prob_PAM50_quintiles, data = CoxData_metabric)

############################################################################
############################################################################
# plot OS

tiff("./figures/PAM50_logreg/PAM50_metabric_os_allpairs.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_PAM50,
           risk.table = FALSE,
           pval = FALSE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           palette = 'jco',
           risk.table.y.text = FALSE, 
           #title = 'THR 50 and METABRIC OS'
)
dev.off()

########
# by quartiles
tiff("./figures/PAM50_logreg/PAM50_metabric_os_quartiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_PAM50_quartiles,
           risk.table = FALSE,
           pval = FALSE,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR 50_1 (logistic regression) and METABRIC OS: quartiles'
)
dev.off()

#############
# by quintiles
tiff("./figures/PAM50_logreg/PAM50_metabric_os_quintiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_os_PAM50_quintiles,
           risk.table = FALSE,
           pval = FALSE,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           #title = 'THR-70 and METABRIC OS: quintiles'
)
dev.off()

######################################
# plot RFS

tiff("./figures/PAM50_logreg/PAM50_metabric_RFS_allpairs.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_PAM50,
           risk.table = FALSE,
           pval = FALSE,
           legend.labs = c('prediction: 0', 'prediction: 1'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           #title = 'THR 50_1 (logistic regression) and METABRIC RFS'
)
dev.off()

########
# by quartiles
tiff("./figures/PAM50_logreg/PAM50_metabric_RFS_quartiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_PAM50_quartiles,
           risk.table = FALSE,
           pval = FALSE,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE, 
           palette = 'jco',
           #title = 'THR 50_1 (logistic regression) and METABRIC RFS: quartiles'
)
dev.off()

########
# by quintiles

tiff("./figures/PAM50_logreg/PAM50_metabric_rfs_quintiles.tiff", width = 3000, height = 3000, res = 300)
ggsurvplot(Fit_sig_metabric_RFS_PAM50_quintiles,
           risk.table = FALSE,
           pval = FALSE,
           legend.labs = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5'),
           ggtheme = theme_survminer(base_size = 30, font.x = c(30, 'bold.italic', 'black'), font.y = c(30, 'bold.italic', 'black'), font.tickslab = c(30, 'plain', 'black'), font.legend = c(30, 'bold', 'black')),
           palette = 'jco',
           risk.table.y.text.col = FALSE,
           risk.table.y.text = FALSE 
           #title = 'THR-70 and METABRIC RFS: quintiles'
)
dev.off()


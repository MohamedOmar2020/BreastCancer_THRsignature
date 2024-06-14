
rm(list = ls())

gc()

library(Seurat)
library(AUCell)
library(pochi)
library(scMiko)
library(UCell)


# Created a vector of labels For cells in the endothelial subset, labelled them as "EC", for others label them as "Other"
labels <- ifelse(cellnames_all %in% cellnames_ec, "EC", "Other")

# Added the labels as a new column to the metadata
everything_ec <- AddMetaData(everything, metadata = labels, col.name = "Labels")

# Checked if the new column was added
head(everything_ec@meta.data, n = 2)


Idents(everything_ec) <- "Labels"

##################
# load the data
sc_epithelium <- readRDS('data/scRNAseq/HBCA_epithelium.rds')
dim(sc_epithelium)

##################
# change the gene ids from ens to symbols
counts <- GetAssayData(sc_epithelium, assay = "RNA", layer = "counts") 
rownames(counts) <- sc_epithelium@assays$RNA@meta.features$feature_name 
#save(counts, file = './data/scRNAseq/b8b5be07-061b-4390-af0a-f9ced877a068_counts.rda')

sc_epithelium <- CreateSeuratObject(counts = counts, meta.data = sc_epithelium@meta.data, min.cells = 3, min.features = 200)
dim(sc_epithelium)

rm(counts)
gc()

####################################
# re-process and re-compute umap
####################################

sc_epithelium <- NormalizeData(sc_epithelium, normalization.method = "LogNormalize", scale.factor = 10000)
sc_epithelium <- FindVariableFeatures(sc_epithelium, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sc_epithelium), 10)

# plot variable features with and without labels
# plot1 <- VariableFeaturePlot(sc_epithelium)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2

sc_epithelium <- ScaleData(sc_epithelium)
gc()

sc_epithelium <- RunPCA(sc_epithelium, features = VariableFeatures(object = sc_epithelium))
sc_epithelium <- FindNeighbors(sc_epithelium, dims = 1:10)
sc_epithelium <- RunUMAP(sc_epithelium, dims = 1:20)
DimPlot(object = sc_epithelium, group.by = "cell_type", label = TRUE)
gc()

save(sc_epithelium, file = './data/scRNAseq/HBCA_epithelium_processed.rda')

# load
load('./data/scRNAseq/HBCA_epithelium_processed.rda')

gc()

#############################
# Define the signatures
THR_signature <- readxl::read_xlsx("./data/THR Signatures_sep23.xlsx")

# get the THR50 signature
THR50 <- THR_signature$`THR-50.1`[!is.na(THR_signature$`THR-50.1`)]
THR50 <- gsub('-', '', THR50)
THR50 <- THR50[ - which(THR50 %in% 'Gene')]
THR50

# get the THR70 signature
THR_signature <- readxl::read_xlsx("./data/THR_Signatures_Jan25_2023.xlsx")
# get the THR50 signature
THR70 <- THR_signature$`THR-70`[!is.na(THR_signature$`THR-70`)]
THR70 <- gsub('-', '', THR70)



# oncotype
oncotype <-  c("MKI67", "AURKA", "BIRC5", "CCNB1", "MYBL2", "MMP11", "CTSL2", 
               "GRB7", "ERBB2", 
               "ESR1", "PGR", "BCL2", "SCUBE2", "GSTM1", "BAG1", "CD68",
               "ACTB", "GAPDH", "GUSB", "RPLP0", "TFRC")


# PAM50
PAM50_signature <- readxl::read_xlsx("./data/THR_Signatures_Jan25_2023.xlsx")
PAM50 <- PAM50_signature$PAM50[!is.na(PAM50_signature$PAM50)]
PAM50 <- gsub('-', '', PAM50)

# MammaPrint
mammaprint <-  c(
  "BBC3", "EGLN1", "ESM1", "IGFBP5", "SCUBE2",
  "TGFB3", "WISP1", "FLT1", "HRASLS", "STK32B", "RASSF7", "DCK",
  "MELK", "EXT1", "GNAZ", "EBF4", "MTDH", "PITRM1", "QSOX2",
  "CCNE2", "ECT2", "CENPA", "LIN9", "NDC80", "MCM6", "NUSAP1",
  "ORC6", "TSPYL5", "RUNDC1", "PRC1", "RFC4", "RECQL5", "CDCA7",
  "DTL", "COL4A2", "GPR180", "MMP9", "GPR126", "RTN4RL1", "DIAPH3",
  "CDC42BPA", "PALM2", "ALDH4A1", "LPCAT1", "OXCT1", "ECI2", "GMPS",
  "GSTM3", "SLC2A3", "FGF18", "COL4A2", "GPR180", "EGLN1",
  "MMP9", "LOC100288906", "C9orf30", "ZNF385B", "C16orf61", "SERF1A",
  "TMEM74B", "LOC730018", "LOC100131053", "AA555029_RC", "DHX58",
  "NMU", "UCHL5", "KDM7A", "AP2B1", "MS4A7", "RAB6B"
)

#########################
# make a list of all signatures
signatures <- list(
  THR50 = THR50,  
  THR70 = THR70,
  OncotypeDx = oncotype,
  MammaPrint = mammaprint,
  PAM50 = PAM50
)



# Prepare gene sets
gene_sets <- lapply(signatures, function(x) toupper(x))
names(gene_sets) <- names(signatures)

###################################################################################
# compute enrichment using UCell

set.seed(123)
sc_epithelium <- AddModuleScore_UCell(sc_epithelium, features = gene_sets, maxRank = 70, ncores = 10, force.gc = TRUE)
head(sc_epithelium[[]])
gc()

signature.names <- paste0(names(signatures), "_UCell")

png('./figures/revision/signature_enr_scRNAseq_HBCA.png', width = 2700, height = 2000, res = 250)
VlnPlot(sc_epithelium, 
        features = signature.names, 
        pt.size = 0,
        group.by = "cell_type", 
        same.y.lims = TRUE, 
)
dev.off()


VlnPlot(sc_epithelium, 
        features = signature.names, 
        pt.size = 0,
        group.by = "cell_type", 
        same.y.lims = TRUE
)

#####################
# signature score smoothing
sc_epithelium <- SmoothKNN(sc_epithelium, signature.names = signature.names, reduction="pca")

signature.names_smooth <- paste0(names(signatures), "_UCell_kNN")

png('./figures/revision/signature_enr_scRNAseq_HBCA_smooth.png', width = 2700, height = 2000, res = 250)
VlnPlot(sc_epithelium, 
        features = signature.names_smooth, 
        pt.size = 0,
        group.by = "cell_type", 
        same.y.lims = TRUE, 
)
dev.off()







# Calculate mean UCell scores for each signature by cell type
mean_scores <- sc_epithelium@meta.data %>%
  group_by(cell_type) %>%
  summarise(across(ends_with("_UCell"), mean, na.rm = TRUE))

library(dplyr)
library(ggpubr)

# Perform ANOVA and Tukey's HSD for each signature
results <- list()
for (signature in signature.names) {
  model <- aov(reformulate('cell_type', signature), data = sc_epithelium@meta.data)
  tukey <- TukeyHSD(model)
  results[[signature]] <- tukey
}
results


# Boxplot to visualize UCell scores by cell type for each signature
sc_epithelium_long <- sc_epithelium@meta.data %>%
  pivot_longer(
    cols = ends_with("_UCell"),
    names_to = "signature",
    values_to = "UCell_score"
  ) %>%
  mutate(signature = gsub("_UCell", "", signature)) # Clean up the signature names if needed

# Boxplot to visualize UCell scores by cell type for each signature
ggplot(sc_epithelium_long, aes(x = cell_type, y = UCell_score, fill = cell_type)) +
  geom_boxplot() +
  facet_wrap(~signature, scales = 'free_y') + # Facet by signature
  theme_bw() +
  labs(title = "Distribution of UCell Scores by Cell Type", y = "UCell Score", x = "Cell Type")


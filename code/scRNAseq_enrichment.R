
rm(list = ls())

gc()

library(Seurat)
library(AUCell)
library(pochi)
library(scMiko)
library(UCell)

##################
# load the data
sc_epithelium <- readRDS('data/scRNAseq/Kumar_epithelium.rds')
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


sc_epithelium <- RunPCA(sc_epithelium, features = VariableFeatures(object = sc_epithelium))
sc_epithelium <- FindNeighbors(sc_epithelium, dims = 1:10)
sc_epithelium <- RunUMAP(sc_epithelium, dims = 1:20)
#DimPlot(object = sc_epithelium, group.by = "author_cell_type", label = TRUE)

save(sc_epithelium, file = './data/scRNAseq/Kumar_epithelium_processed.rda')

# load
load('./data/scRNAseq/Kumar_epithelium_processed.rda')

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

signature.names <- paste0(names(signatures), "_UCell")

png('./figures/revision/signature_enr_scRNAseq_THR70_Kumar.png', width = 2700, height = 2000, res = 250)
VlnPlot(sc_epithelium, 
        features = signature.names, 
        pt.size = 0,
        group.by = "author_cell_type", 
        same.y.lims = TRUE, 
        )
dev.off()


##########################################
# Test if THR-70 is more enriched in LummHR cells compared to other cell types
##########################################
# Define luminal and non-luminal cell types
luminalHR_types <- c("LummHR-SCGB", "LummHR-active", "LummHR-major")
other_types <- setdiff(unique(sc_epithelium$author_cell_type), luminalHR_types)

# Create a new column to distinguish between luminal and non-luminal cell types
sc_epithelium$luminal_vs_other <- ifelse(sc_epithelium$author_cell_type %in% luminalHR_types, "Luminal-HR", "Other epithelial cells")

# Subset the data for relevant comparisons
luminal_scores <- FetchData(sc_epithelium, vars = c("THR70_UCell", "luminal_vs_other"))

# Perform Wilcoxon Rank-Sum test
test_result <- luminal_scores %>%
  group_by(luminal_vs_other) %>%
  summarise(p_value = wilcox.test(THR70_UCell ~ luminal_vs_other, data = .)$p.value,
            .groups = 'drop') # Drop the grouping

# Output the test results
print(test_result)

# Boxplot or Violin plot of UCell scores
ggplot(luminal_scores, aes(x = luminal_vs_other, y = THR70_UCell, fill = luminal_vs_other)) +
  geom_boxplot() +
  labs(title = "Comparison of THR-70 UCell Scores",
       x = "Group",
       y = "THR-70 UCell Score") +
  theme_minimal()


#####################
# signature score smoothing
sc_epithelium <- SmoothKNN(sc_epithelium, signature.names = signature.names, reduction="pca")

signature.names_smooth <- paste0(names(signatures), "_UCell_kNN")

png('./figures/revision/signature_enr_scRNAseq_Kumar_smooth.png', width = 2700, height = 2000, res = 250)
VlnPlot(sc_epithelium, 
        features = signature.names_smooth, 
        pt.size = 0,
        group.by = "author_cell_type", 
        same.y.lims = TRUE, 
)
dev.off()






# Calculate mean UCell scores for each signature by cell type
mean_scores <- sc_epithelium@meta.data %>%
  group_by(author_cell_type) %>%
  summarise(across(ends_with("_UCell"), mean, na.rm = TRUE))

library(dplyr)
library(ggpubr)

# Perform ANOVA and Tukey's HSD for each signature
results <- list()
for (signature in signature.names) {
  model <- aov(reformulate('author_cell_type', signature), data = sc_epithelium@meta.data)
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
ggplot(sc_epithelium_long, aes(x = author_cell_type, y = UCell_score, fill = author_cell_type)) +
  geom_boxplot() +
  facet_wrap(~signature, scales = 'free_y') + # Facet by signature
  theme_bw() +
  labs(title = "Distribution of UCell Scores by Cell Type", y = "UCell Score", x = "Cell Type")



###################################################################################
# AUC
###################################################################################
counts <- GetAssayData(object = sc_epithelium, layer = "counts")

cell_rankings <- AUCell_buildRankings(counts, nCores = 10)
gc()

cells_AUC <- AUCell_calcAUC(gene_sets, cell_rankings)
gc()

cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist = TRUE, assign=TRUE)


## Get assigned cells
cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name = "cell")
colnames(assignmentTable)[2] <- "geneSet"
head(assignmentTable)


assignmentMat <- table(assignmentTable[, "geneSet"], assignmentTable[, "cell"])
assignmentMat[, 1:2]

set.seed(123)
miniAssigMat <- assignmentMat[, sample(1:ncol(assignmentMat), min(500, ncol(assignmentMat)))]
pheatmap(miniAssigMat, cluster_rows = FALSE, show_colnames = FALSE)


assignments <- assignmentTable %>% dplyr::group_by(cell) %>% 
  dplyr::summarize(AUCellType = paste(geneSet, collapse = "/")) %>%
  dplyr::group_by(AUCellType) %>%
  dplyr::mutate(n = length(unique(cell))) %>%
  dplyr::filter(n > 5)
sce$AUCell_celltypes <- assignments$AUCellType[match(sce$cell, assignments$cell)]


plotTSNE(sce, colour_by = "AUCell_celltypes")


colData(sce) <- cbind(colData(sce), DataFrame(t(assay(cell_AUC, "AUC"))))
plotTSNE(sce, colour_by = "B.cells")


plotTSNE(sce, colour_by = "Monocytes")



plotTSNE(sce, colour_by = "NK.cells")



cells_assignment$geneSet$assignment



new_cells <- names(which(getAUC(cells_AUC)["geneSet",]>0.15))


sc_epithelium$Diagonosis <- ifelse(colnames(tiss) %in% new_cells, "Disease", "Normal")


DimPlot(object = sc_epithelium, group.by = "Diagonosis", label = TRUE)

DimPlot(object = sc_epithelium, group.by = "cell.types", label = TRUE)








# Calculate AUC scores
#cells_AUC <- AUCell_run(counts, signatures)
sc_epithelium <- RunAUCell(object = sc_epithelium, genesets = gene_sets, assay = "RNA", slot = "counts")

# View scores
head(sc_epithelium[["AUCell"]][["scores"]])



# Visualize AUC scores on a UMAP plot
sc_epithelium <- RunUMAP(sc_epithelium, dims = 1:20)
DimPlot(sc_epithelium, group.by = "THR50_AUC")  # replace 'THR50_AUC' with the actual name used in your AUCell scores

# Heatmap of scores across clusters
DoHeatmap(sc_epithelium, features = c("THR50_AUC", "THR25_AUC"))  # Adjust according to your AUC score names


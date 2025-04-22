library(data.table)
library(tidyverse)
library(optparse)

#TORUN
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/mergeDNA_RNA_diagnostics.R \
# --mergedPath /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/results/mergeRNADNAme/merged_repeatsDNAme_RNA_table.tsv \
# --metadata /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/results/mergeRNADNAme/merged_repeatsDNAme_RNA_table.tsv \
# --outDir merge_dna_rna_diagnostics

# Define command-line options
option_list <- list(
    make_option(c("-m", "--mergedPath"), type = "character", 
                help = "Path to the merged repeats DNA methylation and RNA table"),
make_option(c("-d", "--metadata"), type = "character", 
                help = "Path to the metadata file"),
    make_option(c("-o", "--outDir"), type = "character", default = ".", 
                help = "Output directory [default %default]")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# opt$mergedPath <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/promoterMethyl/results/mergeRNADNAme/merged_repeatsDNAme_RNA_table.tsv"
# opt$outDir <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/promoterMethyl/results"
# opt$metadata <- "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_recoded.csv"



# Check if mergedPath is provided
if (is.null(opt$mergedPath)) {
    stop("Error: --mergedPath argument is required.")
}
if (is.null(opt$metadata)) {
    stop("Error: --metadata argument is required.")
}

# Load the merged data
metadataPath <- opt$metadata  # Assign opt$metadata to metadataPath
metadataDT <- fread(metadataPath)
mergedDNAmeRNADT <- fread(opt$mergedPath)

# mergedPath <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/results/mergeRNADNAme/merged_repeatsDNAme_RNA_table.tsv"
# Load the merged data

# dt10 <- fread("/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_recoded.csv")
message("Loading merged data...")
message(metadataDT)

dataMerge <- mergedDNAmeRNADT[metadataDT[, c("condition","new_samples_name")], on=c("RNAsampleID"="new_samples_name")]

fwrite(dataMerge, file.path(opt$outDir, "merged_repeatsDNAme_RNA_table_WithConds.tsv"), sep="\t", row.names=FALSE)


splitDT <- split(dataMerge, dataMerge$TE_ID)

pdf(file.path(opt$outDir, "scatterFPKMvsDNAme.pdf"), width = 10, height = 6)
lapply(head(names(splitDT)), function(x) {
  scatterPlot <- ggplot(data = splitDT[[x]]) + labs(title = x) +
    geom_point(aes(x = fracMethyl.geom_mean, y = fpkm, color = condition), size = 4) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
})
dev.off()

# Create the boxplot
pltDNAme <- ggplot(dataMerge, aes(x = filename, y = fracMethyl.geom_mean, fill = filename)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(width = 0.2, alpha = 0.01) +  # Add jitter for individual points
  facet_wrap(~ condition, scales = "free_x") +  # Facet by condition
  labs(
    title = "Boxplot of fracMethyl.geom_mean by Sample",
    x = "Sample",
    y = "Geometric Mean of Fraction Methylation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(pltDNAme, filename = file.path(opt$outDir, "boxplot_fracMethyl_geom_mean_by_sample.png"), width = 10, height = 6)

# ggsave(pltDNAme, filename = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/boxplot_fracMethyl_geom_mean_by_sample.png", width = 10, height = 6)




# Define the order for consensus
L1Order <- c("L1Md_A", "L1Md_T", "L1Md_F", "L1Md_F2", "L1Md_F3", "L1Md_Gf")

# Filter data for the specified consensus values and set the factor levels
dataConsensus <- dataMerge[consensus %in% L1Order & fpkm >= 2, ]
dataConsensus$consensus <- factor(dataConsensus$consensus, levels = L1Order, ordered = TRUE)
dataConsensusL1MdDT <- split(dataConsensus, dataConsensus$TE_ID)

pdf(file.path(opt$outDir, "scatterFPKMvsDNAmeL1MD.pdf"), width = 10, height = 6)
lapply(names(dataConsensusL1MdDT), function(x) {
  scatterPlot <- ggplot(data = splitDT[[x]]) + labs(title = x) +
    geom_point(aes(x = fracMethyl.geom_mean, y = fpkm, color = condition), size = 4) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
})
dev.off()


# Create the boxplot
pltConsensus <- ggplot(dataConsensus, aes(x = consensus, y = fracMethyl.geom_mean, fill = consensus)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(width = 0.2, alpha = 0.01) +  # Add jitter for individual points
  facet_wrap(~ condition + filename, scales = "free_x") +  # Facet by condition
  labs(
    title = "Boxplot of fracMethyl.geom_mean by Consensus",
    x = "Consensus",
    y = "Geometric Mean of Fraction Methylation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Save the plot
# ggsave(pltConsensus, filename = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/boxplot_fracMethyl_geom_mean_by_consensus.png", width = 10, height = 6)




pltConsensusFacetCond <- ggplot(dataConsensus, aes(x = consensus, y = fracMethyl.geom_mean)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(width = 0.2, alpha = 0.01) +  # Add jitter for individual points
  facet_wrap(~ condition , scales = "free_x") +  # Facet by condition
  labs(
    title = "Boxplot of fracMethyl.geom_mean by Consensus",
    x = "Consensus",
    y = "Geometric Mean of Fraction Methylation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Save the plot
# ggsave(pltConsensusFacetCond, filename = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/boxplot_fracMethyl_geom_mean_by_consensusCond.png", width = 10, height = 6)




# Filter data for validReadCounts == TRUE
dataScatter <- dataConsensus[validReadCounts == TRUE]

# Create the scatter plot
pltScatter <- ggplot(dataScatter, aes(x = fpkm, y = fracMethyl.geom_mean)) +
  geom_point(alpha = 0.7, size = 0.01) +  # Scatter points
    facet_wrap(~ condition + filename, scales = "free") +
  labs(
    title = "Scatter Plot of FPKM vs Geometric Mean of Fraction Methylation",
    x = "FPKM",
    y = "Geometric Mean of Fraction Methylation"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 10),
    legend.position = "right"
  )

# Save the plot
# ggsave(pltScatter, filename = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/scatter_fpkm_vs_fracMethyl_geom_mean.png", width = 10, height = 6)

# Save the boxplot for consensus
ggsave(pltConsensus, 
       filename = file.path(opt$outDir, "boxplot_fracMethyl_geom_mean_by_consensus.png"), 
       width = 10, height = 6)

# Save the boxplot for consensus faceted by condition
ggsave(pltConsensusFacetCond, 
       filename = file.path(opt$outDir, "boxplot_fracMethyl_geom_mean_by_consensusCond.png"), 
       width = 10, height = 6)

# Save the scatter plot
ggsave(pltScatter, 
       filename = file.path(opt$outDir, "scatter_fpkm_vs_fracMethyl_geom_mean.png"), 
       width = 10, height = 6)
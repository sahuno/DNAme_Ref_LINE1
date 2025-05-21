#
# merge DNA methylation rates with RNAseq counts
# mkdir -p DNAme_RNA_rln
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/DE_RepeatsHeatmaps.R


library(data.table)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(UpSetR)
library(optparse)
library(pvclust)
library(ggnewscale)
library(viridisLite)
library(cola)
library(knitr)

if(!require('NbClust')) {
  install.packages('NbClust')
  library('NbClust')
}


# library(cowplot)
# library(RColorBrewer)
# library(NbClust)
# library(GenomicRanges)
# library(gUtils)
# library(janitor)
# library(fst)
# library(ggrepel)
# library(ggfortify)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(EnhancedVolcano)
# library(viridis)
# library(magrittr)



cividis_palette <- viridis(100, option = "cividis")
options("width"=200)
# set_fst_threads <- 4
# threads_fst(set_fst_threads)

# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/DE_RepeatsHeatmaps.R \
# --lfcCutoff 0.05 \
# --l1Separately FALSE \
# --paths_rdata "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/LocusSpecific" \
# --FulllengthAnno "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/overlaps75SquireLINE_vs_L1_L1Base_wa.bed"

option_list <- list(
    make_option(c("-l", "--lfcCutoff"), type="numeric", default=0.05, help="Cutoff for adjusted p-value [default %default]"),
    make_option(c("-s", "--l1Separately"), type="logical", default=FALSE, help="run locus specific line-1 separately [default %default]"),
    make_option(c("-d", "--paths_rdata"), type="character", default=NULL,help="Path to Rdata folder"),
    make_option(c("-f", "--FulllengthAnno"), type="character", default=NULL, help="Path to full-length LINE-1 annotation file")
    )

# default="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecific"
parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

# opt$paths_rdata <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/LocusSpecific"
# opt$FulllengthAnno <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/overlaps75SquireLINE_vs_L1_L1Base_wa.bed"

# After parsing, split comma-separated values into vectors
paths_rdata <- opt$paths_rdata


## Therapy list
combinedTherapy <- c("SETDB1i-CKi_SETDB1i" ,"SETDB1i-CKi_CKi","QSTAT-CKi_QSTAT" ,"QSTAT-CKi_CKi","SETDB1i-CKi_QSTAT-CKi")
monTherapy <- c("AZA_DMSO", "CKi_DMSO","QSTAT_DMSO", "QSTAT-CKi_DMSO", "SETDB1i_DMSO", "SETDB1i-CKi_DMSO")

#customize the therapy list
# subSampleConds <-  c("DMSO_AZA", "DMSO_CKi", "DMSO_QSTAT", "DMSO_SETDB1i" , "DMSO_QSTAT-CKi" , "DMSO_SETDB1i-CKi","SETDB1i_SETDB1i-CKi", "QSTAT_QSTAT-CKi")
subSampleConds <-  c("CKi_DMSO", "QSTAT_DMSO", "QSTAT-CKi_DMSO")
# subSampleConds <-  c("AZA_DMSO", "CKi_DMSO", "QSTAT_DMSO", "SETDB1i_DMSO", "QSTAT-CKi_DMSO","SETDB1i-CKi_DMSO")

# subSampleConds <-  c("AZA_DMSO", "CKi_DMSO", "QSTAT_DMSO", "SETDB1i_DMSO", "QSTAT-CKi_DMSO","SETDB1i-CKi_DMSO","SETDB1i-CKi_SETDB1i", "QSTAT-CKi_QSTAT")

pat <- paste(subSampleConds, collapse="|")



if(opt$l1Separately){
  message("running locus specific LINE-1 separately")

### load RNAseq data/ full length active L1s
filesRDATAL1 <- list.files(opt$paths_rdata, pattern = "Differential_locusSpecific_L1ExpressionResults", recursive = TRUE, full.names = TRUE) 


filesRDATAL1 <- str_subset(filesRDATAL1, pat)


# load(filesRDATA[1])

loadedRDataL1 <- lapply(filesRDATAL1, function(x) {
    e <- new.env()
    load(x, envir = e)
    as.list(e)
  })


##rename to list 
baseNamesL1 <- gsub(".*condition_", "", basename(filesRDATAL1))
names(loadedRDataL1) <- gsub(".RData", "", baseNamesL1)

resultsDTL1 <- lapply(names(loadedRDataL1), function(x){
  res <- loadedRDataL1[[x]][["results2SaveLspecL1"]]$resDF
  res <- as.data.frame(res)
  res$RepID <- row.names(res)
  return(res)
})
names(resultsDTL1) <- names(loadedRDataL1)


message("binding results")
resultsDTL1Sample <- rbindlist(resultsDTL1, idcol = "condition")

resultsDTL1Sample <- resultsDTL1Sample %>% mutate(LFClabel = case_when(log2FoldChange > 0 & padj <= opt$lfcCutoff ~ "up",
                                     log2FoldChange < 0 & padj <= opt$lfcCutoff ~ "down",
                                     TRUE ~  "no_change")) 

resultsDTL1Sample <- resultsDTL1Sample %>% filter(condition %in% subSampleConds) #%>% group_by(condition, LFClabel) %>% summarise(count = n(), .groups = "drop") %>% arrange(-count)



resultsDTL1Sample$log2FoldChange[is.na(resultsDTL1Sample$log2FoldChange)] <- 0

resultsDTL1Sample <- resultsDTL1Sample %>% mutate(ClassFamily = str_split_i(RepID, "\\|", 4)) %>% separate(ClassFamily, into = c("subfamily", "family", "class"), sep = ":")
# table(resultsDTSample$class)

fwrite(resultsDTL1Sample, file = file.path("resultsL1LocusSpecificRNA_DE_withMetadata.txt"), sep = "\t")


# Use tidyr to reshape data: pivot so rows are LINE1 elements and columns are samples
heatmap_dataL1 <- resultsDTL1Sample %>%
  select(condition, RepID, log2FoldChange) %>%
  pivot_wider(names_from = condition, values_from = log2FoldChange, values_fill = 0)

# Set the row names to the LINE1 element IDs (RepID) and remove the RepID column from the data frame
heatmap_matrixL1 <- as.data.frame(heatmap_dataL1) %>% 
  column_to_rownames(var = "RepID") %>%
  as.matrix()

##turn thi
heatmap_matrixL1 <- heatmap_matrixL1[, subSampleConds]
fwrite(heatmap_matrixL1, file = file.path("heatmap_matrixL1.txt"), sep = "\t")


result <- pvclust(heatmap_matrixL1, method.hclust="ward.D2", method.dist="euclidean", nboot=1000)
# Plot the dendrogram with AU/BP values
pdf("pvclust_LocusSpecifcL1_dendrogram.pdf", width = 10, height = 6)
plot(result)
pvrect(result, alpha=0.95)
dev.off()

pdf("pvclust_Seplot_LocusSpecifcL1_dendrogram.pdf", width = 10, height = 6)
seplot(result, identify=TRUE)
dev.off()

# run the standard SE plot
# seplot(result, type="au", identify=FALSE,
#        main="SE of AU p-values", xlab="Relative bootstrap sample size",
#        ylab="Standard Error")

#####################

# make a dendrogram and convert to phylo for ggtree
# library(pvclust)
# library(dendextend)
# library(ggtree)

# dend <- as.dendrogram(result$hclust)
# phy  <- ape::as.phylo(dend)

# # plot with ggtree
# p <- ggtree(phy)

# # now annotate nodes with AU bootstrap values
# # pvclust stores AU in result$edges[, "au"]
# au_vals <- data.frame(node = result$hclust$merge %>% 
#                         # get node numbers corresponding to merges
#                         { seq_along(apply(result$hclust$merge,1,sum)) + length(result$hclust$labels) },
#                       au   = result$edges[,"au"])

# p + geom_text2(data=au_vals, aes(subset = TRUE, label=round(au,1)),
#                hjust=-.3, size=3)



#####################
pheatmap(heatmap_matrixL1,
         cluster_rows = TRUE,    # Cluster LINE1 elements
         cluster_cols = TRUE,    # Cluster samples
         scale = "row",          # Scale each row (optional: remove if raw values are preferred)
         show_rownames = FALSE,   # Whether to display row names (can be set to FALSE if too many rows)
         show_colnames = TRUE,   # Display column names
         filename = "heatmap_line1.png",
         color = cividis_palette,
         main = "Heatmap of LINE1 Elements (log2FoldChange)")


#plot lines 
annotation_row <- data.frame(Lineage = sapply(rownames(heatmap_matrixL1), function(x) {
  # First split by "|"
  components <- strsplit(x, split = "\\|")[[1]]
  if (length(components) >= 4) {
    # Now split the 4th component by ":"
    subcomponents <- strsplit(components[4], split = ":")[[1]]
    # Return the 2nd element if available, else NA
    if (length(subcomponents) >= 2) {
      return(subcomponents[1])
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}))
rownames(annotation_row) <- rownames(heatmap_matrixL1)

pheatmap(heatmap_matrixL1[, subSampleConds],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         annotation_row = annotation_row,
         show_rownames = FALSE,
         color = cividis_palette,
         filename = "heatmap_line1_withLineages.png",
         main = "Heatmap of LINE1 Elements")


# table(annotation_row$Lineage)

# "CKi_QSTAT-CKi","CKi_SETDB1i-CKi",,"QSTAT_SETDB1i"         "QSTAT-CKi_SETDB1i-CKi" 

df_filteredL1 <- resultsDTL1Sample %>%
  filter(LFClabel %in% c("up", "down")) %>%
  group_by(condition, LFClabel) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(count = ifelse(LFClabel == "down", -count, count))  # Make "down" counts negative


df_orderedL1 <- df_filteredL1 %>%
  filter(LFClabel == "up") %>%
  arrange(desc(count)) %>%
  pull(condition) %>%
  unique()

df_filteredL1$condition <- factor(df_filteredL1$condition, levels = subSampleConds)


plt <- ggplot(df_filteredL1, aes(x = condition, y = count, fill = LFClabel)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "Up and Down Sig. LINE1 per condition", x = "condition", y = "Count", fill = "Label") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plt, file = paste0("LINE1up_down_counts_per_Condition.png"), width = 10, height = 6)

} else {
  message("skipping locus specific LINE-1 separately")
}






###############################################################
###### Repeat analysis for all repeats
###############################################################
filesRDATAReps <- list.files(opt$paths_rdata, pattern = "Differential_locusSpecific_RepeatsExpressionResults", recursive = TRUE, full.names = TRUE) 

filesRDATAReps <- str_subset(filesRDATAReps, pat)

loadedRDataReps <- lapply(filesRDATAReps, function(x) {
    e <- new.env()
    load(x, envir = e)
    as.list(e)
  })


##rename to list 
# loadedRDataReps[[1]][["results2SaveLspecRepeats"]]
resultsDT <- lapply(loadedRDataReps, function(x){
  res <- x[["results2SaveLspecRepeats"]]$resDF
  res <- as.data.frame(res)
  res$RepID <- row.names(res)
  return(res)
})
names(resultsDT) <- gsub(".RData", "", gsub(".*condition_", "", filesRDATAReps))


resultsDTSample <- rbindlist(resultsDT, idcol = "condition")

# Check if 'fullLengthpath' is a column in resultsDTSample
if ("isfullLengthL1" %in% colnames(resultsDTSample)) {
  message("'isfullLengthL1' column exists in resultsDTSample.")
} else {
  message("'isfullLengthL1' column does not exist in resultsDTSample.")
  message("adding isfullLengthL1 column to resultsDTSample.")

if (file.exists(opt$FulllengthAnno) && !file.info(opt$FulllengthAnno)$isdir) {
  message(paste(opt$FulllengthAnno, "exists and is a file."))
colInterest <- "V11" # this will cause an error
fullDT <- fread(opt$FulllengthAnno, select = colInterest)
# if any entry in fullDT[[colInterest]] matched resultsDTSample$RepID create a new column `isfullLengthL1` in resultsDTSample with TRUE/FALSE
resultsDTSample[, isfullLengthL1 := RepID %in% fullDT[[colInterest]]]
} else {
  message(paste(opt$FulllengthAnno, "does not exist or is not a file."))
}
}

# fullLengthpath
# resultsDTSample[isfullLengthL1==TRUE,]

#Assign labels based on log2FoldChange and adjusted p-value
resultsDTSample <- resultsDTSample %>% mutate(LFClabel = case_when(log2FoldChange > 0 & padj <= opt$lfcCutoff ~ "up",
                                     log2FoldChange < 0 & padj <= opt$lfcCutoff ~ "down",
                                     TRUE ~  "no_change")) 


fwrite(resultsDTSample, file = file.path("resultsAllRepeatsRNA_DE_withMetadata.txt"), sep = "\t")

# awk 'NR==1 || $10 == "TRUE"' resultsAllRepeatsRNA_DE_withMetadata.txt > resultsFromAllRepeatsRNA_DE_FullLengthL1_withMetadata.txt

################################################################################################
######## Repeats up and down counts ###########################################################
################################################################################################
resultsDTSample <- resultsDTSample %>% filter(LFClabel %in% c("up", "down")) %>% 
mutate(ClassFamily = str_split_i(RepID, "\\|", 4)) %>%
separate(ClassFamily, into = c("subfamily", "family", "class"), sep = ":")


head(resultsDTSample)
resultsDTSampleUpDown <- resultsDTSample %>% filter(LFClabel %in% c("up", "down"))

# resultsDTSampleUpDown <- 

resultsDTSampleUpDown %>% filter(str_detect(condition, "SETDB1i")) %>% head() #(condition, LFClabel) %>% summarise(count = n(), .groups = "drop") %>% arrange(-count)

pltUPDown <- ggplot(resultsDTSampleUpDown, aes(x = condition, y = log2FoldChange)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "Up and Down Sig. Repeats per Condition",
       x = "Condition",
       y = "log2FoldChange",
       color = "Label") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~ class, scales = "free_y")    # Facet by repeat class
# ggsave(pltUPDown, file = "RepeatsUp_down_LFC_per_Condition.png", width = 10, height = 6)


library(tidyplots)
pltUPDown <- resultsDTSampleUpDown %>% tidyplot(x = condition, y = log2FoldChange) |> 
    add_boxplot(alpha = 0.4) |> 
  # add_sem_errorbar() |> 
  adjust_x_axis(rotate_labels = 45)|>
  split_plot(by = class, ncol = 2, nrow = 2) |>
  save_plot(filename = "RepeatsUp_down_LFC_per_Condition.pdf")
  # add_data_points_beeswarm() |> 
  # geom_hline(yintercept = 0, color = "black") +
  # labs(title = "Up and Down Sig. Repeats per Condition",
  #      x = "Condition",
  #      y = "log2FoldChange",
  #      color = "Label") +
  # theme_minimal() +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) |>  facet_wrap(vars(class), scales = "free_y")    # Facet by repeat class
# ggsave(pltUPDown, file = "RepeatsUp_down_LFC_per_Condition.png", width = 10, height = 6)
# ncols = 2, nrows = 2, width = 10, height = 6


# group_by(condition, LFClabel) %>% summarise(count = n(), .groups = "drop") %>% arrange(-count)

fullLengthLoc1SETDB1 <- resultsDTSample %>% filter(RepID == "chr17|13838024|13845083|L1Md_T:L1:LINE|32|-") 

plt_chr17_13838024_13845083_L1Md_T <- ggplot(fullLengthLoc1SETDB1, aes(x = condition, y = log2FoldChange, color = LFClabel)) +
  geom_point(size = 8) +
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "Full Length (chr17|13838024|13845083|L1Md_T:L1:LINE|32|-) SETDB1i Group",
       x = "Condition",
       y = "log2FoldChange",
       color = "Label") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plt_chr17_13838024_13845083_L1Md_T, file = "SETDB1i_chr17_13838024_13845083_L1Md_T.png", width = 10, height = 6)

df_filteredCondClass <- resultsDTSample %>% select(condition, RepID, class, family, log2FoldChange, padj, LFClabel) %>%
  arrange(-log2FoldChange) %>%
  group_by(condition, class, LFClabel) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(count = ifelse(LFClabel == "down", -count, count))

pltClass <- ggplot(df_filteredCondClass, aes(x = condition, y = count, fill = LFClabel)) +
  geom_bar(stat = "identity", position = "identity") +
  facet_wrap(~ class, scales = "free_y") +   # Facet by repeat class
  scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "Up and Down Sig. Repeats per Condition by Class",
       x = "Condition",
       y = "Count",
       fill = "Label") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(pltClass, file = "RepeatsUp_down_countsGroups_per_Condition_Class.png", width = 10, height = 6)


# chr17|17554111|17560255|L1Md_T:L1:LINE|35|+
# fullLengthLoc1SETDB1 <- resultsDTSample %>% filter(RepID == "chr17|17554111|17560255|L1Md_T:L1:LINE|35|+") 

plt_qstat <- ggplot(resultsDTSample %>% filter(RepID == "chr17|17554111|17560255|L1Md_T:L1:LINE|35|+"), aes(x = condition, y = log2FoldChange, color = LFClabel)) +
  geom_point(size = 8) +
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "Full Length (chr17|17554111|17560255|L1Md_T:L1:LINE|35|+) QSTATi Group",
       x = "Condition",
       y = "log2FoldChange",
       color = "Label") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plt_qstat, file = "QstatGroup_chr17_17554111_17560255_L1Md_T.png", width = 10, height = 6)



resultsDTSample %>% filter(isfullLengthL1 == TRUE) %>% group_by(condition, LFClabel) %>% summarise(count = n(), .groups = "drop") %>% arrange(-count)
# resultsDTSample$log2FoldChange[is.na(resultsDTSample$log2FoldChange)] <- 0

# Use tidyr to reshape data: pivot so rows are LINE1 elements and columns are samples
heatmap_data <- resultsDTSample %>%
  select(condition, RepID, log2FoldChange) %>%
  pivot_wider(names_from = condition, values_from = log2FoldChange, values_fill = 0)

# Set the row names to the LINE1 element IDs (RepID) and remove the RepID column from the data frame
heatmap_matrix <- as.data.frame(heatmap_data) %>% 
  column_to_rownames(var = "RepID") %>%
  as.matrix()

heatmap_matrixAdjusted <- adjust_matrix(heatmap_matrix, sd_quantile = 0.05, max_na = 0.25, verbose = TRUE)

heatmap_matrixAdjusted_t = t(heatmap_matrixAdjusted)

# res <- NbClust(heatmap_matrixAdjusted, distance = "euclidean", min.nc=2, max.nc=8, 
#             method = "complete", index = "ch")



subSampleCondsDF <-  data.frame(conditions = subSampleConds)
DrugsAnnot <- subSampleCondsDF %>% dplyr::filter(conditions %in% subSampleConds) %>% 
mutate(Combinations = case_when(str_detect(conditions, "-") ~ "combo",TRUE ~ "mono"), 
BroadTargets = case_when(str_detect(conditions, "QSTAT_|SETDB1i_") ~ "Chromatin",
                          str_detect(conditions, "-CKi_") ~ "Chromatin+MEK",
                          str_detect(conditions, "AZA_") ~ "DNAme",
                          str_detect(conditions, "CKi_DMSO") ~ "MEK", TRUE ~ NA),
Action = case_when(str_detect(conditions, "QSTAT_") ~ "HDACi", 
                    str_detect(conditions, "SETDB1i_") ~ "H3K9mei",
                    str_detect(conditions, "QSTAT-CKi_DMSO") ~ "HDACi+MEKi",
                    str_detect(conditions, "SETDB1i-CKi_DMSO") ~ "H3K9mei+MEKi",
                    str_detect(conditions, "AZA_") ~ "DNMTi",
                    str_detect(conditions, "CKi_DMSO") ~ "MEKi", TRUE ~ NA))

DrugsAnnot <- DrugsAnnot %>% mutate(MEKi = case_when(str_detect(Action, "MEKi") ~ "Yes", TRUE ~ "No"))
# )
rownames(DrugsAnnot) <- DrugsAnnot$conditions

k <- pheatmap::pheatmap(heatmap_matrixAdjusted, kmeans_k = 4, annotation = DrugsAnnot[,c("Combinations", "Action", "BroadTargets", "MEKi")], 
# colorRampPalette(c("#5D3A9B", "white", "#E66100"))(100),
main = paste0("abs(LFC) > 1 & padj < 0.05; repeats [NA=0]"), 
scale = "row", 
annotation_names_row=FALSE, 
show_rownames = F, 
filename = paste0("Kmeans_clustering_sig_repeats_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.png"))

# names(k$kmeans)
clusterDF <- as.data.frame(factor(k$kmeans$cluster))
colnames(clusterDF) <- "Cluster"
OrderByCluster <- heatmap_matrixAdjusted[order(clusterDF$Cluster), ]


clusterNum <- NbClust(heatmap_matrixAdjusted, distance = "euclidean", min.nc = 2, max.nc = 12, method = "kmeans", index = "silhouette")
message("optimal number of clusters: ", clusterNum$Best.nc[1], "; with value index: ", clusterNum$Best.nc[2])

orderedCluster <- sort(clusterNum$Best.partition)

heatmap_matrixAdjusted_sorted_NBClust <- heatmap_matrixAdjusted[match(names(orderedCluster), rownames(heatmap_matrixAdjusted)), ]


# colnames(matrixMonoCombo_DMSO_ref_cleaned)

pheatmap::pheatmap(heatmap_matrixAdjusted_sorted_NBClust, 
annotation_row = clusterDF,
scale = "row", show_rownames = FALSE,
    cluster_rows = FALSE,
    annotation = DrugsAnnot[,c("Combinations", "Action", "BroadTargets", "MEKi")], 
    main = paste0("abs(LFC) > 1 & padj < 0.05; repeats [NA=0]"),
    filename = paste0("kmeans_withOptimumClustering_sig_Repeats_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.pdf"))



# Create annotation for treatment
# sample_anno = data.frame(
#   treatments = colnames(heatmap_matrixAdjusted))

# # Add row names to match the sample names
# rownames(sample_anno) = colnames(heatmap_matrixAdjusted)
####################################################################################################################
############################################### Cola consensus clustering ##########################################
####################################################################################################################

################# Begin of subset ################################
# set.seed(42)  # for reproducibility
# subset_matrix <- heatmap_matrixAdjusted[sample(nrow(heatmap_matrixAdjusted), 1000), ]
# feature_ids <- rownames(subset_matrix)
# feature_chr <- sapply(strsplit(feature_ids, "\\|"), `[`, 1)
# feature_names <- sapply(strsplit(feature_ids, "\\|"), `[`, 4)
# feature_class <- sapply(strsplit(feature_names, "\\:"), `[`, 3)
# feature_family <- sapply(strsplit(feature_names, "\\:"), `[`, 2)


# feature_anno <- data.frame(
#   family = feature_class,
#   row.names = feature_ids
# )
# # transpose so that columns = features
# subset_matrix_t <- t(subset_matrix)

# sanity chec
# Subsample 1000 rows
# # subset_matrix <- heatmap_matrixAdjusted[sample(nrow(heatmap_matrixAdjusted), 1000), ]
# message("running cola")
# stopifnot(all(colnames(subset_matrix_t) == rownames(feature_anno)))
# rh = hierarchical_partition(subset_matrix_t, anno= feature_anno, mc.cores = 24, top_n = 500, help = FALSE)
################# End of subset ################################

# feature_ids <- rownames(heatmap_matrixAdjusted)
# feature_chr <- sapply(strsplit(feature_ids, "\\|"), `[`, 1)
# feature_names <- sapply(strsplit(feature_ids, "\\|"), `[`, 4)
# feature_class <- sapply(strsplit(feature_names, "\\:"), `[`, 3)
# feature_family <- sapply(strsplit(feature_names, "\\:"), `[`, 2)

# feature_anno <- data.frame(
#   family = feature_class,
#   row.names = feature_ids
# )
# stopifnot(all(colnames(heatmap_matrixAdjusted_t) == rownames(feature_anno)))
# rh = hierarchical_partition(heatmap_matrixAdjusted_t,
# anno= feature_anno,
# partition_method     = "kmeans",
# min_samples  = 100,
# mc.cores = 8,
# top_n = 500, 
# help = FALSE)

# saveRDS(rh, file = "cola_hclust.rds")
# # set all subsequent figures to be larger and higher-res
# opts_chunk$set(
#   fig.width  = 8,    # width in inches
#   fig.height = 8,    # height in inches
#   dpi        = 300   # dots per inch
# )
# cola_report(rh, "cola_report_basic")
# cola_report(rh)
##########################################################################
##########################################################################
##########################################################################

# ##turn thi
# heatmap_matrix <- heatmap_matrix[, subSampleConds]
# clusterResultReps <- pvclust(heatmap_matrix, method.hclust="ward.D2", method.dist="euclidean", nboot=1000)
# # Plot the dendrogram with AU/BP values
# pdf("pvclust_LocusSpecifcRepeats_dendrogram.pdf", width = 10, height = 6)
# plot(clusterResultReps)           # Dendrogram with AU/BP values shown
# pvrect(result, alpha=0.95)
# dev.off()


# # Calculate the number of unique groups in the annotation
# annotation_row <- data.frame(Lineage = sapply(rownames(heatmap_matrix), function(x) {
#   # First split by "|"
#   components <- strsplit(x, split = "\\|")[[1]]
#   if (length(components) >= 4) {
#     # Now split the 4th component by ":"
#     subcomponents <- strsplit(components[4], split = ":")[[1]]
#     # Return the 2nd element if available, else NA
#     if (length(subcomponents) >= 2) {
#       return(subcomponents[3])
#     } else {
#       return(NA)
#     }
#   } else {
#     return(NA)
#   }
# }))
# rownames(annotation_row) <- rownames(heatmap_matrix)

# # Count the number of unique groups
# num_groups <- length(unique(annotation_row$Lineage[!is.na(annotation_row$Lineage)]))


# message("heatmap of repeats")
# # Add the number of groups to the heatmap title
# mat <- heatmap_matrix[, subSampleConds]
# row_sd <- apply(mat, 1, sd, na.rm = TRUE)
# mat_jit <- mat
# mat_jit[row_sd == 0, ] <- mat_jit[row_sd == 0, ] + rnorm(sum(row_sd==0)*ncol(mat), sd=1e-6)
# # then call pheatmap(mat_jit, …)


# pheatmap(mat_jit[, subSampleConds],
#          cluster_rows = TRUE,
#          cluster_cols = TRUE,
#          scale = "row",
#          annotation_row = annotation_row,
#          show_rownames = FALSE,
#          color = cividis_palette,
#          filename = "heatmap_Repeats.png",
#          main = paste("Heatmap Repeat Elements (", num_groups, " groups)", sep = ""))


# # Filter out rows where Lineage contains "SINE?", "DNA?", or "LTR?"
# lineSineLTR <- str_detect(annotation_row$Lineage, "SINE|LINE|LTR|DNA") & !str_detect(annotation_row$Lineage, "LINE.|SINE.|DNA.|LTR.")
# # Subset the annotation and heatmap matrix
# annotation_rowNew <- annotation_row[lineSineLTR, , drop = FALSE]
# heatMap_mat <- mat_jit[lineSineLTR, subSampleConds]
# table(annotation_rowNew$Lineage)

# # Generate the heatmap
# pheatmap(heatMap_mat,
#          cluster_rows = TRUE,
#          cluster_cols = TRUE,
#          scale = "row",
#          annotation_row = annotation_rowNew,
#          show_rownames = FALSE,
#          color = cividis_palette,
#          filename = "heatmap_Repeats_LINE_SINE_LTR_DNA.png",
#          main = "Heatmap Repeat Elements")



message("done running DE_RepeatsHeatmaps.R")
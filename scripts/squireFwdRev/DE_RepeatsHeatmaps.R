#
# merge DNA methylation rates with RNAseq counts
# mkdir -p DNAme_RNA_rln
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/mergeDNAme_and_RNASeq.R


library(data.table)
library(DESeq2)
library(tidyverse)
# library(EnhancedVolcano)
# library(viridis)
# library(magrittr)
library(pheatmap)
# library(ggrepel)
# library(ggfortify)
# library(clusterProfiler)
# library(org.Hs.eg.db)
library(ggnewscale)
# library(cowplot)
# library(RColorBrewer)
# library(NbClust)
# library(GenomicRanges)
# library(gUtils)
# library(janitor)
# library(fst)
library(UpSetR)
library(optparse)

library(viridisLite)

cividis_palette <- viridis(100, option = "cividis")



options("width"=200)
# set_fst_threads <- 4
# threads_fst(set_fst_threads)



option_list <- list(
    make_option(c("-l", "--lfcCutoff"), type="numeric", default=0.5, help="Cutoff for adjusted p-value [default %default]"),
    make_option(c("-d", "--paths_rdata"), type="character", default="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/LocusSpecific",help="Path to Rdata folder"))

parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

# After parsing, split comma-separated values into vectors
paths_rdata <- opt$paths_rdata

## Therapy list
combinedTherapy <- c("SETDB1i-CKi_SETDB1i" ,"SETDB1i-CKi_CKi","QSTAT-CKi_QSTAT" ,"QSTAT-CKi_CKi","SETDB1i-CKi_QSTAT-CKi")
monTherapy <- c("AZA_DMSO", "CKi_DMSO","QSTAT_DMSO", "QSTAT-CKi_DMSO", "SETDB1i_DMSO", "SETDB1i-CKi_DMSO")

#customize the therapy list
# subSampleConds <-  c("DMSO_AZA", "DMSO_CKi", "DMSO_QSTAT", "DMSO_SETDB1i" , "DMSO_QSTAT-CKi" , "DMSO_SETDB1i-CKi","SETDB1i_SETDB1i-CKi", "QSTAT_QSTAT-CKi")
subSampleConds <-  c("AZA_DMSO", "CKi_DMSO", "QSTAT_DMSO", "SETDB1i_DMSO", "QSTAT-CKi_DMSO","SETDB1i-CKi_DMSO","SETDB1i-CKi_SETDB1i", "QSTAT-CKi_QSTAT")

### load RNAseq data/ full length active L1s
filesRDATA <- list.files(opt$paths_rdata, pattern = "Differential_locusSpecific_L1ExpressionResults", recursive = TRUE, full.names = TRUE) 


# load(filesRDATA[1])

loadedRData <- lapply(filesRDATA, function(x) {
    e <- new.env()
    load(x, envir = e)
    as.list(e)
  })


##rename to list 
baseNames <- gsub(".*condition_", "", basename(filesRDATA))
names(loadedRData) <- gsub(".RData", "", baseNames)

# names(loadedRData[[1]])
# loadedRData[[1]][["results2SaveLspecL1"]]$resDF

#subSample list
# loadedRData <- loadedRData[subSampleConds]
# loadedRData[["AZA_DMSO"]][["results2SaveLspecL1"]]


resultsDT <- lapply(names(loadedRData), function(x){
  res <- loadedRData[[x]][["results2SaveLspecL1"]]$resDF
  res <- as.data.frame(res)
  res$RepID <- row.names(res)
  return(res)
})
names(resultsDT) <- names(loadedRData)



resultsDTSample <- rbindlist(resultsDT, idcol = "condition")

resultsDTSample <- resultsDTSample %>% mutate(LFClabel = case_when(log2FoldChange > 0 & padj <= opt$lfcCutoff ~ "up",
                                     log2FoldChange < 0 & padj <= opt$lfcCutoff ~ "down",
                                     TRUE ~  "no_change")) 



resultsDTSample$log2FoldChange[is.na(resultsDTSample$log2FoldChange)] <- 0

# Use tidyr to reshape data: pivot so rows are LINE1 elements and columns are samples
heatmap_data <- resultsDTSample %>%
  select(condition, RepID, log2FoldChange) %>%
  pivot_wider(names_from = condition, values_from = log2FoldChange, values_fill = 0)

# Set the row names to the LINE1 element IDs (RepID) and remove the RepID column from the data frame
heatmap_matrix <- as.data.frame(heatmap_data) %>% 
  column_to_rownames(var = "RepID") %>%
  as.matrix()

##turn thi
heatmap_matrix <- heatmap_matrix[, subSampleConds]


pheatmap(heatmap_matrix,
         cluster_rows = TRUE,    # Cluster LINE1 elements
         cluster_cols = TRUE,    # Cluster samples
         scale = "row",          # Scale each row (optional: remove if raw values are preferred)
         show_rownames = FALSE,   # Whether to display row names (can be set to FALSE if too many rows)
         show_colnames = TRUE,   # Display column names
         filename = "heatmap_line1.png",
         color = cividis_palette,
         main = "Heatmap of LINE1 Elements (log2FoldChange)")


#plot lines 
annotation_row <- data.frame(Lineage = sapply(rownames(heatmap_matrix), function(x) {
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
rownames(annotation_row) <- rownames(heatmap_matrix)

pheatmap(heatmap_matrix[, subSampleConds],
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

df_filtered <- resultsDTSample %>%
  filter(LFClabel %in% c("up", "down")) %>%
  group_by(condition, LFClabel) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(count = ifelse(LFClabel == "down", -count, count))  # Make "down" counts negative


df_ordered <- df_filtered %>%
  filter(LFClabel == "up") %>%
  arrange(desc(count)) %>%
  pull(condition) %>%
  unique()

df_filtered$condition <- factor(df_filtered$condition, levels = subSampleConds)


plt <- ggplot(df_filtered, aes(x = condition, y = count, fill = LFClabel)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "Up and Down Sig. LINE1 per condition", x = "condition", y = "Count", fill = "Label") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plt, file = paste0("LINE1up_down_counts_per_Condition.png"), width = 10, height = 6)


###############################################################
###### Repeat analysis for all repeats
filesRDATAReps <- list.files(opt$paths_rdata, pattern = "Differential_locusSpecific_RepeatsExpressionResults", recursive = TRUE, full.names = TRUE) 

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


#Assign labels based on log2FoldChange and adjusted p-value
resultsDTSample <- resultsDTSample %>% mutate(LFClabel = case_when(log2FoldChange > 0 & padj <= opt$lfcCutoff ~ "up",
                                     log2FoldChange < 0 & padj <= opt$lfcCutoff ~ "down",
                                     TRUE ~  "no_change")) 



resultsDTSample$log2FoldChange[is.na(resultsDTSample$log2FoldChange)] <- 0

# Use tidyr to reshape data: pivot so rows are LINE1 elements and columns are samples
heatmap_data <- resultsDTSample %>%
  select(condition, RepID, log2FoldChange) %>%
  pivot_wider(names_from = condition, values_from = log2FoldChange, values_fill = 0)

# Set the row names to the LINE1 element IDs (RepID) and remove the RepID column from the data frame
heatmap_matrix <- as.data.frame(heatmap_data) %>% 
  column_to_rownames(var = "RepID") %>%
  as.matrix()

##turn thi
heatmap_matrix <- heatmap_matrix[, subSampleConds]

# Calculate the number of unique groups in the annotation
annotation_row <- data.frame(Lineage = sapply(rownames(heatmap_matrix), function(x) {
  # First split by "|"
  components <- strsplit(x, split = "\\|")[[1]]
  if (length(components) >= 4) {
    # Now split the 4th component by ":"
    subcomponents <- strsplit(components[4], split = ":")[[1]]
    # Return the 2nd element if available, else NA
    if (length(subcomponents) >= 2) {
      return(subcomponents[3])
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}))
rownames(annotation_row) <- rownames(heatmap_matrix)

# Count the number of unique groups
num_groups <- length(unique(annotation_row$Lineage[!is.na(annotation_row$Lineage)]))



# Add the number of groups to the heatmap title
pheatmap(heatmap_matrix[, subSampleConds],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         annotation_row = annotation_row,
         show_rownames = FALSE,
         color = cividis_palette,
         filename = "heatmap_Repeats.png",
         main = paste("Heatmap Repeat Elements (", num_groups, " groups)", sep = ""))


# Filter out rows where Lineage contains "SINE?", "DNA?", or "LTR?"
lineSineLTR <- str_detect(annotation_row$Lineage, "SINE|LINE|LTR|DNA") & !str_detect(annotation_row$Lineage, "LINE.|SINE.|DNA.|LTR.")
# Subset the annotation and heatmap matrix
annotation_rowNew <- annotation_row[lineSineLTR, , drop = FALSE]
heatMap_mat <- heatmap_matrix[lineSineLTR, subSampleConds]
table(annotation_rowNew$Lineage)

# Generate the heatmap
pheatmap(heatMap_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         annotation_row = annotation_rowNew,
         show_rownames = FALSE,
         color = cividis_palette,
         filename = "heatmap_Repeats_LINE_SINE_LTR_DNA.png",
         main = "Heatmap Repeat Elements")





resultsDTSample <- resultsDTSample %>% filter(LFClabel %in% c("up", "down")) %>% 
mutate(ClassFamily = str_split_i(RepID, "\\|", 4)) %>%
separate(ClassFamily, into = c("subfamily", "family", "class"), sep = ":")


df_filteredCondClass <- resultsDTSample %>% select(condition, RepID, class, family, log2FoldChange, padj, LFClabel) %>%
  arrange(-log2FoldChange) %>%
  group_by(condition, class, LFClabel) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(count = ifelse(LFClabel == "down", -count, count))
  # arrange(-count)



# df_filtered <- resultsDTSample %>%
#   filter(LFClabel %in% c("up", "down")) %>%
#   group_by(condition, LFClabel) %>%
#   summarise(count = n(), .groups = "drop") %>%
#   mutate(count = ifelse(LFClabel == "down", -count, count))  # Make "down" counts negative


# df_ordered <- df_filtered %>%
#   filter(LFClabel == "up") %>%
#   arrange(desc(count)) %>%
#   pull(condition) %>%
#   unique()

# df_filtered$condition <- factor(df_filtered$condition, levels = subSampleConds)


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
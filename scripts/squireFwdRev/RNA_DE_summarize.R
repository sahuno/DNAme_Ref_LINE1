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
library(fst)
library(UpSetR)


options("width"=200)
set_fst_threads <- 4
threads_fst(set_fst_threads)


library(optparse)

option_list <- list(
    make_option(c("-l", "--lfcCutoff"), type="numeric", default=0.5, help="Cutoff for adjusted p-value [default %default]"),
    make_option(c("-d", "--paths_rdata"), type="character", default="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/FullLengthActive",help="Path to Rdata folder"))

parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

# After parsing, split comma-separated values into vectors
paths_rdata <- opt$paths_rdata

## Therapy list
combinedTherapy <- c("SETDB1i-CKi_SETDB1i" ,"SETDB1i-CKi_CKi","QSTAT-CKi_QSTAT" ,"QSTAT-CKi_CKi","SETDB1i-CKi_QSTAT-CKi")
monTherapy <- c("AZA_DMSO", "CKi_DMSO","QSTAT_DMSO", "QSTAT-CKi_DMSO", "SETDB1i_DMSO", "SETDB1i-CKi_DMSO")


# setdiff(c(monTherapy, combinedTherapy), colnames(log2FoldChangeDF))
# monTherapy <- c("DMSO_AZA","DMSO_CKi", "DMSO_QSTAT", "DMSO_SETDB1i", "QSTAT_SETDB1i")    
# combinedTherapy <- c("QSTAT_QSTAT-CKi","CKi_QSTAT-CKi" ,"CKi_SETDB1i-CKi","DMSO_SETDB1i-CKi","DMSO_QSTAT-CKi", "QSTAT-CKi_SETDB1i-CKi", "SETDB1i_SETDB1i-CKi"  )                  

### load RNAseq data/ full length active L1s
filesRDATA <- list.files(opt$paths_rdata, pattern = "ActiveL1ExpressionResults_", recursive = TRUE, full.names = TRUE) 

loadedRData <- lapply(filesRDATA, function(x) {
    e <- new.env()
    load(x, envir = e)
    as.list(e)
  })


  ##rename to list 
baseNames <- gsub(".*condition_", "", basename(filesRDATA))
names(loadedRData) <- gsub(".RData", "", baseNames)

# names(loadedRData[[1]])
# loadedRData[[1]][["results2SaveActiveL1"]]$resDF

resultsDT <- lapply(names(loadedRData), function(x){
  res <- loadedRData[[x]][["results2SaveActiveL1"]]$resDF
  res$L1_ID <- row.names(res)
  return(res)
})
names(resultsDT) <- names(loadedRData)

resultsDTSample <- rbindlist(resultsDT, idcol = "sample")

resultsDTSample <- resultsDTSample %>% mutate(label = case_when(log2FoldChange > 0 & padj <= opt$lfcCutoff ~ "up",
                                     log2FoldChange < 0 & padj <= opt$lfcCutoff ~ "down",
                                     TRUE ~  "no_change")) 


# loadedRData[[3]][["results2SaveActiveL1"]][["ddsFiltered"]]
# resultsNames(loadedRData[[1]][["results2SaveActiveL1"]][["ddsFiltered"]])

#make matrix for heatmap
df_filtered <- resultsDTSample %>%
  filter(label %in% c("up", "down")) %>%
  group_by(sample, label) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(count = ifelse(label == "down", -count, count))  # Make "down" counts negative

# Plot with ggplot2
plt<-ggplot(df_filtered, aes(x = sample, y = count, fill = label)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "Up and Down Counts per Sample", x = "Sample", y = "Count", fill = "Label") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plt, file = paste0(paths_rdata, "up_down_counts_per_sample.png"), width = 10, height = 6)


resultsDTSample %>% filter(label %in% c("up", "down")) %>% group_by(sample) %>%
arrange(-pvalue, log2FoldChange)


##Heatmap
log2FoldChangeDF <- resultsDTSample %>% filter(label %in% c("up", "down")) %>% select(!c("lfcSE","stat","pvalue","padj","label", "baseMean")) %>% pivot_wider(names_from = sample, values_from = log2FoldChange) #%>% 
log2FoldChangeDF <- log2FoldChangeDF %>% column_to_rownames("L1_ID") 
log2FoldChangeDF[is.na(log2FoldChangeDF)] = 0

pdf(file = paste0("clusteringlog2FoldChangeConditions.pdf"), width=12, height=9)
pheatmap::pheatmap(log2FoldChangeDF, main = paste0("abs(LFC) > 1 & padj < 0.05"), scale = "row", annotation_names_row=FALSE, show_rownames = TRUE, na_col = "black")
dev.off()



##########################
# Convert to binary matrix for UpSet plot
df_binary <- resultsDTSample %>% filter(label %in% c("up", "down")) %>% select(!c("lfcSE","stat","pvalue","padj","label", "baseMean", "log2FoldChange")) %>% 
  mutate(value = 1) %>%
  pivot_wider(names_from = sample, values_from = value, values_fill = 0)

# rowSums(df_binary[, -1])  # Check the row sums
# df_binary
# Plot UpSet
# pltUset <- upset(df_binary, 
#       sets = colnames(df_binary)[-1],  # Use sample names as sets
#       nsets = length(unique(df$sample)), 
#       nintersects = 10,  # Adjust as needed
#       keep.order = TRUE)

# Ensure L1_ID is not used as a set
df_matrix <- df_binary %>%
  column_to_rownames("L1_ID")  # Convert L1_ID to row names

# Convert to a standard dataframe (UpSetR does not like tibbles)
df_matrix <- as.data.frame(df_matrix)

# Ensure all columns are numeric (0 or 1)
# df_matrix[] <- lapply(df_matrix, function(x) as.numeric(as.character(x)))

# Generate the UpSet plot
# pltUset <- upset(df_matrix, 
#       sets = colnames(df_matrix), 
#       nsets = length(colnames(df_matrix)), 
#       nintersects = 10, 
#       keep.order = TRUE)

# Open a PNG graphics device
png("up_down_counts_per_sample.png", width = 13, height = 6, units = "in", res = 300)
# Generate and display the UpSet plot
upset(df_matrix, 
      sets = colnames(df_matrix), 
      nsets = length(colnames(df_matrix)), 
      nintersects = 40, 
      keep.order = TRUE)
dev.off()

# Save the plot
# ggsave(filename = "up_down_counts_per_sample.png", plot = pltUset, width = 10, height = 6, dpi = 300)

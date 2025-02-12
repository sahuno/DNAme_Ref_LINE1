library(data.table)
library(tidyverse)
library(optparse)

# import polars as pl
# import argpase as arg

library("optparse")
option_list <- list(make_option(c("-f", "--files"), type="character", help="files to process"))

parse_args(OptionParser(option_list=option_list))

# 
# mamba activate R
# list_paths_overlaps <- c('/data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/results/DNAme_overlaps/D-0-1_5000_4000/D-0-1_5000_4000_DNAme_RepEOverlaps.bed', '/data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/results/DNAme_overlaps/D-S-1_5000_4000/D-S-1_5000_4000_DNAme_RepEOverlaps.bed')
# meta_samples <- fread("/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq//metadata_triplicates_recoded.csv")

# list_paths_overlaps <- list.files("/data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/results", recursive = TRUE, full.names = TRUE, pattern="_Overlaps.bed")
# l /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/results/DNAme_overlaps/D-C-1_4000/overlaps/*_minCov15_CpGIslands.bed

# list_paths_overlaps <- list.files("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/results/DNAme_overlaps/", recursive = TRUE, full.names = TRUE, pattern="_DNAme_mmflil1_8438_Overlaps_minCov10_CpGIslands.bed")
list_paths_overlaps <- list.files("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/results/DNAme_overlaps/", recursive = TRUE, full.names = TRUE, pattern="_DNAme_mmflil1_8438_Overlaps_minCov10.bed")
PlotTag = "minCov10_noCGI"

# PlotTag = "minCov10_CpGIslands"

#remove unmerged sampels
#if you dont want don;t process
# removeSamplesOut <- c("D-0-1_4000_5mCpG_5hmCpG_DNAme_mmflil1_8438_", "D-0-1_5000_5mCpG_5hmCpG_DNAme_mmflil1_8438_","D-0-2_4000_5mCpG_5hmCpG_DNAme_mmflil1_8438_", "D-0-2_5000_5mCpG_5hmCpG_DNAme_mmflil1_8438_", "D-Q-1_5000_5mCpG_5hmCpG_DNAme_mmflil1_8438_", "D-Q-1_4000_5mCpG_5hmCpG_DNAme_mmflil1_8438_", "D-S-1_5000_5mCpG_5hmCpG_DNAme_mmflil1_8438_", "D-S-1_4000_5mCpG_5hmCpG_DNAme_mmflil1_8438_")
# list_paths_overlaps <- list_paths_overlaps[!grepl(paste(removeSamplesOut, collapse = "|"),list_paths_overlaps)]



DNAmeOverlaps <- lapply(list_paths_overlaps, function(x){fread(x)})
names(DNAmeOverlaps) <- basename(list_paths_overlaps)


##convert list to df for ploting 
DNAmeOverlaps_repeats <- rbindlist(DNAmeOverlaps, idcol = "samples")
names(DNAmeOverlaps_repeats) <- c("samples", "chrom", "start", "end", "RepeatID", "score", "strand", "start2", "end2", "color", "min", "max", "mean", "median", "count")
DNAmeOverlaps_repeats$samples <- gsub("_[0-9].*","", DNAmeOverlaps_repeats$samples)

#convert `.` to NA
DNAmeOverlaps_repeats <- DNAmeOverlaps_repeats[, lapply(.SD, function(x) {
  if (is.character(x)) {
    x[x == "."] <- NA
  }
  return(x)
}), .SDcols = names(DNAmeOverlaps_repeats)]


cols2Convert <- c("min", "max", "mean", "median", "count")
DNAmeOverlaps_repeats[, (cols2Convert):=lapply(.SD, as.numeric), .SDcols=cols2Convert]

mouseTriEpi_metadata <- read_csv("/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_recoded.csv")
mouseTriEpi_metadata$samples <- gsub("R-","D-", mouseTriEpi_metadata$new_samples_name)
mouseTriEpi_metadata$condition <- gsub("AZA","5-AZA", mouseTriEpi_metadata$condition)


DNAmeOverlaps_repeats2 <- DNAmeOverlaps_repeats %>% left_join(mouseTriEpi_metadata)

DNAmeOverlaps_repeats2_standChrom <- DNAmeOverlaps_repeats2 %>% filter(!str_detect(chrom, "chrM|chrY|chrX")) 

DNAmeOverlaps_repeats2_standChrom %>% group_by(samples) %>% summarise(n())
# table(DNAmeOverlaps_repeats2_standChrom$samples, DNAmeOverlaps_repeats2_standChrom$chrom)
#at least 5 valid cpgs
DNAmeOverlaps_repeats2_standChromFiltered <- DNAmeOverlaps_repeats2_standChrom %>% filter(count > 5)

DNAmeOverlaps_repeats2_standChromFiltered %>% group_by(samples) %>% summarise(n())

# DNAmeOverlaps_repeats[V1 == "chr",]
boxPlot_DNAmeMean <- ggplot(data= DNAmeOverlaps_repeats2_standChromFiltered, aes(samples, mean)) + geom_boxplot() + 
# geom_jitter() + 
labs(y = "Mean DNAme" , x = "samples", title = "full length L1 - mmflil1_8438") + 
facet_wrap(~condition, scales = "free_x") + 
theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(boxPlot_DNAmeMean, filename = paste0("boxPlot_DNAmeMean_", PlotTag,".png"))
# ncol(DNAmeOverlaps_repeats) - 3


histogram_DNAmeMean <- ggplot(data= DNAmeOverlaps_repeats2_standChromFiltered, aes(mean)) + geom_histogram() + 
labs(y = "Mean DNAme" , x = "samples", title = "full length L1 - mmflil1_8438") + 
facet_wrap(~condition + samples, scales = "free") + 
theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(histogram_DNAmeMean, filename = paste0("histogram_DNAmeMean_", PlotTag,".png"))
# ncol(DNAmeOverlaps_repeats) - 3

boxPlot_DNAmeMedian <- ggplot(data= DNAmeOverlaps_repeats2_standChromFiltered, aes(samples, median)) + geom_boxplot() + 
labs(y = "DNAme Median" , x = "samples", title = "full length L1 - mmflil1_8438") + 
facet_wrap(~condition, scales = "free_x") + 
theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(boxPlot_DNAmeMedian, filename = paste0("boxPlot_DNAmeMedian_",PlotTag , ".png"))


boxPlot_ValidCpGs <- ggplot(data= DNAmeOverlaps_repeats2_standChromFiltered, aes(samples, count)) + geom_boxplot() + 
labs(y = "counts ValidCpGs" , x = "samples", title = "full length L1 - mmflil1_8438") + 
facet_wrap(~condition, scales = "free_x") + 
theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(boxPlot_ValidCpGs, filename = paste0("boxPlot_ValidCpGs_", PlotTag, ".png"))
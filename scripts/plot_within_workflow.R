library(data.table)
library(tidyverse)
library(optparse)
library(ggrepel)
# import polars as pl
# import argpase as arg
message("begining ploting")
library("optparse")
option_list <- list(make_option(c("-f", "--files"), type="character", default = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/results/DNAme_overlaps//D-0-1_5000_4000/overlaps/D-0-1_5000_4000_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov10.bed /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/results/DNAme_overlaps//D-0-2_5000_4000/overlaps/D-0-2_5000_4000_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov10.bed /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/results/DNAme_overlaps//D-0-3_4000/overlaps/D-0-3_4000_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov10.bed /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/results/DNAme_overlaps//D-A-1_4000/overlaps/D-A-1_4000_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov10.bed", help="files to process"),
make_option(c("-m", "--metadata"), type="character", default= "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_DNArecoded.csv"))

args <- parse_args(OptionParser(option_list=option_list))

list_paths_overlaps <- str_split(args$files, " ")[[1]]

mouseTriEpi_metadata <- read_csv(args$metadata)

# list_paths_overlaps <- list.files("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/results/DNAme_overlaps/", recursive = TRUE, full.names = TRUE, pattern="_DNAme_mmflil1_8438_Overlaps_minCov10_CpGIslands.bed")

message("reading bed files")
DNAmeOverlaps <- lapply(list_paths_overlaps, function(x){fread(x)})

message("done reading bed files")
names(DNAmeOverlaps) <- basename(list_paths_overlaps)

#get plot tag
PlotTag <- gsub(".bed","",str_extract(names(DNAmeOverlaps), "5mCpG_[^_]+.*"))[1]

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


#join line 1 and metadata
DNAmeOverlaps_repeats2 <- DNAmeOverlaps_repeats %>% left_join(mouseTriEpi_metadata)

DNAmeOverlaps_repeats2_standChrom <- DNAmeOverlaps_repeats2 %>% filter(!str_detect(chrom, "chrM|chrY|chrX")) 

DNAmeOverlaps_repeats2_standChrom %>% group_by(samples) %>% summarise(n())
# table(DNAmeOverlaps_repeats2_standChrom$samples, DNAmeOverlaps_repeats2_standChrom$chrom)
#at least 5 valid cpgs
DNAmeOverlaps_repeats2_standChromFiltered <- DNAmeOverlaps_repeats2_standChrom %>% filter(count > 5)

# DNAmeOverlaps_repeats2_standChromFiltered %>% group_by(samples) %>% summarise(n())

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


# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/plot_within_workflow.R

##### plot individual line1 ### 
# split(DNAmeOverlaps_repeats2_standChromFiltered, by = RepeatID)
# l1 <- unique(DNAmeOverlaps_repeats2_standChromFiltered$RepeatID)

# RepeatIDGroups <- DNAmeOverlaps_repeats2_standChromFiltered %>% group_by(RepeatID) %>% group_split()
RepeatIDGroups <- split(DNAmeOverlaps_repeats2_standChromFiltered, as.factor(DNAmeOverlaps_repeats2_standChromFiltered$RepeatID))

# RepeatIDGroups[[1]]

#Que: is there any differnce between line 1 in any condition. that coudl be answered with dmr
#what is 
plotL1perCondition <- function(l1Name){
plotOut <- ggplot(RepeatIDGroups[[l1Name]], aes(x=condition, y=mean, fill=condition) ) + geom_boxplot() + 
geom_jitter() + 
geom_text_repel(aes(label=samples)) + labs(title = paste(l1Name, " minCov"), y = "mean 5mCpG_5mCpG") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
}

# plots_plotL1perCondition.pdf
message("creating dir results/boxplotsRepPerCond/ to save combined plots")
dir.create("results/boxplotsRepPerCond/")
pdf(paste0("results/boxplotsRepPerCond/combinedBoxPlot_", PlotTag, ".pdf"))
lapply(names(RepeatIDGroups), plotL1perCondition)
dev.off()

# lapply()

# gTest <- ggplot(RepeatIDGroups[[3]], aes(x=condition, y=mean, fill=condition) ) + geom_boxplot() + 
# geom_jitter() + 
# geom_text_repel(aes(label=samples)) + labs(title = "", y = "mean 5mCpG_5mCpG")

# ggsave(gTest, filename="testPerGenePlot.png")

# lapply(l1, function(x){
#   ggplot(data=DNAmeOverlaps_repeats2_standChromFiltered, aes())
# })
# if (str_detect(PlotTag, "_CpGIslands") == TRUE){
  
# }
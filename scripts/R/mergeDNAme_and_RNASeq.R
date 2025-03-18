#
# merge DNA methylation rates with RNAseq counts

library(data.table)
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(viridis)
library(magrittr)
library(pheatmap)
library(ggrepel)
library(ggfortify)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggnewscale)
library(cowplot)
library(RColorBrewer)
library(NbClust)
library(GenomicRanges)
library(gUtils)
library(janitor)
library(fst)


options("width"=200)

## To  run
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/locusSpecificRNA_DNAmeRate.R

#### program options/settings
ref_variable = "DMSO"
MIN_ReadsCounts = 10
smallestGroupSize <- 3
set_fst_threads <- 4
path_DNAmeRates <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/mmflil1_entireLength/resultsOverlapsStats_5mCpG_5hmCpG_sortedBed_minCov10.fst"
squireCounts_path <- "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/CT/squire_te_fwd.tsv"
counts_annot_path <- '/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/CT/counts_annot.tsv'
metadata_path <- "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_recoded.csv"
l1_path <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed"

threads_fst(set_fst_threads)


#read rates
DNAme_rates_dt <- read_fst(path_DNAmeRates, as.data.table = TRUE)
# DNAme_rates_dt[variable=="Freq_5mCpG.N" & value >= 3.0,]
DNAme_rates_geom_mean_dt_ <- DNAme_rates_dt[variable=="Freq_5mCpG.geom_mean",]

DNAme_rates_dt_statsWider <- dcast(DNAme_rates_dt, ... ~ variable, value.var = "value")[`Freq_5mCpG.N` >= 3,]
DNAme_rates_dt_statsWiderGmean <- DNAme_rates_dt_statsWider[,!c("Freq_5mCpG.N", "Freq_5mCpG.mean", "Freq_5mCpG.median", "Freq_5mCpG.entropy_in_log2")]
# summary($)

DNAme_rates_wide_dt <- dcast(DNAme_rates_dt_statsWiderGmean[,!c("id")], ... ~ samples, value.var = "Freq_5mCpG.geom_mean")


### load RNAseq data/ full length active L1s
paths_rdata <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs"
filesRDATA <- list.files(paths_rdata, pattern = "ActiveL1ExpressionResults_", recursive = TRUE, full.names = TRUE) 

loadedRData <- lapply(filesRDATA, function(x) {
    e <- new.env()
    load(x, envir = e)
    as.list(e)
  })


##rename to list 
baseNames <- gsub(".*condition_", "", basename(filesRDATA))
names(loadedRData) <- gsub(".RData", "", baseNames)
# load(filesRDATA[[1]])

#join rna & DNA
# loadedRData[["DMSO_CKi"]][["results2SaveActiveL1"]]$results2SaveActiveL1$cpmDF

cpmDF_list <- lapply(loadedRData, function(x) pivot_longer(x$results2SaveActiveL1$cpmDF, !genomicElementID, names_to = "RNAsamples", values_to = "cpm"))


cpmDF <- plyr::ldply(cpmDF_list, data.frame)
cpmDF <- mutate(cpmDF, sampleID = gsub("R.", "sN.", RNAsamples))
unique(cpmDF$sampleID)

DNAme_rates_dt_statsWiderGmean[, sampleID:=gsub("D-", "sN.", samples)][,sampleID:=gsub("-", ".", sampleID)]
unique(DNAme_rates_dt_statsWiderGmean$sampleID)

#merge
cpmDF_joined <- left_join(cpmDF, DNAme_rates_dt_statsWiderGmean, by = c("sampleID","genomicElementID"="RepeatID"))

rnaDNA_DMSO_AZA <- cpmDF_joined %>% filter(`.id` == "DMSO_AZA") 

# plot_DNARNA_L1 <- ggplot(rnaDNA_DMSO_AZA, aes(x = Freq_5mCpG.geom_mean, y = cpm)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~sampleID, scale = "free") + theme_minimal()
plot_DNARNA_L1 <- ggplot(rnaDNA_DMSO_AZA, aes(x = Freq_5mCpG.geom_mean, y = cpm)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~sampleID, scale = "free") + theme_minimal()
ggsave(plot_DNARNA_L1, filename = "DNA_and_RNA_L1_geomMean.png", width = 13, height = 10, units = "cm")
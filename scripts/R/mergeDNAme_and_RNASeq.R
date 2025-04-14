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
# library(pheatmap)
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


options("width"=200)
set_fst_threads <- 4
threads_fst(set_fst_threads)


library(optparse)

option_list <- list(
    make_option(c("-r", "--pathsDNAmeRates"), type="character", default=paste(
      "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAme_L1Base_SummaryStats/nonCpGIslands/full_Length_L1DNAme_Stats_100bpPromoter_minCov10.fst",
      "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAme_L1Base_SummaryStats/nonCpGIslands/full_Length_L1DNAme_Stats_400_600bpPromoter_minCov10.fst",
      "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAme_L1Base_SummaryStats/nonCpGIslands/full_Length_L1DNAme_Stats_wholeLengthPromoter_minCov10.fst",
      sep=","),
      help="Comma-separated list of paths for DNA methylation rates [default paths]"),

    make_option(c("-i", "--pathsDNAmeRates_id"), type="character", default="100bp5UTR,400bp600bp5UTR,entireLength",
      help="Comma-separated list of identifiers for DNA methylation rates paths [default identifiers]"),

    make_option(c("-d", "--paths_rdata"), type="character",
      default="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DE_repeats/ActiveL1/data",
      help="Path to Rdata folder")
)

parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

# After parsing, split comma-separated values into vectors
pathsDNAmeRates <- strsplit(opt$pathsDNAmeRates, ",")[[1]]
pathsDNAmeRates_id <- strsplit(opt$pathsDNAmeRates_id, ",")[[1]]
paths_rdata <- opt$paths_rdata

# Check results
print(pathsDNAmeRates)
print(pathsDNAmeRates_id)
print(paths_rdata)

## "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/DNAme_bed/version2/results/BedGraphsBigWigs/nonCpGIslands/D-0-1_5000_4000/D-0-1_5000_4000_CpGs_sortedBed_minCov15.bw"

## To  run
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/locusSpecificRNA_DNAmeRate.R
# pathsDNAmeRates <- c("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/mmflil1_100bp5UTR/resultsOverlapsStats_5mCpG_5hmCpG_sortedBed_minCov10.fst",
# "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/mmflil1_400bp600bp5UTR/resultsOverlapsStats_5mCpG_5hmCpG_sortedBed_minCov10.fst",
# "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/mmflil1_entireLength/resultsOverlapsStats_5mCpG_5hmCpG_sortedBed_minCov10.fst")

# pathsDNAmeRates <- c("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAme_L1Base_SummaryStats/nonCpGIslands/full_Length_L1DNAme_Stats_100bpPromoter_minCov10.fst",
# "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAme_L1Base_SummaryStats/nonCpGIslands/full_Length_L1DNAme_Stats_400_600bpPromoter_minCov10.fst",
# "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAme_L1Base_SummaryStats/nonCpGIslands/full_Length_L1DNAme_Stats_wholeLengthPromoter_minCov10.fst")

# pathsDNAmeRates_id <- c("100bp5UTR", "400bp600bp5UTR", "entireLength")
# paths_rdata <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DE_repeats/ActiveL1/data"


#### program options/settings
# path_DNAmeRates <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/mmflil1_100bp5UTR/resultsOverlapsStats_5mCpG_5hmCpG_sortedBed_minCov10.fst"
# path_DNAmeRates <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/mmflil1_entireLength/resultsOverlapsStats_5mCpG_5hmCpG_sortedBed_minCov10.fst"

# squireCounts_path <- "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/CT/squire_te_fwd.tsv"
# counts_annot_path <- '/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/CT/counts_annot.tsv'
# metadata_path <- "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_recoded.csv"
# l1_path <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed"
# paths_rdata <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs"


# #read rates
# DNAme_rates_dt <- read_fst(path_DNAmeRates, as.data.table = TRUE)
# # DNAme_rates_dt[variable=="Freq_5mCpG.N" & value >= 3.0,]
# DNAme_rates_geom_mean_dt_ <- DNAme_rates_dt[variable=="Freq_5mCpG.geom_mean",]

# DNAme_rates_dt_statsWider <- dcast(DNAme_rates_dt, ... ~ variable, value.var = "value")[`Freq_5mCpG.N` >= 3,]
# DNAme_rates_dt_statsWiderGmean <- DNAme_rates_dt_statsWider[,!c("Freq_5mCpG.N", "Freq_5mCpG.mean", "Freq_5mCpG.median", "Freq_5mCpG.entropy_in_log2")]
# # summary($)

# DNAme_rates_wide_dt <- dcast(DNAme_rates_dt_statsWiderGmean[,!c("id")], ... ~ samples, value.var = "Freq_5mCpG.geom_mean")


# ### load RNAseq data/ full length active L1s
# filesRDATA <- list.files(paths_rdata, pattern = "ActiveL1ExpressionResults_", recursive = TRUE, full.names = TRUE) 

# loadedRData <- lapply(filesRDATA, function(x) {
#     e <- new.env()
#     load(x, envir = e)
#     as.list(e)
#   })


# ##rename to list 
# baseNames <- gsub(".*condition_", "", basename(filesRDATA))
# names(loadedRData) <- gsub(".RData", "", baseNames)
# # load(filesRDATA[[1]])

# #join rna & DNA
# # loadedRData[["DMSO_CKi"]][["results2SaveActiveL1"]]$results2SaveActiveL1$cpmDF

# cpmDF_list <- lapply(loadedRData, function(x) pivot_longer(x$results2SaveActiveL1$cpmDF, !genomicElementID, names_to = "RNAsamples", values_to = "cpm"))


# cpmDF <- plyr::ldply(cpmDF_list, data.frame)
# cpmDF <- mutate(cpmDF, sampleID = gsub("R.", "sN.", RNAsamples))
# unique(cpmDF$sampleID)
# unique(cpmDF$`.id`)
# AzaDmsoDF <- cpmDF %>% dplyr::filter(`.id`=="DMSO_AZA")
# write_tsv(AzaDmsoDF, file = "FullLengthActiveL1_CPM.tsv")


# DNAme_rates_dt_statsWiderGmean[, sampleID:=gsub("D-", "sN.", samples)][,sampleID:=gsub("-", ".", sampleID)]
# unique(DNAme_rates_dt_statsWiderGmean$sampleID)

# #merge
# cpmDF_joined <- left_join(cpmDF, DNAme_rates_dt_statsWiderGmean, by = c("sampleID","genomicElementID"="RepeatID"))

# rnaDNA_DMSO_AZA <- cpmDF_joined %>% filter(`.id` == "DMSO_AZA") 

# # plot_DNARNA_L1 <- ggplot(rnaDNA_DMSO_AZA, aes(x = Freq_5mCpG.geom_mean, y = cpm)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~sampleID, scale = "free") + theme_minimal()
# plot_DNARNA_L1 <- ggplot(rnaDNA_DMSO_AZA, aes(x = Freq_5mCpG.geom_mean, y = cpm)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~sampleID, scale = "free") + theme_minimal()
# ggsave(plot_DNARNA_L1, filename = "DNA_and_RNA_L1_geomMean.png", width = 13, height = 10, units = "cm")





##########################################################################################################
##############################################################################################################################
# runnig in a function
##############################################################################################################################

# pathsDNAme <- list.files(dir_DNAme, recursive=TRUE, pattern = "")

runDNAmeRNA <-function(RNA_path, DNAme_path, plt_title = ""){

#read rates
DNAme_rates_dt <- read_fst(DNAme_path, as.data.table = TRUE)
# DNAme_rates_dt[variable=="Freq_5mCpG.N" & value >= 3.0,]
DNAme_rates_geom_mean_dt_ <- DNAme_rates_dt[variable=="Freq_5mCpG.geom_mean",]

DNAme_rates_dt_statsWider <- dcast(DNAme_rates_dt, ... ~ variable, value.var = "value")[`Freq_5mCpG.N` >= 3,]
DNAme_rates_dt_statsWiderGmean <- DNAme_rates_dt_statsWider[,!c("Freq_5mCpG.N", "Freq_5mCpG.mean", "Freq_5mCpG.median", "Freq_5mCpG.entropy_in_log2")]
# summary($)

DNAme_rates_wide_dt <- dcast(DNAme_rates_dt_statsWiderGmean[,!c("id")], ... ~ samples, value.var = "Freq_5mCpG.geom_mean")


### load RNAseq data/ full length active L1s
filesRDATA <- list.files(RNA_path, pattern = "ActiveL1ExpressionResults_", recursive = TRUE, full.names = TRUE) 

loadedRData <- lapply(filesRDATA, function(x) {
    e <- new.env()
    load(x, envir = e)
    as.list(e)
  })


##rename to list 
baseNames <- gsub(".*condition_", "", basename(filesRDATA))
names(loadedRData) <- gsub(".RData", "", baseNames)

cpmDF_list <- lapply(loadedRData, function(x) pivot_longer(x$results2SaveActiveL1$cpmDF, !genomicElementID, names_to = "RNAsamples", values_to = "cpm"))


cpmDF <- plyr::ldply(cpmDF_list, data.frame)
cpmDF <- mutate(cpmDF, sampleID = gsub("R.", "sN.", RNAsamples))
unique(cpmDF$sampleID)

datAzaDMSO <- cpmDF %>% dplyr::filter(`.id`=="DMSO_AZA")
unique(datAzaDMSO$genomicElementID)
write_tsv(datAzaDMSO, file = "FullLengthActiveL1_CPM.tsv")
# unique(cpmDF[`.id`=="DMSO_AZA", c(".id","genomicElementID")])

DNAme_rates_dt_statsWiderGmean[, sampleID:=gsub("D-", "sN.", samples)][,sampleID:=gsub("-", ".", sampleID)]
# unique(DNAme_rates_dt_statsWiderGmean$sampleID)

#merge
cpmDF_joined <- left_join(cpmDF, DNAme_rates_dt_statsWiderGmean, by = c("sampleID","genomicElementID"="RepeatID"))

rnaDNA_DMSO_AZA <- cpmDF_joined %>% filter(`.id` == "DMSO_AZA") 

# plot_DNARNA_L1 <- ggplot(rnaDNA_DMSO_AZA, aes(x = Freq_5mCpG.geom_mean, y = cpm)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~sampleID, scale = "free") + theme_minimal()
plot_DNARNA_L1 <- ggplot(rnaDNA_DMSO_AZA, aes(x = Freq_5mCpG.geom_mean, y = cpm)) + 
geom_point() + 
geom_smooth(method = "lm") + 
labs(title = plt_title) +
facet_wrap(~sampleID, scale = "free_y") + theme_minimal()

# ggsave(plot_DNARNA_L1, filename = paste0("DNA_and_RNA_L1_geomMean_",plt_title,".svg"), width = 13, height = 10, units = "cm")
tryCatch({
    ggsave(plot_DNARNA_L1, 
           filename = paste0("DNA_and_RNA_L1_geomMean_", plt_title, ".svg"), 
           width = 13, height = 10, units = "cm")
}, error = function(e) {
    message("ggsave failed with error: ", e$message)
    # Optionally log the error or continue without stopping
})

}


# runDNAmeRNA <-function(RNA_path, DNAme_path, plt_title = ""){
# x=3

lapply(seq_along(pathsDNAmeRates_id), function(x) {runDNAmeRNA(RNA_path=paths_rdata, DNAme_path=pathsDNAmeRates[x], plt_title = pathsDNAmeRates_id[x])})

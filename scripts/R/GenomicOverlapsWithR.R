#author: find overlaps with R and compute other stats
#sanity checks
library("data.table")
library(GenomicRanges)
library(entropy)
library(tidyverse)
library("fst")

# source("/data1/greenbab/users/ahunos/methylONT/utils_onts_downstream.R")


# ssend r-env overlaps /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R

# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R

# source("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R")

#read bed files
paths_bed = list.files("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/prepareBedFiles/", pattern="*minCov10.bed$",recursive=TRUE, full.names=TRUE)
# paths_bed <- paths_bed[1:2]
# args(write_fst)

# fst_files <- multi_mpileup_fst(paths=paths_bed, colnames)
# fst::read_fst(paths_bed[[1]]) 

######### read data
message("reading samples bed files")
ls_data = lapply(paths_bed, fread, col.names = c("chrom", "start", "end", "name", "score", "strand", "tstart", "tend", "color", "coverage", "Freq_5mCpG", "mod", "canon", "other", "del", "fail", "diff", "nocall")) 
names(ls_data) <- str_extract(basename(paths_bed), "^[^_]+")
message("done reading files")
gr_data <- lapply(ls_data, makeGRangesFromDataFrame, keep.extra.columns = TRUE)

# paths_bed
PlotTag <- gsub(".bed","",str_extract(paths_bed, "5mCpG_[^_]+.*"))[1]


########## read the active l1 files from l1base
mm10_fli <- '/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed'

dt_mfli <- fread(mm10_fli, col.names = c("chrom", "start", "end", "RepeatID", "score", "strand", "start2", "end2", "color"))
dr_mfli <- makeGRangesFromDataFrame(dt_mfli, keep.extra.columns = TRUE)
### make id's out of line1
dr_mfli$id <- paste0(seqnames(dr_mfli), ":", paste0(ranges(dr_mfli)))

dt_mfliUnique <- unique(as.data.table(mcols(dr_mfli)[, c("RepeatID","id")]))

######################### define list of functions to use
#entropy calculations
entropy_cal <- function(x, method2Use = "ML", na.rm = TRUE, units2Use = "log2"){
  res <- entropy::entropy(x, method=method2Use, unit=units2Use) 
  as.numeric(res)
}

EntropyCS <- function(x, method2Use = "CS", na.rm = TRUE, units2Use = "log2"){
  res <- entropy::entropy(x, method=method2Use, unit=units2Use) 
  as.numeric(res)
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

category = c("RepeatID")
ColChoice = c("Freq_5mCpG")
my.summary = function(x) list(N = length(x), mean = mean(x, na.rm =TRUE), median = median(x, na.rm =TRUE), entropy_in_log2=entropy_cal(x), geom_mean=gm_mean(x))
####################################################################################################

####putting all in a fucntion
Overlaps <- function(x){
    # replace `gr_data[[1]]` with `x`
ovs1 <- gUtils::gr.findoverlaps(x, dr_mfli, scol = c("RepeatID", "id"), qcol="Freq_5mCpG", return.type = "data.table")

ovs1_stats <- ovs1[!is.na('Freq_5mCpG'), unlist(lapply(.SD, my.summary), recursive = FALSE), .SDcols = ColChoice, by = category]

ovs1_stats <- ovs1_stats[dt_mfliUnique, on = c("RepeatID"), nomatch = NULL]

ovs1_stats_melted <- melt(ovs1_stats, id.vars=c("id", "RepeatID"))
return(ovs1_stats_melted)
}

# computeStats <- function(x){
# ovs1_stats <- x[!is.na('Freq_5mCpG'), unlist(lapply(.SD, my.summary), recursive = FALSE), .SDcols = ColChoice, by = category]
# ovs1_stats <- ovs1_stats[dt_mfliUnique, on = c("RepeatID"), nomatch = NULL]
# ovs1_stats_melted <- melt(ovs1_stats, id.vars=c("id", "RepeatID") )
# return(ovs1_stats_melted)
# }

# gr_data[!is.na(`Freq_5mCpG`)
# ovs1 <- gUtils::gr.findoverlaps(gr_data[[1]], dr_mfli, scol = c("RepeatID", "id"), qcol="coverage", return.type = "data.table")

message("done computing stats")
resultsOverlapsStats <- lapply(gr_data, function(x){
  dt_ovlps <- Overlaps(x)
})
    # message("done overlaps, repeats & cpG sites")
# return(dt_ovlps)

message("bind results")
resultsOverlapsStats_dt <- rbindlist(resultsOverlapsStats, idcol = "samples")

# message("done computing stats")
# resultsOverlapsStats <- lapply(gr_data, function(x){
#     dt_ovlps_stats <- computeStats(dt_ovlps)
#     message("done overlaps, repeats & cpG sites")
#     return(dt_ovlps)
# })


##save a fst copy for fast io 
# ovs1_stats[,summarise(`Freq_5mCpG.N`)]
message("save .fst results to disk")
write_fst(resultsOverlapsStats_dt, path = paste0("resultsOverlapsStats_",PlotTag,".fst"), compress = 100)


message("plotting")
plotAllStats_perSample <- function(x){
distri_plots <- ggplot(resultsOverlapsStats[[x]][!grepl("chrX|chrY|chrM", id), ], aes(value)) + geom_histogram() + facet_wrap(~variable, scale = "free") + 
labs(title = paste0(x, ": 5mCpG methylation mm10 Active full length L1 (mmflil1_8438)"), x="values", y="counts") + theme_minimal()
}

#plot stats
list_stats_dt <- split(resultsOverlapsStats_dt, as.factor(resultsOverlapsStats_dt$variable))

plotStats_perSample <- function(x){
distri_plots <- ggplot(list_stats_dt[[x]][!grepl("chrX|chrY|chrM", id), ], aes(value)) + geom_histogram() + facet_wrap(~samples, scale = "free") + 
labs(title = paste0(x, ": 5mCpG methylation mm10 Active full length L1 (mmflil1_8438)"), x=x, y="counts") + theme_minimal()
}


message("saving plots")
# ovs1_stats_melted[[1]][!grepl("chrX|chrY|chrM", id), ]
# distri_plots <- ggplot(ovs1_stats_melted %>% dplyr::filter(!str_detect(id, "chrX|chrY|chrM")), aes(value)) + geom_histogram() + facet_wrap(~variable, scale = "free") + theme_minimal()

# ggsave(distri_plots, filename = "D01_distributions_methylation.png")

pdf(paste0("full_length_L1RNA_DNAme_DMSO_1samplePerPage",PlotTag,".pdf"), width=9, height=9)
lapply(names(resultsOverlapsStats), plotAllStats_perSample)
dev.off()

pdf(paste0("full_length_L1RNA_DNAme_DMSO_stats_perPage",PlotTag,".pdf"), width=9, height=7)
lapply(names(list_stats_dt), plotStats_perSample)
dev.off()

message("done running analysis")
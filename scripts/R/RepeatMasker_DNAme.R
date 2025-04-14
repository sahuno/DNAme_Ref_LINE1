#author: DNAme of repeat masker elements

#load libraries
library("data.table")
library(GenomicRanges)
library(entropy)
library(tidyverse)
library(fst)
library(optparse)
library(gUtils)
####READ


###get inputs from user
option_list <- list(
    make_option(c("-r", "--repMasker"), type="character", default="", 
        help="repeat Masker annnotations .bed [default %default]"),
        make_option(c("-b", "--DNAme"), type="character", 
        help="path to DNA methyllation bed files"),
        make_option(c("-k", "--key"), type="character", default="_minCov10.bed", 
        help="suffix or type of bed file used [default %default]"))

opt <- parse_args(OptionParser(option_list=option_list))



chrs <- paste0("chr", c(1:19))
# chrs <- paste0("chr", c(1:22, "X", "Y", "MT"))



bedFile <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAmeQC_Overlaps_L1Base/results/prepareBedFiles/nonCpGIslands/D-0-2_5000_4000/D-0-2_5000_4000_CpGs_sortedBed_minCov5.bed"

bedDNAmeDT <- fread(bedFile)
sortedRepMasker <- fread("/data1/greenbab/database/RepeatMaskerDB/repeatmasker_dot_org/mm10/rmsk405_20140131/rm_mm10_rm_sorted.bed")

setnames(bedDNAmeDT, old = c("V1", "V2", "V3", "V4",  "V5", "V6", "V11"), new = c("chrom", "start", "end", "name", "score", "strand", "fracMethyl") )
setnames(sortedRepMasker, old = c("V1", "V2", "V3", "V4",  "V5", "V6", "V7", "V8"), new = c("chrom", "start", "end", "name", "score", "strand", "family", "subFamily") )

#keep standard chrs only
bedDNAmeDT <- bedDNAmeDT[chrom %in% chrs, ]
# sum(duplicated(sortedRepMasker$name))

sortedRepMasker <- sortedRepMasker[chrom %in% chrs, ]
# table(sortedRepMasker$chrom)

sortedRepMasker[,ID := .I,]

#convert to grnages for overlaps
rmGr <- dt2gr(sortedRepMasker)
DNAmeGr <- dt2gr(bedDNAmeDT)

seqnames(DNAmeGr)
seqnames(rmGr)

                                   
#overlaps between dna methylation and 
ovs1 <- gUtils::gr.findoverlaps(DNAmeGr, rmGr, qcol=c("fracMethyl"), scol = c("name", "family","subFamily", "ID"))


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

category = c("ID")
ColChoice = c("fracMethyl")
computeSummary = function(x) list(N = length(x), mean = mean(x, na.rm =TRUE), median = median(x, na.rm =TRUE), entropy_in_log2=entropy_cal(x), geom_mean=gm_mean(x))
####################################################################################################

#convert back to data.table for computations
# ovs1 <- sortSeqlevels(ovs1)
# ovs1 <- sort(ovs1)

ovlapsDT <- gr2dt(ovs1)

#compute stats methylation
ovlapsStatsDT <- ovlapsDT[!is.na('fracMethyl'), unlist(lapply(.SD, computeSummary), recursive = FALSE), .SDcols = ColChoice, by = category]

# get annotations of REPEATS with DNAme Summary stats
colsInterest <-  c("seqnames", "start", "end", "strand", "width", "name", "family","subFamily","ID")
ovlapsStatsDTFull <- ovlapsStatsDT[ovlapsDT[ , ..colsInterest], on = "ID" ]

# ovlapsStatsDTFull[family == "LINE", .N, by=subFamily]
# table(ovlapsStatsDTFull$seqnames)
# table(ovlapsStatsDTFull$family)
# table(ovlapsStatsDTFull$family)
# unique(ovlapsStatsDTFull$family)

# c("SINE", "Simple_repeat", "Satellite", "scRNA", "LINE")
plt_methyl <- ggplot(data = ovlapsStatsDTFull[family %in% c("SINE", "Simple_repeat", "Satellite", "scRNA", "LINE")], aes(x=family, y=fracMethyl.mean)) + geom_boxplot()
ggsave(plt_methyl, filename = paste0("D02_repeats_methylation.pdf"))

# L1DNAme <- ovlapsStatsDTFull[subFamily %in% c("L1")][,.N, by=name]


plt_methylL1 <- ggplot(data = ovlapsStatsDTFull[subFamily == "L1"], aes(x=name, y=fracMethyl.mean)) + 
geom_boxplot() + 
labs(title = "Reference LINE-1 SubFamilies") + 
theme(axis.text.x = element_text(angle = 90))
ggsave(plt_methylL1, filename = paste0("D02_L1_methylation.pdf"), width = 18, height = 7)

write.fst(ovlapsStatsDTFull, path=paste0("DNAme_Stats_RepeatMasker_D02.fst"))

####################################
## message("done running analysis")
####################################
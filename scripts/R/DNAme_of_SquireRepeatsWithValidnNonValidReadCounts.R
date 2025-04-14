#author: DNAme of SquireElements

#load libraries
library("data.table")
library(GenomicRanges)
library(entropy)
library(tidyverse)
library(fst)
library(optparse)
library(gUtils)
library(ggpubr)

####Test Run

# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/DNAme_of_SquireRepeatsWithValidReadCounts.R \
# --RepeatRegions "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/SQuireRepeatsValidReadCounts_Aggregated_FWD_REV.fst" \
# --FileDNAme "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAmeQC_Overlaps_L1Base/results/prepareBedFiles/nonCpGIslands/D-0-2_5000_4000/D-0-2_5000_4000_CpGs_sortedBed_minCov5.bed" \
# --sampleName "D-D-1" \
# --chrs NULL

###get inputs from user
option_list <- list(
    make_option(c("-r", "--RepeatRegions"), type="character", default="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/SQuireRepeatsValidnNonValidReadCounts_Aggregated_FWD_REV.fst", 
        help="repeat annnotations .bed [default %default]"),
        make_option(c("-f", "--FileDNAme"), type="character",
        help="path to DNA methylation bed files"),
        make_option(c("-s", "--sampleName"), type="character",
        help="sampleName; check snakemake wildcards"),
        make_option(c("-c", "--chrs"), type="character",
        help="chromosomes to analyse"))

opt <- parse_args(OptionParser(option_list=option_list))
 
message("parsed all aguments")

if(is.null(opt$chrs)){
opt$chrs <- paste0("chr", c(1:19))
}

# chrs <- paste0("chr", c(1:22, "X", "Y", "MT"))



# opt$FileDNAme <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAmeQC_Overlaps_L1Base/results/prepareBedFiles/nonCpGIslands/D-0-2_5000_4000/D-0-2_5000_4000_CpGs_sortedBed_minCov5.bed"

bedDNAmeDT <- fread(opt$FileDNAme)
# RepDT <- fread(opt$RepeatRegions)
RepDT <- read_fst(path = opt$RepeatRegions, as.data.table = TRUE)


# RepDT <- fread(opt$RepeatRegions)

# RepDT <- fread("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/SQuireRepeatsValidReadCounts_Aggregated_FWD_REV.fst")

setnames(bedDNAmeDT, old = c("V1", "V2", "V3", "V4",  "V5", "V6", "V11"), new = c("chrom", "start", "end", "name", "score", "strand", "fracMethyl") )
# setnames(sortedRepMasker, old = c("V1", "V2", "V3", "V4",  "V5", "V6", "V7", "V8"), new = c("chrom", "start", "end", "name", "score", "strand", "family", "subFamily") )

#keep standard chrs only
bedDNAmeDT <- bedDNAmeDT[chrom %in% opt$chrs, ]
# sum(duplicated(sortedRepMasker$name))

RepDT <- RepDT[chr %in% opt$chrs, ]


# sum(duplicated(RepDT$TE_ID))
# RepDT[duplicated(RepDT$TE_ID)]
RepDT[TE_ID == "chr1|4773954|4774028|MIRc:MIR:SINE|192|-",]

RepDT[,ID := .I,]
# anyNA(RepDT)
# RepDT[,.N,by="chr" ]

RepGr <- makeGRangesFromDataFrame(as.data.frame(RepDT), keep.extra.columns = TRUE)
# table(RepDT$chr)

# RepDT[,chrom:=chr]
#convert to grnages for overlaps
# RepGr <- gUtils::dt2gr(RepDT)
DNAmeGr <- dt2gr(bedDNAmeDT)

DNAmeGr <- sortSeqlevels(DNAmeGr)
RepGr <- sortSeqlevels(RepGr)
seqnames(DNAmeGr)
seqnames(RepGr)

                                   
#overlaps between dna methylation and 
# RepGr_list <- split(RepGr, RepGr$TE_ID)
# seqinfo("mm10")
# genome(RepGr) <- "mm10"
# genome(DNAmeGr) <- "mm10"


colsInterested <- c("name", "number" , "consensus", "clade",  "class", "validReadCounts", "TE_ID", "ID")
ovs1 <- gUtils::gr.findoverlaps(DNAmeGr, RepGr, qcol=c("fracMethyl"), scol = colsInterested)


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
colsInterest <-  c("seqnames", "start", "end", "strand", "width", colsInterested)
# dim(ovlapsDT[ , ..colsInterest])
# dim(unique(ovlapsDT[ , ..colsInterested]))
# dim(unique(RepDT[ , ..colsInterested]))

# dim(unique(RepDT[TE_ID == "chr9|123364995|123365889|Lx5c:L1:LINE|149|+", c("TE_ID", "ID")]))
# RepDT[TE_ID == "chr9|123364995|123365889|Lx5c:L1:LINE|149|+", ..colsInterested]
# unique(RepDT[TE_ID == "chr9|123364995|123365889|Lx5c:L1:LINE|149|+" , ..colsInterested])

ovlapsStatsDTFull <- ovlapsStatsDT[unique(ovlapsDT[ , c("validReadCounts","TE_ID","ID")]), on = c("ID") ]
# ovlapsStatsDTFullNew <- copy(ovlapsStatsDTFull)[,c("chr", "start", "end", "consensusCladeClass","number","strand") := tstrsplit(`TE_ID`, "|", fixed = TRUE)] ##split the `TE_ID` column
ovlapsStatsDTFull[,c("chr", "start", "end", "name","number","strand") := tstrsplit(`TE_ID`, "|", fixed = TRUE)][, c("consensus","clade", "class"):=tstrsplit(name, ":", fixed = TRUE)] ##split the `TE_ID` column


# consensusCladeClass

# ovlapsStatsDTFull <- ovlapsStatsDT[ovlapsDT[ , ..colsInterest], on = "ID" ]
L1MD <- ovlapsStatsDTFull[validReadCounts == TRUE & clade == "L1" & grepl("L1Md",consensus)]

# ovlapsStatsDTFull[family == "LINE", .N, by=subFamily]
# table(ovlapsStatsDTFull$seqnames)
# table(ovlapsStatsDTFull$class)
# table(ovlapsStatsDTFull$clade)
# unique(ovlapsStatsDTFull$family)

fwrite(ovlapsStatsDTFull, file = paste0(opt$sampleName,"_repeatsDNAme.tsv"))
# write.fst(ovlapsStatsDTFull, path=paste0("DNAme_Stats_RepeatMasker_D02.fst"))



# c("SINE", "Simple_repeat", "Satellite", "scRNA", "LINE")
plt_methyl <- ggboxplot(ovlapsStatsDTFull[validReadCounts == TRUE & class %in% c("SINE", "scRNA", "LINE")], x = "class", y = "fracMethyl.mean", palette = "jco")+
  stat_compare_means(method = "anova", vjust=0.5) +      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")  +
labs(title = "Reference Repeats") + 
theme(axis.text.x = element_text()) + theme_pubclean()
ggsave(plt_methyl, filename = paste0(opt$sampleName, "_repeats_methylation.pdf"), width = 18, height = 7)


# L1DNAme <- ovlapsStatsDTFull[subFamily %in% c("L1")][,.N, by=name]

plt_methylL1 <- ggplot(data = ovlapsStatsDTFull[validReadCounts == TRUE & clade == "L1"], aes(x=consensus, y=fracMethyl.mean)) + 
geom_boxplot() + 
labs(title = "Reference LINE-1 SubFamilies") + 
theme(axis.text.x = element_text(angle = 90))
ggsave(plt_methylL1, filename = paste0(opt$sampleName,"_L1_methylation.pdf"), width = 18, height = 7)


# plt_methylL1mD <- ggplot(data = L1MD, aes(x=consensus, y=fracMethyl.mean)) + 
# geom_boxplot() + 
# labs(title = "Reference LINE-1 SubFamilies") + 
# theme(axis.text.x = element_text(angle = 90))
# ggsave(plt_methylL1mD, filename = paste0("D02_L1MD_methylation.pdf"), width = 18, height = 7)


# Visualize
plt_methylL1mD <- ggboxplot(L1MD, x = "consensus", y = "fracMethyl.mean", palette = "jco")+
  stat_compare_means(method = "anova", vjust=0.5) +      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")  +
labs(title = "Reference LINE-1 SubFamilies") + 
theme(axis.text.x = element_text()) + theme_pubclean()
ggsave(plt_methylL1mD, filename = paste0(opt$sampleName,"_L1MD_methylation_pval.pdf"), width = 18, height = 7)

####################################
message("\n\n done running analysis\n")
####################################
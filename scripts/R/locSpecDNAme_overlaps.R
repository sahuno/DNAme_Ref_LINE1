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
library(GGally)

###get inputs from user
option_list <- list(
    make_option(c("-r", "--RepeatRegions"), type="character", default="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LINE1_promoters/l1_300bp_5UTR_stranded.bed", 
        help="repeat annnotations .bed [default %default]"),
        make_option(c("-f", "--FileDNAme"), type="character",
        help="path to DNA methylation bed files"),
        make_option(c("-s", "--sampleName"), type="character",
        help="sampleName; check snakemake wildcards"),
        make_option(c("-c", "--chrs"), type="character",
        help="chromosomes to analyse"),
        make_option(c("-p", "--runPlots"), type="logical", default=TRUE,
        help="whether to run plotting sections [default %default]"),
        make_option(c("-o", "--outDir"), type="character", default=".", help = "output directory [default %default]"),
        make_option(c("-R", "--RNASeqCountsPath"), type="character", default=NULL, help = "path to RNAseq counts table [default %default]")
        )

opt <- parse_args(OptionParser(option_list=option_list))
 
message("parsed all aguments")

if(is.null(opt$chrs)){
opt$chrs <- paste0("chr", c(1:19))
}

# chrs <- paste0("chr", c(1:22, "X", "Y", "MT"))



# opt$FileDNAme <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAmeQC_Overlaps_L1Base/results/prepareBedFiles/nonCpGIslands/D-0-2_5000_4000/D-0-2_5000_4000_CpGs_sortedBed_minCov5.bed"

bedDNAmeDT <- fread(opt$FileDNAme)
# RepDT <- fread(opt$RepeatRegions)
RepDT <- fread(opt$RepeatRegions)
# RepDT <- read_fst(path = opt$RepeatRegions, as.data.table = TRUE)



# RepDT <- fread(opt$RepeatRegions)

# RepDT <- fread("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/SQuireRepeatsValidReadCounts_Aggregated_FWD_REV.fst")

setnames(bedDNAmeDT, old = c("V1", "V2", "V3", "V4",  "V5", "V6", "V11"), new = c("chrom", "start", "end", "name", "score", "strand", "fracMethyl") )
# setnames(sortedRepMasker, old = c("V1", "V2", "V3", "V4",  "V5", "V6", "V7", "V8"), new = c("chrom", "start", "end", "name", "score", "strand", "family", "subFamily") )


if(ncol(RepDT) > 10){
    setnames(RepDT, old = c("V1", "V2", "V3", "V4",  "V5", "V6", "V7", "V8", "V9", "V10", "V11"), 
    new = c("chrom", "start", "end", "name", "score", "strand","consensus", "clade", "class", "validReadCounts", "TE_ID") )
}

#keep standard chrs only
bedDNAmeDT <- bedDNAmeDT[chrom %in% opt$chrs, ]
# sum(duplicated(sortedRepMasker$name))

RepDT <- RepDT[chrom %in% opt$chrs, ]


# sum(duplicated(RepDT$TE_ID))
# RepDT[duplicated(RepDT$TE_ID)]
# RepDT[TE_ID == "chr1|4773954|4774028|MIRc:MIR:SINE|192|-",]
RepDT[,ID := .I,]
RepGr <- makeGRangesFromDataFrame(as.data.frame(RepDT), keep.extra.columns = TRUE)

# RepDT[,chrom:=chr]
#convert to grnages for overlaps
# RepGr <- gUtils::dt2gr(RepDT)
DNAmeGr <- dt2gr(bedDNAmeDT)

DNAmeGr <- sortSeqlevels(DNAmeGr)
RepGr <- sortSeqlevels(RepGr)
# seqnames(DNAmeGr)
# seqnames(RepGr)

                                   
#overlaps between dna methylation and 
# RepGr_list <- split(RepGr, RepGr$TE_ID)
# seqinfo("mm10")
# genome(RepGr) <- "mm10"
# genome(DNAmeGr) <- "mm10"


colsInterested <- c("name", "consensus", "clade",  "class", "validReadCounts", "TE_ID", "ID")
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
# table(L1MD$consensus)
L1Order <- c("L1Md_A", "L1Md_T", "L1Md_F", "L1Md_F2", "L1Md_F3", "L1Md_Gf")
L1MD$consensus <- factor(L1MD$consensus, levels = L1Order, ordered = TRUE)

# ovlapsStatsDTFull[family == "LINE", .N, by=subFamily]
# table(ovlapsStatsDTFull$seqnames)
# table(ovlapsStatsDTFull$class)
# table(ovlapsStatsDTFull$clade)
# unique(ovlapsStatsDTFull$family)
ovlapsStatsDTFull$filename <-  opt$sampleName
ovlapsStatsDTFull$sample <-  gsub("_.*","",ovlapsStatsDTFull$filename)
ovlapsStatsDTFull$RNAsampleID <- gsub("^D-","R-",ovlapsStatsDTFull$sample) #needed to merge with RNAseq data


if(!is.null(opt$RNASeqCountsPath)){
# RNASeqCountsPath <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/squireRepeats_AggregatedFWDnREV_masterTable.fst"
RNAdt <- read_fst(opt$RNASeqCountsPath, columns=c("TE_ID", "tot_counts", "fpkm", "sample"), as.data.table = TRUE)
RNADNA_merged <- ovlapsStatsDTFull[RNAdt, on = c("TE_ID", "RNAsampleID" = "sample"), nomatch = 0,]

fwrite(RNADNA_merged, file = file.path(opt$outDir, paste0(opt$sampleName, "_repeatsDNAme.tsv")))
}else{
  fwrite(ovlapsStatsDTFull, file = file.path(opt$outDir, paste0(opt$sampleName, "_repeatsDNAme.tsv")))
}


if (opt$runPlots) {
  message("Running plotting sections...")

# c("SINE", "Simple_repeat", "Satellite", "scRNA", "LINE")
tryCatch({
  # Attempt to create the plot
  plt_methyl <- ggboxplot(ovlapsStatsDTFull[validReadCounts == TRUE & class %in% c("SINE", "scRNA", "LINE")], 
                          x = "class", y = "fracMethyl.mean", palette = "jco") +
    stat_compare_means(method = "anova", vjust = 0.5) +      # Add global p-value
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
    labs(title = "Reference Repeats") + 
    theme(axis.text.x = element_text()) + theme_pubclean()
  
  # Save the plot
  ggsave(plt_methyl, filename = file.path(opt$outDir, paste0(opt$sampleName, "_repeats_methylation.pdf")), width = 18, height = 7)
  
  message("Plot successfully created and saved.")
}, error = function(e) {
  # Handle the error
  message("An error occurred while creating the plot: ", e$message)
})


# L1DNAme <- ovlapsStatsDTFull[subFamily %in% c("L1")][,.N, by=name]

plt_methylL1 <- ggplot(data = ovlapsStatsDTFull[validReadCounts == TRUE & clade == "L1"], aes(x=consensus, y=fracMethyl.mean)) + 
geom_boxplot() + 
labs(title = "Reference LINE-1 SubFamilies") + 
theme(axis.text.x = element_text(angle = 90))
ggsave(plt_methylL1, filename = file.path(opt$outDir, paste0(opt$sampleName, "_L1_methylation.pdf")), width = 18, height = 7)


# plt_methylL1mD <- ggplot(data = L1MD, aes(x=consensus, y=fracMethyl.mean)) + 
# geom_boxplot() + 
# labs(title = "Reference LINE-1 SubFamilies") + 
# theme(axis.text.x = element_text(angle = 90))
# ggsave(plt_methylL1mD, filename = paste0("D02_L1MD_methylation.pdf"), width = 18, height = 7)


# Visualize this
# 1. This adds a global ANOVA test p-value.
# 2. This performs pairwise t-tests comparing each subfamily against the pooled group (all data combined).
plt_methylL1mD <- ggboxplot(L1MD, x = "consensus", y = "fracMethyl.mean")+
  stat_compare_means(method = "anova", vjust=0.5) +      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")  +
labs(title = "Reference LINE-1 SubFamilies") + 
theme(axis.text.x = element_text()) + theme_pubclean()
ggsave(plt_methylL1mD, filename = file.path(opt$outDir, paste0(opt$sampleName, "_L1MD_methylation_pval.pdf")), width = 18, height = 7)




##compare all metrics of methylation
# Assume your data is already loaded as a data.table
# Example: ovlapsStatsDTFull <- fread("your_file.tsv")

# Select only columns that match 'fracMethyl.' prefix
frac_cols <- grep("^fracMethyl\\.", names(ovlapsStatsDTFull), value = TRUE)

# Filter rows with no NA in selected columns
# First, subset the full data including fracMethyl columns + grouping variable
keep_cols <- c(frac_cols, "TE_ID", "validReadCounts", "consensus")
plot_data <- ovlapsStatsDTFull[, ..keep_cols][complete.cases(ovlapsStatsDTFull[, ..frac_cols])]


L1Order <- c("L1Md_A", "L1Md_T", "L1Md_F", "L1Md_F2", "L1Md_F3", "L1Md_Gf")
plot_data$consensus <- factor(plot_data$consensus, levels = L1Order, ordered = TRUE)


# plot_data$validReadCounts <- factor(plot_data$validReadCounts,
#                                      levels = c(TRUE, FALSE),
#                                      labels = c("Valid", "NotValid"))


# Plot pairwise scatter plots
message("ggpairs for either sufficient or non-sufficient read RNA counts colored by validReadCounts")
pltGGpairs <- ggpairs(
  data = plot_data,
  columns = 1:length(frac_cols),
  title = "sufficient and non-sufficient read RNA counts",
  aes(color = validReadCounts, alpha = 0.6),
  upper = list(continuous = wrap("cor", size = 4)),
  lower = list(continuous = wrap("points", size = 0.7, alpha = 0.1)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5))
) + 
# scale_color_manual(values = c("Valid" = "gray70", "NotValid" = "dodgerblue3")) + 
theme(legend.position = "bottom") +
theme_minimal()

ggsave(pltGGpairs,
       filename = file.path(opt$outDir, paste0(opt$sampleName, "_DNA_MethylStats_ggpairs.pdf")),
       width = 18, height = 7)

ggsave(pltGGpairs,
       filename = file.path(opt$outDir, paste0(opt$sampleName, "_DNA_MethylStats_ggpairs.png")),
       width = 18, height = 9, dpi = 300)

message("ggpairs for either sufficient or non-sufficient read RNA counts")
pltGGpairsAllReadCounts <- ggpairs(
  data = plot_data,
  columns = 1:length(frac_cols),
  upper = list(continuous = wrap("cor", size = 4)),
  lower = list(continuous = wrap("points", size = 0.7, alpha = 0.1)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5))
) + 
# scale_color_manual(values = c("Valid" = "gray70", "NotValid" = "dodgerblue3")) + 
theme(legend.position = "bottom") +
theme_minimal()

ggsave(pltGGpairsAllReadCounts,
       filename = file.path(opt$outDir, paste0(opt$sampleName, "_DNA_MethylStats_ggpairsCombineValidNonValidReadCounts.png")),
       width = 18, height = 9, dpi = 300)

message("ggpairs for LINE-1 subfamilies")
all_cols <- c("consensus", frac_cols)

pltGGpairsFacetConsensus <- ggpairs(
  data = plot_data,
  columns = 1:length(frac_cols),
  aes(color = validReadCounts, alpha = 0.6),
  # upper = list(continuous = wrap("cor", size = 4)),
  upper = list(continuous = wrap("points", size = 0.7, alpha = 0.1)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5))
) + facet_wrap(~consensus) + 
theme(legend.position = "right")
#+ theme_minimal()
# scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "dodgerblue3")) + 


ggsave(pltGGpairsFacetConsensus,
       filename = file.path(opt$outDir, paste0(opt$sampleName, "_DNA_MethylStats_ggpairsConsensus.png")),
       width = 18, height = 11, dpi = 300)
} else {
    message("Skipping plotting sections as --runPlots is set to FALSE.")
}


####################################
message("\n\n done running analysis\n")
####################################
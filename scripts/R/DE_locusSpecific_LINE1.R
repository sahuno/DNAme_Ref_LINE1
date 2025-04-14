#author: Samuel ahuno
#date: march 5th 2025
#purpose: differential expression analysis of active full length L1 & locus specific L1
#1: merge readcounts from l1 and protein coding genes: out_L1_and_proteinCoding_genes_readcounts.csv
#2. use `out_L1_and_proteinCoding_genes_readcounts.csv` as input for active L1, locus specific L1
# install.packages("svglite")

# mkdir -p DE_repeats

library(optparse)
library(data.table)
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(viridis)
library(magrittr)
# library(pheatmap)
library(ggrepel)
library(ggfortify)
# library(clusterProfiler)
# library(org.Hs.eg.db)
library(ggnewscale)
library(cowplot)
library(RColorBrewer)
# library(NbClust)
# library(GenomicRanges)
# library(gUtils)
# library(janitor)
# library(fst)
library(svglite)


##QUE: do youy need different full length L1 datasets for this to run? No
##n
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/DE_locusSpecific_LINE1.R --pathMapL1BaseRepMasker "../filterSquireL1_L1Base/mapped_repeatMasker_L1Base.tsv" --FullLengthL1 TRUE --LocusSpecific TRUE


##old way
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/DE_locusSpecific_LINE1.R --pathMapL1BaseRepMasker "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mapped_repeatMasker_L1Base.tsv" --FullLengthL1 FALSE --LocusSpecific TRUE

####begin getting user inputs 
option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
        help="Print extra output [default]"),
    make_option(c("-q", "--quietly"), action="store_false",
        dest="verbose", help="Print little output"),
    make_option(c("-c", "--count"), type="integer", default=10,
        help="minimum RNA read counts [default %default]",
        metavar="number"),
    make_option(c("-L", "--LocusSpecific"), default=FALSE,
        help="run Locus Specific L1 RNA [default]"),
    make_option(c("-s", "--StatsMultipleL1"), type = "character", default="max",
        help="aggregate multiple locus specific L1 by `max`` or `geometric mean` [default]"),
    make_option(c("-F", "--FullLengthL1"), action="store_true", default=TRUE, help="run Full length L1 [default]"),
    make_option(c("-p", "--pathMapL1BaseRepMasker"), type = "character",
        help="path to bed file that maps active L1 (L1Base) with locus specific L1(SQuIRE/Repeat masker), default won't be set to avoid errors")
    )

#pathMapL1BaseRepMasker
# path_map_locusL1_activeFullLen <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mapped_repeatMasker_L1Base.tsv"

    # make_option(c("-x", "--xyz"), action="store_true", default=FALSE,
    #     help="description of option [default]")
# runLocusSpec
arguments <- parse_args(OptionParser(option_list=option_list))




options("width"=200)
#input options:
# ref_variable = "DMSO"
MIN_ReadsCounts = arguments$count
smallestGroupSize <- 3
# set_fst_threads <- 4
squireCounts_path <- "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/CT/squire_te_fwd.tsv"
counts_annot_path <- '/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/CT/counts_annot.tsv'
metadata_path <- "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_recoded.csv"
# path_map_locusL1_activeFullLen11 <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mapped_repeatMasker_AssignedTo_L1Base_mm10_mmflil1_8438.tsv"
path_map_locusL1_activeFullLen <- arguments$pathMapL1BaseRepMasker

sampleNames <- c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3","R.C.1", "R.C.2", "R.C.3", "R.Q.1", "R.Q.2", "R.Q.3","R.QC.1", "R.QC.2", "R.QC.3", "R.S.1", "R.S.2", "R.S.3","R.SC.1", "R.SC.2", "R.SC.3")
RefAlt_samples <- c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3")
DE_conditions <- data.frame(condition = c("DMSO;AZA", "DMSO;CKi", "DMSO;QSTAT", "DMSO;QSTAT-CKi"))
# threads_fst(set_fst_threads)

methodComputeMultipleLocusSpecL1ReadCounts = arguments$StatsMultipleL1
runLocusSpec = arguments$LocusSpecific
runFullLengthActiveL1 = arguments$FullLengthL1

##make ref and alt list to run a loop
REF <- "DMSO"
ALT <- "AZA"
REF_LIST <- c("DMSO", "DMSO", "DMSO", "DMSO", "DMSO","DMSO", "DMSO", "QSTAT", "CKi","SETDB1i","CKi", "QSTAT-CKi", "QSTAT")
ALT_LIST <- c("AZA", "QSTAT", "CKi", "QSTAT", "QSTAT-CKi", "SETDB1i","SETDB1i-CKi", "QSTAT-CKi", "QSTAT-CKi","SETDB1i-CKi","SETDB1i-CKi", "SETDB1i-CKi","SETDB1i")

# REF_LIST <- c("DMSO", "DMSO", "DMSO")
# ALT_LIST <- c("AZA", "QSTAT", "CKi") 

MultiDE_conditions <- data.frame(REF = REF_LIST, ALT = ALT_LIST)



##########################################################################
######## Start Analysis #################
##########################################################################

#create folder for full lengths
dir.create("LocusSpecific/figures", recursive = TRUE)
dir.create("LocusSpecific/data", recursive = TRUE)

message(paste0("\n\nrunning locusSpecific? ",arguments$LocusSpecific, "\n\n"))
# dir.create("ActiveL1/figures", recursive = TRUE)
# dir.create("ActiveL1/data", recursive = TRUE)

#read squire data
counts_Sq_dt <- fread(squireCounts_path)
counts_annot <- read_tsv(paste0(counts_annot_path))
counts_annot_pCoding <- counts_annot[counts_annot$`gene.type`=="protein_coding", ]
counts_annot_pCoding <-  counts_annot_pCoding %>% mutate_if(is.double, ~as.integer(.))

#read the metadata
metadata_df <- read.csv(file = metadata_path, sep="," ,header = TRUE)

#read overlaps data to get cordiantes of active full length L1
ovlps_counts_Sq_dt <- fread(path_map_locusL1_activeFullLen, select = c("L1UID_seqnames", "L1UID_start", "L1UID_end", "L1UID_name", "L1UID_score", "L1UID_strand", "L1UID_thickStart", "L1UID_thickEnd", "L1UID_itemRgb", "rm_name", "rm_score"))

counts_Sq_dtL1 <- counts_Sq_dt[grep(":L1:LINE",`te.name`),]
dimsL1=dim(counts_Sq_dtL1)

message(paste0("there are ",dimsL1[1], " SquireL1"))
#interger counts of repearts
intsSquire <- counts_Sq_dt[, lapply(.SD, as.integer), .SDcols = !c("te.id", "te.name")]

# countsSepNameCols <- counts_Sq_pCoding_df_df %>% separate_wider_delim(te.name, "|" , names = c("chrom", "start", "end", "name","number","strand")) %>% separate_wider_delim(name, ":" ,names = c("consensus","clade", "class"))


##### remove all rows with less than 10 reads, 
keepDF <- rowSums(intsSquire >= MIN_ReadsCounts) >= smallestGroupSize
keepCodingDF <- rowSums(counts_annot_pCoding[, sampleNames] >= MIN_ReadsCounts) >= smallestGroupSize

PercentageRepKept <- (sum(keepDF)/length(keepDF))*100
PercentageCodingKept <- (sum(keepCodingDF)/length(keepCodingDF))*100

message("Percentage of repeats kept after filtering: ", PercentageRepKept, "%")
message("Percentage of protein coding genes kept after filtering: ", PercentageCodingKept, "%")


#remove rows with less than 10 reads
counts_Sq_ints_dt <- cbind(counts_Sq_dt[keepDF,c("te.id", "te.name")],intsSquire[keepDF,])
counts_Sq_ints_dt[,`gene.type` := "repeats"] ## repeats id
counts_Sq_ints_dt[,c("chr", "start", "end", "name","number","strand") := tstrsplit(`te.name`, "|", fixed = TRUE)][, c("consensus","clade", "class"):=tstrsplit(name, ":", fixed = TRUE)] ##split the te.name column


##Raw Values of valid L1
# table(counts_Sq_ints_dt$clade) # L1(N=3652)

counts_Sq_ints_dt[, c("start", "end"):= lapply(.SD, as.integer), .SDcols = c("start", "end")]
##merge protein coding data and  repeats, all counts matrix become numeric 
counts_Sq_pCoding_df <- rbind(counts_Sq_ints_dt, counts_annot_pCoding[keepCodingDF,], fill=TRUE)

#create ColumnKey for repeats and protein coding 
counts_Sq_pCoding_df[,genomicElementID:=(fcase(!is.na(te.name), te.name,
                                                    !is.na(gene.id), gene.id ) )]

counts_Sq_pCoding_df <- counts_Sq_pCoding_df %>% relocate(all_of(sampleNames), .after = last_col())


#create ColumnKey for repeats and protein coding 


##########################################
#### DE; all locus specific repeats
##########################################
#step 1: convert to dataframe
#step 2: metadata
#step 3: DE analysis
counts_Sq_pCoding_df_df <- as.data.frame(counts_Sq_pCoding_df)
rownames(counts_Sq_pCoding_df_df) <- counts_Sq_pCoding_df_df$genomicElementID

# counts_Sq_pCoding_df[grepl("protein_coding", `gene.type`) & grepl("L1", genomicElementID),] #is there L1 names in protein coding genes?  


#needed functions
############################
#volcano plots of line 1
##FUN
plotEVolcano <- function(df_input, title_in, subtitle_in, selectLab_in=NULL, labs=NULL, REF, ALT){
    if(is.null(labs)){
      labs = rownames(df_input)
    }
    volcano_ggplot <- EnhancedVolcano(df_input,
    #lab = rownames(df_input),
    lab = labs,
    #selectLab = selectLab_in,
    x = paste0('log2FoldChange(',REF,"/", ALT,")"),
    y = 'padj',
    pCutoff = 10e-2,
    title = title_in,
    subtitle = subtitle_in,
    col=c('grey', 'grey', 'grey', 'blue'),
    ylab = bquote(~-Log[10]~ 'padjust'),
    gridlines.major = FALSE,
    gridlines.minor = FALSE)
}
#subtitle_in="Differential locus Specific L1 expression"
##FUN
makeCPM <- function(DeSeq2_res, prior_cnts=2, keyColName){
#get expression matrix
DeSeq2_res_cpm <- edgeR::cpm(DeSeq2_res, log = TRUE, prior.count = 2) #prior.count = 2 is added to avoid log(0) error
#message("computed cpm")
DeSeq2_res_df <- as.data.frame(DeSeq2_res_cpm)
DeSeq2_res_df <- DeSeq2_res_df %>% tibble::rownames_to_column(var = keyColName)
return(DeSeq2_res_df)
}

## compute the geometric mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}







# #drop samples if necesssary
# conditions_Default <- c("DMSO", "AZA")
# cond_collapse <- paste0(c("condition",conditions_Default), collapse="_")

# metadata_df_short <- metadata_df %>% filter(str_detect(samples, "^R.A|^R.0")) #%>% mutate(condition)
# metadata_df_short <- metadata_df_short %>% dplyr::filter(condition %in% conditions_Default)

# ##change this for multiple conditions
# conditions_inUse <- conditions_Default
# #################################
# #################################
# # run differenrial gene expression
# dds <- DESeqDataSetFromMatrix(countData = data.matrix(counts_Sq_pCoding_df_df[,RefAlt_samples]),
#                               colData = metadata_df_short,
#                               design = ~ condition)

# rowData(dds) <- counts_Sq_pCoding_df_df[, c("gene.type","length","gene.symbol","ensembl.gene.id", "gene.id", "name","consensus")]

# #normalise repeats with non-repeats
# dds <- DESeq(dds)

# ### filter L1 only for DE analysis
# dds_L1 <- dds[str_detect(rownames(dds),"L1"),]
# dds_Repeats <- dds[str_detect(rowData(dds)$`gene.type`,"repeats"),]
# dds_pCoding <- dds[str_detect(rowData(dds)$`gene.type`,"protein_coding"),]


# ddsL1Results <- results(dds_L1)
# ddsRepeatsResults <- results(dds_Repeats)
# ddsProteinCodingResults <- results(dds_pCoding)


# #summary of DE results
# ddsL1ResultsDf <- as.data.frame(ddsL1Results)
# summary(ddsResults)

# #shrink results for volcano
# results_L1Shrink <- lfcShrink(dds_L1, contrast = c("condition", conditions_inUse), res=ddsL1Results, type = 'normal')



# # volcano_ggplot <- EnhancedVolcano(results_L1Shrink,
# #     lab = rownames(results_L1Shrink),
# #     x = 'log2FoldChange',
# #     y = 'padj',
# #     title = '5-AZA versus DMSO',
# #     subtitle = "Differential locus Specific L1 expression",
# #     col=c('grey', 'grey', 'grey', 'blue'),
# #     ylab = bquote(~-Log[10]~ 'padjust'))


# ggsave(plotEVolcano(results_L1Shrink, title_in='5-AZA vrs DMSO'), filename = "volcano_ggplot_locusSpecificLINE1.png")






# # #get expression matrix
# # dds_sf_Disp_L1_cpm <- edgeR::cpm(dds_L1, log = TRUE, prior.count = 2) #prior.count = 2 is added to avoid log(0) error
# # dds_sf_Disp_L1_cpm_df <- as.data.frame(dds_sf_Disp_L1_cpm)
# # dds_sf_Disp_L1_cpm_df <- dds_sf_Disp_L1_cpm_df %>% rownames_to_column(var = "te.id")

# dds_sf_Disp_L1_cpm_df <- makeCPM(DeSeq2_res=dds_L1, prior_cnts=2, keyColName="te.id")

# #save results
# results2Save <- list(dds, dds_L1, ddsL1Results, results_L1Shrink, dds_sf_Disp_L1_cpm_df)
# save(results2Save, file=paste0("Differential_locusSpecific_L1ExpressionResults_",cond_collapse,".RData"))  # Till here everything is ok. 

######################
##do for
######################
# function()





##########################################
#### DE; function to run loop
##########################################
Repeat_DE <- function(data_in, metadata_in, REF, ALT){
cond_collapsed <- paste0(c("condition",REF, ALT), collapse="_")
#get metadata
metadata_short <- metadata_in %>% filter(condition %in% c(REF, ALT))
metadata_short <- metadata_short %>% mutate(condition = factor(condition, levels = c(REF, ALT))) #factor for DESeq2

ddsInFunc <- DESeqDataSetFromMatrix(countData = data.matrix(data_in[,metadata_short$samples]),
                              colData = metadata_short,
                              design = ~ condition)

rowData(ddsInFunc) <- data_in[, c("gene.type","length","gene.symbol","ensembl.gene.id", "gene.id", "name","consensus")]
#normalise repeats with non-repeats
dds_main <- DESeq(ddsInFunc)

### filter out repeats of interest
dds_L1 <- dds_main[str_detect(rownames(dds_main),"L1"),]
dds_Repeats <- dds_main[str_detect(rowData(dds_main)$`gene.type`,"repeats"),]
dds_pCoding <- dds_main[str_detect(rowData(dds_main)$`gene.type`,"protein_coding"),]

ddsL1Results <- results(dds_L1)
ddsRepeatsResults <- results(dds_Repeats)
ddsProteinCodingResults <- results(dds_pCoding)

#summary of DE results
ddsL1ResultsDf <- as.data.frame(ddsL1Results)
ddsRepeatsResultsDf <- as.data.frame(ddsRepeatsResults)

message("summary of locus-Spec L1 DE results")
message(summary(ddsL1Results))
message("/n")

#shrink results for volcano
resultsL1Shrink <- lfcShrink(dds_L1, contrast = c("condition", REF, ALT), res=ddsL1Results, type = 'normal')
resultsRepeatsShrink <- lfcShrink(dds_Repeats, contrast = c("condition", REF, ALT), res=ddsRepeatsResults, type = 'normal')
# resultsProteinCodingShrink <- lfcShrink(dds_pCoding, contrast = c("condition", REF, ALT), res=ddsProteinCodingResults, type = 'normal')

message("making volcano plots")
volL1NoShrink <- plotEVolcano(ddsL1Results, title_in=cond_collapsed, subtitle_in="Differential Locus Specific L1 Expression")
volL1Shrink <- plotEVolcano(resultsL1Shrink, title_in=cond_collapsed, subtitle_in="Differential Locus Specific L1 Expression")

volRepeatsShrink <- plotEVolcano(resultsRepeatsShrink, title_in=cond_collapsed, subtitle_in="Differential Locus Specific Repeats Expression")
ggsave(volL1Shrink, filename = paste0("LocusSpecific/figures/volcanoShrink_locusSpecificLINE1_",cond_collapsed,".svg"))
ggsave(volL1NoShrink, filename = paste0("LocusSpecific/figures/volcanoNoShrink_locusSpecificLINE1_",cond_collapsed,".svg"))
ggsave(volRepeatsShrink, filename = paste0("LocusSpecific/figures/volcanoShrink_locusSpecificRepeats_",cond_collapsed,".svg"))
# ggsave(plotEVolcano(resultsProteinCodingShrink, title_in=cond_collapsed), filename = "volcano_proteinCoding_",cond_collapsed,".png")

#make CPM
message("making CPM")
dds_locusSpecL1_cpm_df <- makeCPM(DeSeq2_res=dds_L1, prior_cnts=2, keyColName="te.id")
dds_repeats_cpm_df <- makeCPM(DeSeq2_res=dds_Repeats, prior_cnts=2, keyColName="te.id")

#save results
message("saving results to disk")
results2SaveLspecL1 <- list(ddsMain=dds_main, ddsFiltered=dds_L1, resDF=ddsL1Results, resShrinkedDF=resultsL1Shrink, cpmDF=dds_locusSpecL1_cpm_df, volcano_ggplot = volL1Shrink)
results2SaveLspecRepeats <- list(ddsMain=dds_main, ddsFiltered=dds_Repeats, resDF=ddsRepeatsResultsDf, resShrinkedDF=resultsRepeatsShrink, cpmDF=dds_repeats_cpm_df, volcano_ggplot = volRepeatsShrink)
save(results2SaveLspecL1, file=paste0("LocusSpecific/data/Differential_locusSpecific_L1ExpressionResults_",cond_collapsed,".RData"))  # Till here everything is ok. 
# return(ddsRepeatsResults)
}


#### run the function
# metadata_df

# ddsRepeatsResults <- Repeat_DE(data_in = counts_Sq_pCoding_df_df, metadata_in = metadata_df, REF = REF, ALT = ALT)
# Repeat_DE(data_in = counts_Sq_pCoding_df_df, metadata_in = metadata_df, REF = REF, ALT = ALT)

# MultiDE_conditions; locus specific analysis
if(runLocusSpec == TRUE){
lapply(1:nrow(MultiDE_conditions), function(x) {Repeat_DE(data_in = counts_Sq_pCoding_df_df, metadata_in = metadata_df, REF = MultiDE_conditions[x,1], ALT = MultiDE_conditions[x,2])})
}else{
message("skipping locus specific analysis")
}
# DeSeq2res1 <- edgeR::cpm(ddsRepeatsResults, log = TRUE, prior.count = 2) #prior.count = 2 is added to avoid log(0) error


##########################################
#### DE; full length L1
##########################################
# mapSquireL1L1Base <- fread(path_map_locusL1_activeFullLen)
message("\n\n preparing data for full length L1 analysis\n\n")
dfRepL1 <- fread("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mapped_overlaps_repeatMasker_L1Base_mm10.tsv")
# dfRepL1[L1UID_name == "UID-828"]
# dfRepL1[L1UID_name == "UID-827"] 
# countsSepNameCols <- counts_Sq_pCoding_df_df %>% separate_wider_delim(te.name, "|" , names = c("chrom", "start", "end", "name","number","strand")) %>% separate_wider_delim(name, ":" ,names = c("consensus","clade", "class"))
# countsSepNameCols <- counts_Sq_pCoding_df_df %>% separate_wider_delim(te.name, "|" , names = c("chrom", "start", "end", "name","number","strand")) %>% separate_wider_delim(name, ":" ,names = c("consensus","clade", "class"))

# dfRepL1[, paste0("seqnames", "start", "end",  "rm_name")]
dfRepL1[, idCol:=paste0("seqnames", "start", "end",  "rm_name", "rm_score", "strand")]
dfRepL1 <- dfRepL1[, idCol:=paste0(seqnames, "|",start,"|" ,end, "|", rm_name, "|",rm_score, "|" ,strand)][,.SD, .SDcols = c("idCol","L1UID_name")]

# dfRepL1[idCol == "chr15|55677600|55684869|L1Md_T:L1:LINE|30|-", ]



###### get active full length L1 and merge with repeatmasker
#filter out repeats `counts_Sq_pCoding_df` 
# merge l1 active and repeatmasker.
# re-add merged l1 active and rm to `counts_Sq_pCoding_df`

# copy(counts_Sq_pCoding_df)[str_detect(name, "L1"),][mapSquireL1L1Base, on = .(name=rm_name)]
# mapSquireL1L1Base

# counts_Sq_pCoding_df[name == "L1_Mus3:L1:LINE"]
# dfRepL1[str_detect(idCol, "L1_Mus3:L1:LINE"),]

# separate out L1
validReadCountsRepeatsOnlyDT <- counts_Sq_pCoding_df[grepl("repeats", `gene.type`), ]
validReadCountsLociSpecL1RepeatsOnlyDT <- counts_Sq_pCoding_df[grepl("L1", name),]
# sum(duplicated(validReadCountsLociSpecL1RepeatsOnlyDT$genomicElementID))
# [dfRepL1, on = .(`te.name`=`idCol`), nomatch=0L]
validReadCountsRepeatsOnlyDT[genomicElementID == "chr15|55677600|55684869|L1Md_T:L1:LINE|30|-", ]

validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT <- validReadCountsLociSpecL1RepeatsOnlyDT[dfRepL1, on = .(`te.name`=`idCol`), nomatch=0L]

#save a copy for review
validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT_withCordinates <- validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT[ovlps_counts_Sq_dt, on = .(L1UID_name), nomatch=0L]
fwrite(validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT_withCordinates, file="validReadCountsLociSpecL1_overlapping_withFullLengthL1.tsv", sep="\t", col.names = TRUE)
#################
#plot #lsL1 per full length
numberLsL1_perFullLengthDT <- validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT[ , .N , by = L1UID_name]
plotLS_number <- ggplot() + geom_bar(data = numberLsL1_perFullLengthDT, aes(x = N), fill = "blue") + theme_minimal() + labs(title = "number of locus specific L1 per full length L1", x = "Active Full Length L1", y = "Frequency of locus specific L1")
ggsave(plotLS_number, filename = "plot_number_LocusSpecL1_perFulllength.svg")


#duplicate: chr15|55677600|55684869|L1Md_T:L1:LINE|30|-
##add width  of TE
validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT[, L1width := (end - start) + 1]
validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT[, idx := .I]
validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT[, isMaxWidth:=(L1width == max(L1width)), by = L1UID_name]
validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT[,`gene.type` := "activeL1"] #for filtering, do it earlier


if(methodComputeMultipleLocusSpecL1ReadCounts == "max") {
# To get other columns corresponding to these maximum rows, join back to the original data.table:
gmMean_validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT <- validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT[isMaxWidth == TRUE,]
} else if (methodComputeMultipleLocusSpecL1ReadCounts == "geometericMean") {
   #compute geometeric mean of read counts to combine multiple L1 that overlaps full length L1
gmMean_validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT <- validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT[, lapply(.SD, function(x){(as.integer(gm_mean(x)))}),by = L1UID_name, .SDcols = sampleNames]
}

#sanity checks, is active l1 id duplicated?
sum(duplicated(gmMean_validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT$L1UID_name))

####TO-DO
# merge 1. protein coding, non LINE1, 3. geometeric mean of line1
counts_Sq_pCoding_df[!grepl("protein_coding", `gene.type`),] #protein coding only  
counts_Sq_pCoding_df[!grepl("protein_coding", `gene.type`) & !grepl("L1", name),] #not protein coding but some other repeats except L1  
counts_Sq_pCoding_df[grepl("protein_coding", `gene.type`) & grepl("L1", `te.name`),] #is there L1 names in protein coding genes?  
counts_Sq_pCoding_df[grepl("protein_coding", `gene.type`) & grepl("L1", genomicElementID),] #is there L1 names in protein coding genes?  

# gmMean_validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT
# table(counts_Sq_pCoding_df$`gene.type`)
ptCoding_allRepeats_exceptL1 <- rbind(counts_Sq_pCoding_df[!grepl("repeats", `gene.type`), ], counts_Sq_pCoding_df[!grepl("protein_coding", `gene.type`) & !grepl("L1", name),])
sum(duplicated(ptCoding_allRepeats_exceptL1$genomicElementID))
sum(is.na(ptCoding_allRepeats_exceptL1$genomicElementID))

# table(gmMean_validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT$genomicElementID) #these are the LS_L1 that make up the UID
###this is what i needed for the final table; merge protin coding, all Repeats(except locus specific L1) & active L1 for DeSeq2 Normalization and later filtering
ptCoding_allRepeats_activeL1 <- rbind(ptCoding_allRepeats_exceptL1, gmMean_validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT, fill=TRUE)
sum(duplicated(ptCoding_allRepeats_activeL1$genomicElementID) == TRUE)
table(gmMean_validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT$`gene.type`)

idxDuplicatedRows <- which(duplicated(ptCoding_allRepeats_activeL1$genomicElementID)==TRUE)
idDuplicatedRows <- ptCoding_allRepeats_activeL1[idxDuplicatedRows, genomicElementID]
ptCoding_allRepeats_activeL1[genomicElementID %in% idDuplicatedRows, ]

validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT[genomicElementID %in% idDuplicatedRows,]

# make unique genomicElementID for DESeq2 rownames
ptCoding_allRepeats_activeL1[, genomicElementID := make.unique(genomicElementID)]
# ptCoding_allRepeats_activeL1[, genomicElementIDUnique := make.unique(genomicElementID)]
message("below are the same RepmaskerL1 that maps to full Length L1")
ptCoding_allRepeats_activeL1[genomicElementID %in% idDuplicatedRows,]

#get actual rows with duplicated id
ptCoding_allRepeats_activeL1[duplicated(ptCoding_allRepeats_activeL1$genomicElementID),]
# table(ptCoding_allRepeats_activeL1$`gene.type`)

# sum(is.na(ptCoding_allRepeats_activeL1$genomicElementID))

# ptCoding_allRepeats_activeL1[is.na(`gene.type`), `gene.type` := "activeL1"] #for filtering #change this!
sum(is.na(ptCoding_allRepeats_activeL1$`gene.type`))
sum(is.na(ptCoding_allRepeats_activeL1$genomicElementID))
#if `gene.type` := "activeL1" then Use the UID-naming convention
ptCoding_allRepeats_activeL1[, genomicElementID := fcase(`gene.type` == "activeL1", L1UID_name, default = genomicElementID)]#[`gene.type` == "activeL1",]
# ptCoding_allRepeats_activeL1[is.na(genomicElementID), genomicElementID := L1UID_name] #if use genomicElementID is na replace with L1UID for deseq rownmaes `genomicElementID` for rownames
# ptCoding_allRepeats_activeL1[duplicated(ptCoding_allRepeats_activeL1$genomicElementID),] #

sum(duplicated(ptCoding_allRepeats_activeL1$genomicElementID)) #this should be unique
sum(is.na(ptCoding_allRepeats_activeL1$genomicElementID))
sum(!is.na(ptCoding_allRepeats_activeL1$L1UID_name))

###save a copy of main data frame used for deseq
fwrite(ptCoding_allRepeats_activeL1, file=paste0("validReadCounts_ptCoding_LociSpecL1_ActiveL1_allOtherRepeats.tsv"), sep="\t", col.names = TRUE)
# fwrite(validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT_withCordinates, file="validReadCountsLociSpecL1_overlapping_withFullLengthL1.tsv", sep="\t", col.names = TRUE)

###prepare data for active line 1 DE run
message("make matrix data for active L1")
colsDE <- c("genomicElementID", "gene.type", "L1UID_name",sampleNames) #use this for the final table
ptCoding_allRepeats_activeL1_DF <- as.data.frame(ptCoding_allRepeats_activeL1[, ..colsDE])
rownames(ptCoding_allRepeats_activeL1_DF) <- ptCoding_allRepeats_activeL1_DF$genomicElementID
# head(ptCoding_allRepeats_activeL1_DF)

#sanity checks, could you filter
#sum(str_detect(ptCoding_allRepeats_activeL1_DF$genomicElementID,"UID-"))#
#sanity checks
# ptCoding_allRepeats_activeL1[!is.na(genomicElementID),]

# ptCoding_allRepeats_activeL1[is.na(`gene.type`), `gene.type` := "activeL1"] #for filtering
# counts_Sq_pCoding_df[,genomicElementID:=(fcase(!is.na(te.name), te.name,
                                                    # !is.na(gene.id), gene.id ) )]

#merge locus specific l1 
# validReadCountsRepeatsOnlyDT[!grepl("L1", name), .(idCol = genomicElementID, L1UID_name = NA),]
#rbind(validReadCountsRepeatsOnlyDT[!grepl("L1", name), .(idCol = genomicElementID, L1UID_name = NA),], validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT,fill=TRUE)
# validReadCountsRepeatsOnlyDT[grepl("L1", name),][dfRepL1,on = .(genomicElementID=idCol) ]

# validReadCountsRepeatsOnlyDT[grepl("L1", name), .(idCol = genomicElementID, L1UID_name = NA),][dfRepL1,on = .(genomicElementID=idCol) ]


# counts_Sq_pCoding_active_pre <- copy(counts_Sq_pCoding_df)[str_detect(name, "L1"),][dfRepL1, on = .(`te.name`=`idCol`), nomatch=0L]
# copy(counts_Sq_pCoding_df)[str_detect(name, "L1"),][dfRepL1[,c("seqnames", "start", "end", "rm_name", "L1UID_name")], on = .(chr=seqnames, start=start, end=end)]
# counts_Sq_pCoding_active_pre[, lapply(.SD, function(x){as.integer(mean(x))}),by = L1UID_name, .SDcols = c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3")]
# counts_Sq_pCoding_active_pre[, lapply(.SD, function(x){(as.integer(mean(x)))}),by = L1UID_name, .SDcols = c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3")]

# counts_Sq_pCoding_active_pre[, lapply(.SD, function(x){(as.integer(gm_mean(x)))}),by = L1UID_name, .SDcols = c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3")]

# "R.0.1 R.0.2 R.0.3 R.A.1 R.A.2 R.A.3"
# dtDeseq <- ptCoding_allRepeats_activeL1_DF[, c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3")]
# data_in <- ptCoding_allRepeats_activeL1_DF[, c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3")]
# any(is.na(dtDeseq))

activeL1DE <- function(data_in, metadata_in, REF, ALT){
cond_collapsed <- paste0(c("condition", REF, ALT), collapse="_")
#get metadata
metadata_short <- metadata_in %>% filter(condition %in% c(REF, ALT))
metadata_short <- metadata_short %>% mutate(condition = factor(condition, levels = c(REF, ALT))) #factor for DESeq2

message("making DESeq2 object")
ddsInFunc <- DESeqDataSetFromMatrix(countData = data.matrix(data_in[,metadata_short$samples]),
                              colData = metadata_short,
                              design = ~ condition)

rowData(ddsInFunc) <- data_in[, c("L1UID_name", "gene.type")]

message("checking for NAs")
message(any(is.na(counts(ddsInFunc))))

message("normalise repeats with non-repeats")
dds_main <- DESeq(ddsInFunc)
### filter out repeats of interest

message("filter out repeats of interest")
dds_activeL1 <- dds_main[str_detect(rownames(dds_main),"UID-"),]
resActiveL1 <- results(dds_activeL1)
resActiveL1Df <- as.data.frame(resActiveL1)

message("summary of active L1 DE results")
message(summary(resActiveL1))
message("/n")

# shrink results for volcano
resultsActiveL1Shrink <- lfcShrink(dds_activeL1, contrast = c("condition", REF, ALT), res=resActiveL1, type = 'normal')

message("making volcano plots")
volActiveL1Shrink <- plotEVolcano(resultsActiveL1Shrink, title_in=cond_collapsed, subtitle_in="Differential Active L1 Expression")
volActiveL1NoShrink <- plotEVolcano(resActiveL1, title_in=cond_collapsed, subtitle_in="Differential Active L1 Expression")

ggsave(volActiveL1Shrink, filename = paste0("ActiveL1/figures/volcanoShrink_ActiveLINE1_",cond_collapsed,".svg"))
ggsave(volActiveL1NoShrink, filename = paste0("ActiveL1/figures/volcanoNoShrink_ActiveLINE1_",cond_collapsed,".svg"))


# make CPM
message("making CPM")
ActiveL1_cpm_df <- makeCPM(DeSeq2_res=dds_activeL1, prior_cnts=2, keyColName="genomicElementID")

#save results
message("saving results to disk")
results2SaveActiveL1 <- list(ddsMain=dds_main, ddsFiltered=resActiveL1, resDF=resActiveL1Df, resShrinkedDF=resultsActiveL1Shrink, cpmDF=ActiveL1_cpm_df, volcano_ggplot = volActiveL1Shrink)
save(results2SaveActiveL1, file=paste0("ActiveL1/data/Differential_ActiveL1ExpressionResults_",cond_collapsed,".RData"))  # Till here everything is ok. 
}



# tstrsplit("chr1|100|200|L1|1|+", "|", fixed = TRUE)[4]
# tstrsplit(head(counts_Sq_pCoding_df$genomicElementID, 10), "|", fixed = TRUE)[4][[1]]

# counts_Sq_pCoding_df$genomicElementID[grepl("L1", counts_Sq_pCoding_df$genomicElementID)]
if(runFullLengthActiveL1 == TRUE){
  message("running locus specific l1")
  dir.create("ActiveL1/figures", recursive = TRUE)
  dir.create("ActiveL1/data", recursive = TRUE)
lapply(1:nrow(MultiDE_conditions), function(x) {activeL1DE(data_in = ptCoding_allRepeats_activeL1_DF, metadata_in = metadata_df, REF = MultiDE_conditions[x,1], ALT = MultiDE_conditions[x,2])})
}else{
  message("skipping full length L1 analysis")
}


###############################
##### detiled 
##############################
#how many full lenghts have >1 locus-specific l1 mapping
numberLsL1_perFullLengthDT[N > 1,] 
#
ptCoding_allRepeats_activeL1_DF[str_detect(rownames(ptCoding_allRepeats_activeL1_DF),"UID-"),]

multipleLocusSpecReadsDT <- validReadCountsLociSpecL1RepeatsOnly_withL1BaseID_DT[L1UID_name %in% numberLsL1_perFullLengthDT[N > 1, L1UID_name],]
multipleLocusSpecReadsList <- split(multipleLocusSpecReadsDT[,..colsDE ], multipleLocusSpecReadsDT$L1UID_name)

multipleLocusSpecReadsDT_long <- melt(multipleLocusSpecReadsDT, id.vars = c("genomicElementID", "L1UID_name", "gene.type"), measure.vars = sampleNames, variable.name = "samples", value.name = "readCounts")
multipleLocusSpecReadsDT_longCnds <- multipleLocusSpecReadsDT_long[metadata_df, on=.(samples), nomatch=0L]
setnames(multipleLocusSpecReadsDT_longCnds, "genomicElementID", "ID")



multipleLocusSpecReadsDT_longCndsLs <- split(multipleLocusSpecReadsDT_longCnds, multipleLocusSpecReadsDT_longCnds$L1UID_name)
plot_readCounts <- ggplot(data = multipleLocusSpecReadsDT_longCndsLs[['UID-828']], aes(x=samples, y=readCounts)) + geom_point(aes(color=ID)) + facet_wrap(~condition, scales = "free")  + theme_minimal() + theme(legend.position = "bottom")+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_color_discrete(labels = function(x) stringr::str_wrap(x, width = 20))
ggsave(plot_readCounts, filename = "readCounts_activeL1.svg")

# 
plot_multipleReadCounts <- function(namesList){
  plot_readCounts <- ggplot(data = multipleLocusSpecReadsDT_longCndsLs[[namesList]], aes(x=samples, y=readCounts)) + geom_point(aes(color=ID)) + facet_wrap(~condition, scales = "free")  + theme_minimal() + theme(legend.position = "bottom")+labs(title = namesList) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_color_discrete(labels = function(x) stringr::str_wrap(x, width = 20))
return(plot_readCounts)
}


pdf("readCounts_activeL1_list.pdf")
lapply(names(multipleLocusSpecReadsDT_longCndsLs), function(x) {plot_multipleReadCounts(x)})
dev.off()



#######################################################
#### get and plot stats
############################################################
# multipleLocusSpecReadsDT_longCndsLs[["UID-181"]][,.(readCounts_mean = mean(readCounts), readCounts_var = var(readCounts), readCounts_geomMean=gm_mean(readCounts), readCounts_SD=sd(readCounts) ), by=samples]

multipleLocusSpecReadsDT_stats_longCndsLs <- lapply(multipleLocusSpecReadsDT_longCndsLs, function(x){
  x[,.(readCounts_mean = mean(readCounts), readCounts_var = var(readCounts), readCounts_geomMean=gm_mean(readCounts), readCounts_SD=sd(readCounts) ), by=samples][metadata_df, on=.(samples), nomatch=0L]
})

plot_StatsMultipleReadCounts <- function(namesList){
  plot_readCountsStats <- ggplot(data = multipleLocusSpecReadsDT_stats_longCndsLs[[namesList]], aes(x=readCounts_mean, y=readCounts_SD)) + geom_point(aes(color=samples)) + facet_wrap(~condition, scales = "free")  + theme_minimal() + theme(legend.position = "bottom")+ 
  labs(title = namesList) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_color_discrete(labels = function(x) stringr::str_wrap(x, width = 20))
return(plot_readCountsStats)
}

pdf("stats_readCounts_activeL1_list.pdf")
lapply(names(multipleLocusSpecReadsDT_stats_longCndsLs), function(x) {plot_StatsMultipleReadCounts(x)})
dev.off()


# multipleLocusSpecReadsDT_longCndsLs[["UID-181"]][,.(readCounts_mean = mean(readCounts), readCounts_var = var(readCounts), readCounts_geomMean=gm_mean(readCounts), readCounts_SD=sd(readCounts) ), by=samples]



# ddsActiveL1Main <- activeL1DE(data_in=ptCoding_allRepeats_activeL1_DF, metadata_in=metadata_df, REF=REF, ALT=ALT)
# sum(duplicated(ptCoding_allRepeats_activeL1$genomicElementID))

# dds_activeL1 <- ddsActiveL1Main[str_detect(rownames(ddsActiveL1Main),"UID-"),]
# resActiveL1 <- results(dds_activeL1)
# resActiveL1Df <- as.data.frame(resActiveL1)

# ptCoding_allRepeats_activeL1[!is.na(L1UID_name), c("R.0.1", "L1UID_name")]
# resActiveL1$rowname <- rownames(resActiveL1)
# merge(resActiveL1, ptCoding_allRepeats_activeL1[!is.na(L1UID_name), c("R.0.1", "L1UID_name")], by.x = "rowname", by.y = "L1UID_name")

# message("summary of active L1 DE results")
# message(summary(resActiveL1))
# message("/n")

#shrink results for volcano
# resultsActiveL1Shrink <- lfcShrink(dds_activeL1, contrast = c("condition", REF, ALT), res=resActiveL1, type = 'normal')

# message("making volcano plots")
# volActiveL1Shrink <- plotEVolcano(resultsActiveL1Shrink, title_in=cond_collapsed)
# ggsave(volActiveL1Shrink, filename = paste0("volcano_ActiveLINE1_",cond_collapsed,".png"))

#make CPM
# message("making CPM")
# ActiveL1_cpm_df <- makeCPM(DeSeq2_res=dds_activeL1, prior_cnts=2, keyColName="genomicElementID")

# #save results
# message("saving results to disk")
# results2SaveActiveL1 <- list(ddsMain=ddsActiveL1Main, ddsFiltered=resActiveL1, resDF=resActiveL1Df, resShrinkedDF=resultsActiveL1Shrink, cpmDF=ActiveL1_cpm_df, volcano_ggplot = volActiveL1Shrink)

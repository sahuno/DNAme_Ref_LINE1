#here, scripts locus specfic line1 Analysis
#input: total counts table, metadata, 
#outputs: volcano plots, rdata DEResults, pca of repeat 



#author: Samuel ahuno
#date: April 3rd 2025
#purpose: differential expression analysis of locus specific L1
#1: merge readcounts from l1 and protein coding genes: out_L1_and_proteinCoding_genes_readcounts.csv
#2. use `out_L1_and_proteinCoding_genes_readcounts.csv` as input for active L1, locus specific L1


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
library(fst)
library(svglite)


##QUE: do youy need different full length L1 datasets for this to run? No
##n
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/DE_locusSpecificLINE1.R --LocusSpecific TRUE


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
        help="aggregate multiple locus specific L1 by `max`` or `geometric mean` [default]")
    )

#pathMapL1BaseRepMasker
# path_map_locusL1_activeFullLen <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mapped_repeatMasker_L1Base.tsv"

# dt1 <- fread("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/SQuireRepeatsValidReadCounts_Aggregated_FWD_REV_totCounts.bed")
arguments <- parse_args(OptionParser(option_list=option_list))

options("width"=500)

#input options:
# ref_variable = "DMSO"
MIN_ReadsCounts = arguments$count
smallestGroupSize <- 3
# set_fst_threads <- 4
readCountsPath <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/squireRepeatsAnnot_Aggregated_FWD_REV_totalCounts.fst"
counts_annot_path <- '/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/CT/counts_annot.tsv'
metadata_path <- "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_recoded.csv"
path_map_locusL1_activeFullLen <- arguments$pathMapL1BaseRepMasker

sampleNames <- c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3","R.C.1", "R.C.2", "R.C.3", "R.Q.1", "R.Q.2", "R.Q.3","R.QC.1", "R.QC.2", "R.QC.3", "R.S.1", "R.S.2", "R.S.3","R.SC.1", "R.SC.2", "R.SC.3")
RefAlt_samples <- c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3")
DE_conditions <- data.frame(condition = c("DMSO;AZA", "DMSO;CKi", "DMSO;QSTAT", "DMSO;QSTAT-CKi"))
# threads_fst(set_fst_threads)

methodComputeMultipleLocusSpecL1ReadCounts = arguments$StatsMultipleL1
runLocusSpec = arguments$LocusSpecific


##make ref and alt list to run a loop
REF <- "DMSO"
ALT <- "AZA"
REF_LIST <- c("DMSO", "DMSO", "DMSO",  "DMSO",   "DMSO",       "DMSO",        "SETDB1i",     "QSTAT")
ALT_LIST <- c("AZA",  "CKi",  "QSTAT", "SETDB1i", "QSTAT-CKi", "SETDB1i-CKi",  "SETDB1i-CKi", "QSTAT-CKi")

# REF_LIST <- c("DMSO", "DMSO", "DMSO", "DMSO", "DMSO","DMSO", "DMSO", "QSTAT", "CKi","SETDB1i","CKi", "QSTAT-CKi", "QSTAT")
# ALT_LIST <- c("AZA", "QSTAT", "CKi", "QSTAT", "QSTAT-CKi", "SETDB1i","SETDB1i-CKi", "QSTAT-CKi", "QSTAT-CKi","SETDB1i-CKi","SETDB1i-CKi", "SETDB1i-CKi","SETDB1i")


# REF_LIST <- c("DMSO", "DMSO")
# ALT_LIST <- c("AZA", "QSTAT")

# REF_LIST <- c("DMSO", "DMSO", "DMSO")
# ALT_LIST <- c("AZA", "QSTAT", "CKi") 

MultiDE_conditions <- data.frame(REF = REF_LIST, ALT = ALT_LIST)



##########################################################################
######## Start Analysis #################
##########################################################################
#read in the data
countsMatrixDT <- read_fst(readCountsPath, as.data.table = TRUE)
names(countsMatrixDT) <- gsub("-", ".", names(countsMatrixDT))

validCountsMatrixRepsDT <- countsMatrixDT[validReadCounts == TRUE,]
validCountsMatrixRepsDT[, `gene.type` := class] #for id when merging with protein coding genes
#create ColumnKey for repeats and protein coding 
validCountsMatrixRepsDT[,RNAelementID := `TE_ID`] #for id when merging with protein coding genes

validCountsMatrixL1DT <- countsMatrixDT[validReadCounts == TRUE & clade == "L1",]

RepStatsDF <- as.data.frame(table(validCountsMatrixRepsDT$`class`), stringsAsFactors = FALSE) %>% mutate(percentage = round(Freq/sum(Freq)*100, 2))
LINECladeStatsDF <- as.data.frame(table(validCountsMatrixRepsDT[`class` == "LINE"]$`clade`), stringsAsFactors = FALSE) %>% mutate(percentage = round(Freq/sum(Freq)*100, 2))
LINE1SubFamStatsDF <- as.data.frame(table(validCountsMatrixRepsDT[`class` == "LINE" & clade == "L1"]$`consensus`), stringsAsFactors = FALSE) %>% mutate(percentage = round(Freq/sum(Freq)*100, 2))

#Function for repeats distributions
plotFunc <- function(dataIn, titleIn, xlabIn, ylabIn, angles=45, textAngle = 45, vjustIn = 1){
    plt <- ggplot(data = dataIn, aes(x = reorder(Var1, -percentage), y = percentage)) + 
        geom_bar(stat = "identity", fill = "steelblue") +
        geom_text(aes(label = Freq), vjust = vjustIn, size = 4, angle = textAngle) + 
        ylim(0, 100) +
        theme_minimal() +
        labs(x = "Repeat Element") +
        theme(axis.text.x = element_text(angle = angles, hjust = 1))
}

pltReps <- plotFunc(RepStatsDF, ylabIn = "Percentage")
pltLineClade <- plotFunc(LINECladeStatsDF, ylabIn = "Percentage")
pltL1 <- plotFunc(LINE1SubFamStatsDF, ylabIn = "Percentage",  angles=90, textAngle=90, vjustIn =0.5)


message(paste0("\n\nrunning locusSpecific? ", arguments$LocusSpecific, "\n\n"))
#create folder for locus-specific lengths
dir.create("LocusSpecific/figures", recursive = TRUE)
dir.create("LocusSpecific/data", recursive = TRUE)

ggsave(pltReps, filename = "LocusSpecific/figures/RepStatsDF.png", width = 10, height = 5)
ggsave(pltLineClade, filename = "LocusSpecific/figures/LINEClades.png", width = 10, height = 5)
ggsave(pltL1, filename = "LocusSpecific/figures/LINE1SubFamilies.png", width = 18, height = 5)

# names(validCountsMatrixL1DT) <- gsub("-", ".", names(validCountsMatrixL1DT))



#read squire data
counts_annot <- read_tsv(paste0(counts_annot_path))
counts_annot_pCoding <- counts_annot[counts_annot$`gene.type`=="protein_coding", ]
counts_annot_pCoding <-  counts_annot_pCoding %>% mutate_if(is.double, ~as.integer(.))

message("reading the metadata")
metadata_df <- read.csv(file = metadata_path, sep="," ,header = TRUE)

message(paste0("filtering  protein coding genes with", MIN_ReadsCounts, " reads"))  
counts_annot_pCoding <- counts_annot_pCoding %>% mutate(sumVar = select(., all_of(sampleNames)) %>% rowSums()) %>% filter(sumVar > MIN_ReadsCounts)


dimsL1=dim(validCountsMatrixL1DT)
message(paste0("there are ",dimsL1[1], " SquireL1"))


#remove rows with less than 10 reads
counts_Sq_ints_dt <- rbind(validCountsMatrixRepsDT,counts_annot_pCoding, fill=TRUE)
# counts_Sq_ints_dt[,RNAelementID := `gene.id`]
counts_Sq_ints_dt[, RNAelementID := fcase(`gene.type` == "protein_coding", `gene.id`, default = RNAelementID)] 

message(paste0("length of final Reps & Protein coding count matrix ",dim(counts_Sq_ints_dt)[1]))
# counts_Sq_ints_dt[is.na(`gene.type`),]
# table(counts_Sq_ints_dt$`gene.type`)

counts_Sq_pCoding_df <- counts_Sq_ints_dt %>% relocate(all_of(sampleNames), .after = last_col())

#need this for full length L1
fwrite(counts_Sq_pCoding_df, file = paste0("LocusSpecific/data/CountMatrix_squireRepeats_AggregatedFWDnREV_ProteinCoding.bed"), sep="\t", col.names = TRUE)

##################################################################
############ DE; all locus specific repeats. ################
##################################################################

#step 1: convert to dataframe
#step 2: metadata
#step 3: DE analysis

counts_Sq_pCoding_df_df <- as.data.frame(counts_Sq_pCoding_df)
rownames(counts_Sq_pCoding_df_df) <- counts_Sq_pCoding_df_df$RNAelementID
message("making countMatrix Table")

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
    x = 'log2FoldChange',
    y = 'padj',
    xlab = bquote(~Log[2]~ "FC(" * .(ALT) * "/" * .(REF) * ")"),
    pCutoff = 0.05,
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


#####DE all
# metadata_Copydf <- metadata_df
# metadata_Copydf$condition <- factor(metadata_Copydf$condition)
# ddsAll <- DESeqDataSetFromMatrix(countData = data.matrix(counts_Sq_pCoding_df_df[,sampleNames]),
#                               colData = metadata_Copydf,
#                               design = ~ condition)

# rowData(ddsAll) <- counts_Sq_pCoding_df_df[, c("gene.type","length","gene.symbol","ensembl.gene.id", "gene.id", "name","consensus", "clade")]

# ddsAll$condition <- relevel(ddsAll$condition, ref = "DMSO")
# keep <- rowSums( counts(ddsAll) >= 10 ) >= 3
# ddsAll <- ddsAll[keep,]
# ddsAll <- DESeq(ddsAll)



##################################################
#### DE; function to run loop ####################
##################################################
Repeat_DE <- function(data_in, metadata_in, REF, ALT){
cond_collapsed <- paste0(c("condition",ALT,REF), collapse="_")
#get metadata
metadata_short <- metadata_in %>% filter(condition %in% c(REF, ALT))
metadata_short <- metadata_short %>% mutate(condition = factor(condition, levels = c(REF, ALT))) #factor for DESeq2

ddsInFunc <- DESeqDataSetFromMatrix(countData = data.matrix(data_in[,metadata_short$samples]),
                              colData = metadata_short,
                              design = ~ condition)

rowData(ddsInFunc) <- data_in[, c("gene.type","length","gene.symbol","ensembl.gene.id", "gene.id", "name","consensus", "clade")]

keep <- rowSums( counts(ddsInFunc) >= 10 ) >= 3
ddsInFunc <- ddsInFunc[keep,]

#normalise repeats with non-repeats
dds_main <- DESeq(ddsInFunc)

### filter out repeats of interest
# dds_L1 <- dds_main[str_detect(rownames(dds_main),":L1:"),]
dds_L1 <- dds_main[!is.na(rowData(dds_main)$clade) & rowData(dds_main)$clade == "L1", ]
#not protein coding, not L1 is all repeats
dds_Repeats <- dds_main[!is.na(rowData(dds_main)$`gene.type`) & rowData(dds_main)$`gene.type` != "protein_coding" & rowData(dds_main)$`gene.type` != "L1", ] 

#just protein coding for ssanity checks
dds_pCoding <- dds_main[!is.na(rowData(dds_main)$`gene.type`) & rowData(dds_main)$`gene.type` == "protein_coding",]
dds_ERVs <- dds_main[!is.na(rowData(dds_main)$clade) & str_detect(rowData(dds_main)$clade, "ERV"),]

# table(rowData(dds_main)$`gene.type`)

ddsL1Results <- results(dds_L1, contrast=c("condition",ALT, REF))
ddsERVsResults <- results(dds_ERVs, contrast=c("condition",ALT, REF))
ddsRepeatsResults <- results(dds_Repeats, contrast=c("condition",ALT, REF))
ddsProteinCodingResults <- results(dds_pCoding, contrast=c("condition",ALT, REF))

ddsL1Results$sig <- ifelse(!is.na(ddsL1Results$padj) & abs(ddsL1Results$log2FoldChange) > 1 & is.na(ddsL1Results$padj) < 0.05, "sig", "notSig")
ddsL1Results$UpDown <- ifelse(ddsL1Results$sig == "sig" & ddsL1Results$log2FoldChange > 0, "up", ifelse(ddsL1Results$sig == "sig" & ddsL1Results$log2FoldChange < 0, "down", "notSig"))
proportions(table(ddsL1Results$UpDown))

# equalsRownames <- equals(rownames(rowData(dds_Repeats)), rownames(ddsRepeatsResults))
dt <- data.frame(ddsObj=equals(rownames(rowData(dds_Repeats)), rownames(ddsRepeatsResults)))
# dt[ddsObj==FALSE]#
# dt[!equalsRownames]
sigFilter <- which(!is.na(ddsRepeatsResults$padj) & abs(ddsRepeatsResults$log2FoldChange) > 1 & is.na(ddsRepeatsResults$padj) < 0.05)
ddsRepeatsResults$sig <- ifelse(!is.na(ddsRepeatsResults$padj) & abs(ddsRepeatsResults$log2FoldChange) > 1 & is.na(ddsRepeatsResults$padj) < 0.05, "sig", "notSig")
# ddsRepeatsResults[sigFilter, ]
ddsRepeatsSig <- ddsRepeatsResults[sigFilter, ]

#summary of DE results
ddsL1ResultsDf <- as.data.frame(ddsL1Results)
ddsRepeatsResultsDf <- as.data.frame(ddsRepeatsResults)


# Transform counts for data visualization
rld <- rlog(dds_Repeats, blind=TRUE)
# Plot PCA 
pltPCA_samples <- plotPCA(rld, intgroup=c("condition", "samples"))
pltPCA_ConditionsOnly <- plotPCA(rld, intgroup=c("condition"))
ggsave(pltPCA_samples, file = paste0("LocusSpecific/figures/PCA_by_condition_and_samples_",paste0(collapse="_"),cond_collapsed,".pdf"), width = 12, height = 12)
ggsave(pltPCA_ConditionsOnly, file = paste0("LocusSpecific/figures/PCA_by_conditionOnly_",paste0(collapse="_"),cond_collapsed,".pdf"), width = 12, height = 12)

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
rld_df <- rld_mat %>% as.data.frame() %>% rownames_to_column("RNAelementID")

#save a copy of the rlog matrix
fwrite(rld_df, file=paste0("LocusSpecific/data/DESeq2_rlog_Transform_Blind",cond_collapsed,".tsv"), sep="\t")


message("summary of locus-Spec L1 DE results")
message(summary(ddsL1Results))
message("/n")

#shrink results for volcano
resultsL1Shrink <- lfcShrink(dds_L1, contrast=c("condition",ALT, REF), res=ddsL1Results, type = 'normal')
resultsRepeatsShrink <- lfcShrink(dds_Repeats, contrast=c("condition",ALT, REF), res=ddsRepeatsResults, type = 'normal')
# resultsProteinCodingShrink <- lfcShrink(dds_pCoding, contrast = c("condition", REF, ALT), res=ddsProteinCodingResults, type = 'normal')

# row.names(ddsRepeatsResults)

# str_split_i(row.names(ddsRepeatsResults)[1:2], "\\|", 4)

message("making volcano plots")
volL1NoShrink <- plotEVolcano(ddsL1Results, title_in=cond_collapsed, subtitle_in="Differential Locus Specific L1 Expression", REF=REF, ALT=ALT, labs=str_split_i(row.names(ddsL1Results), "\\|", 4))
volL1Shrink <- plotEVolcano(resultsL1Shrink, title_in=cond_collapsed, subtitle_in="Differential Locus Specific L1 Expression", REF=REF, ALT=ALT, labs=str_split_i(row.names(resultsL1Shrink), "\\|", 4))

volERVsNoShrink <- plotEVolcano(ddsERVsResults, title_in=cond_collapsed, subtitle_in="Differential ERVs Expression", REF=REF, ALT=ALT, labs=str_split_i(row.names(ddsERVsResults), "\\|", 4))


volRepeatsShrink <- plotEVolcano(resultsRepeatsShrink, title_in=cond_collapsed, subtitle_in="Differential Locus Specific Repeats Expression", REF=REF, ALT=ALT, labs=str_split_i(row.names(ddsRepeatsResults), "\\|", 4))
volRepeatsNoShrink <- plotEVolcano(ddsRepeatsResults, title_in=cond_collapsed, subtitle_in="Differential Locus Specific Repeats Expression", REF=REF, ALT=ALT, labs=str_split_i(row.names(ddsRepeatsResults), "\\|", 4))

message("save volcano plots")
ggsave(volL1Shrink, filename = paste0("LocusSpecific/figures/volcanoShrink_locusSpecificLINE1_",cond_collapsed,".svg"))
ggsave(volL1NoShrink, filename = paste0("LocusSpecific/figures/volcanoNoShrink_locusSpecificLINE1_",cond_collapsed,".svg"))
ggsave(volL1NoShrink, filename = paste0("LocusSpecific/figures/volcanoNoShrink_locusSpecificLINE1_",cond_collapsed,".pdf"))

ggsave(volERVsNoShrink, filename = paste0("LocusSpecific/figures/volcanoNoShrink_locusERVs_",cond_collapsed,".pdf"))
ggsave(volRepeatsShrink, filename = paste0("LocusSpecific/figures/volcanoShrink_locusSpecificRepeats_",cond_collapsed,".svg"))
ggsave(volRepeatsNoShrink, filename = paste0("LocusSpecific/figures/volcanoNoShrink_locusSpecificRepeats_",cond_collapsed,".svg"))
ggsave(volRepeatsNoShrink, filename = paste0("LocusSpecific/figures/volcanoNoShrink_locusSpecificRepeats_",cond_collapsed,".pdf"))
# ggsave(plotEVolcano(resultsProteinCodingShrink, title_in=cond_collapsed), filename = "volcano_proteinCoding_",cond_collapsed,".png")

#make CPM
message("making CPM")
dds_locusSpecL1_cpm_df <- makeCPM(DeSeq2_res=dds_L1, prior_cnts=2, keyColName="te.id")
dds_repeats_cpm_df <- makeCPM(DeSeq2_res=dds_Repeats, prior_cnts=2, keyColName="te.id")

#save results
message("saving results to disk")
results2SaveLspecL1 <- list(ddsMain=dds_main, ddsFiltered=dds_L1, resDF=ddsL1Results, resShrinkedDF=resultsL1Shrink, cpmDF=dds_locusSpecL1_cpm_df, volcanoShrink_ggplot = volL1Shrink, volcanoNoShrink_ggplot = volL1NoShrink, metaData = metadata_short, normalizedCounts=rld_df)
results2SaveLspecRepeats <- list(ddsMain=dds_main, ddsFiltered=dds_Repeats, resDF=ddsRepeatsResultsDf, resShrinkedDF=resultsRepeatsShrink, cpmDF=dds_repeats_cpm_df, volcano_ggplot = volRepeatsShrink, metaData = metadata_short)
save(results2SaveLspecL1, file=paste0("LocusSpecific/data/Differential_locusSpecific_L1ExpressionResults_",cond_collapsed,".RData"))  # Till here everything is ok. 
save(results2SaveLspecRepeats, file=paste0("LocusSpecific/data/Differential_locusSpecific_RepeatsExpressionResults_",cond_collapsed,".RData"))
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

message("done with locuspecific line1 DESeq2 analysis")

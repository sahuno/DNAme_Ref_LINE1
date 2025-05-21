#here, scripts locus specfic line1 Analysis
#input: total counts table, metadata, 
#outputs: volcano plots, rdata DEResults, pca of repeat 



#author: Samuel ahuno
#date: April 3rd 2025
#purpose: differential expression analysis of locus specific L1
#1: merge readcounts from l1 and protein coding genes: out_L1_and_proteinCoding_genes_readcounts.csv
#2. use `out_L1_and_proteinCoding_genes_readcounts.csv` as input for active L1, locus specific L1


# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/DE_locusSpecificLINE1.R --LocusSpecific FALSE

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
metadata_path <- "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_LowQualSamplesDropped_recoded.csv"
path_map_locusL1_activeFullLen <- arguments$pathMapL1BaseRepMasker

sampleNames <- c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3","R.C.1", "R.C.2", "R.Q.1", "R.Q.2", "R.Q.3","R.QC.1", "R.QC.2", "R.QC.3", "R.S.1", "R.S.3","R.SC.1", "R.SC.2", "R.SC.3")
RefAlt_samples <- c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3")
DE_conditions <- data.frame(condition = c("DMSO;AZA", "DMSO;CKi", "DMSO;QSTAT", "DMSO;QSTAT-CKi"))
# threads_fst(set_fst_threads)

sampleNames_focus_Qstats <- c("R.0.1", "R.0.2", "R.0.3", "R.C.1", "R.C.2", "R.Q.1", "R.Q.2", "R.Q.3","R.QC.1", "R.QC.2", "R.QC.3")


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

# head(pltL1)
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
counts_Sq_ints_dt <- rbind(validCountsMatrixRepsDT, counts_annot_pCoding, fill=TRUE)
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

ddsL1Results <- results(dds_L1, contrast=c("condition", ALT, REF))
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
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################




# head(counts_Sq_pCoding_df_df)
# tail(counts_Sq_pCoding_df_df)

message("creating normalized counts all comditions")
# important_types <- c("SINE",  "LTR", "LINE", "DNA", "protein_coding")
important_types <- c("SINE",  "LTR", "LINE", "DNA")

df_proteinCoding_repeats <- counts_Sq_pCoding_df_df %>% filter(`gene.type` %in% important_types)


mat_df <- df_proteinCoding_repeats[, sampleNames_focus_Qstats]
rownames(head(df_proteinCoding_repeats))

table(df_proteinCoding_repeats$gene.type)
table(counts_Sq_pCoding_df_df$gene.type)



metadata_short <- metadata_df %>% filter(condition %in% c("DMSO", "QSTAT", "CKi", "QSTAT-CKi"))
metadata_short <- metadata_short %>% mutate(condition = factor(condition, levels = c("DMSO", "QSTAT", "CKi", "QSTAT-CKi"))) #factor for DESeq2

ddsInFunc <- DESeqDataSetFromMatrix(countData = data.matrix(mat_df),
                              colData = metadata_short,
                              design = ~ condition)

rowData(ddsInFunc) <- df_proteinCoding_repeats[, c("gene.type","length","gene.symbol","ensembl.gene.id", "gene.id", "name","consensus", "clade")]
head(counts(ddsInFunc))

pdf("hist_reads.pdf")
hist(log10(rowSums(counts(ddsInFunc))))
dev.off()
# head(counts(ddsInFunc))

# rownames(counts(ddsInFunc))[(rownames(counts(ddsInFunc)) == "chr1|4558487|4558680|IAPEz-int:ERVK:LTR|10|-")]

# df_counts <- data.frame(counts(ddsInFunc)) %>% rownames_to_column("RNAelementID") %>% filter(RNAelementID == "chr1|4558487|4558680|IAPEz-int:ERVK:LTR|10|-")
# df_counts <- data.frame(mat_df) %>% rownames_to_column("RNAelementID") %>% filter(RNAelementID == "chr1|4558487|4558680|IAPEz-int:ERVK:LTR|10|-")

# rownames(counts(ddsInFunc))
# counts(ddsInFunc) 
keep <- rowSums(counts(ddsInFunc) >= 3 ) >= 10
ddsInFunc <- ddsInFunc[keep,]

#normalise repeats with non-repeats
dds_main <- DESeq(ddsInFunc)
# dim(dds_main)
# head(counts(dds_main))

# results(dds_main)

#not protein coding, not L1 is all repeats
dds_Repeats <- dds_main[!is.na(rowData(dds_main)$`gene.type`) & rowData(dds_main)$`gene.type` != "protein_coding" & rowData(dds_main)$`gene.type` != "L1", ] 
#just protein coding for ssanity checks
# dds_pCoding <- dds_main[!is.na(rowData(dds_main)$`gene.type`) & rowData(dds_main)$`gene.type` == "protein_coding",]
# dim(dds_pCoding)
#line1
dds_L1 <- dds_main[!is.na(rowData(dds_main)$clade) & rowData(dds_main)$clade == "L1", ]
dds_ERVs <- dds_main[!is.na(rowData(dds_main)$clade) & str_detect(rowData(dds_main)$clade, "ERV"),]



# dds_subsetType <- dds_main[!is.na(rowData(dds_main)$`gene.type`) & rowData(dds_main)$`gene.type` == uniq_genomicElements[1], ]

# # Transform counts for data visualization
# rld <- rlog(dds_main, blind=FALSE)
# # Plot PCA 
# pltPCA_samples <- plotPCA(rld, intgroup=c("condition", "samples"), ntop = Inf)
# pltPCA_ConditionsOnly <- plotPCA(rld, intgroup=c("condition"))
# ggsave(pltPCA_samples, file = paste0("LocusSpecific/figures/PCA_all_conditions_with_samples",paste0(collapse="_"),".pdf"), width = 12, height = 12)
# ggsave(pltPCA_ConditionsOnly, file = paste0("LocusSpecific/figures/PCA_all_conditions_colored",paste0(collapse="_"),".pdf"), width = 12, height = 12)

# # Extract the rlog matrix from the object and compute pairwise correlation values
# rld_mat <- assay(rld)
# rld_cor <- cor(rld_mat)
# rld_df <- rld_mat %>% as.data.frame() %>% rownames_to_column("RNAelementID")
# #save a copy of the rlog matrix
# fwrite(rld_df, file=paste0("LocusSpecific/data/DESeq2_rlog_Transform_Blind",cond_collapsed,".tsv"), sep="\t")



#### diagnotic plots for siyu
# at top of your script, make sure to load:
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tibble)
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(cola)

#-------------------------------------------------------------------------------
#' Run rlog, PCA plotting, export rlog matrix, plot sample similarity,
#' and adjust matrix for downstream activities
#'
#' @param dds              A DESeqDataSet (already run through DESeq())
#' @param cond_label       A short label for this dataset (used in file names)
#' @param groups           Character vector of colData columns for the full PCA
#' @param primary_group    Single colData column for the second PCA
#' @param output_prefix    Base directory under which "figures/" and "data/" live
#' @return                 Invisibly, a list with elements:
#'                         rld (RLog object),
#'                         matrix (rlog assay matrix),
#'                         correlation (matrix),
#'                         sampleDist (distance matrix),
#'                         adjustedMatrix (adjusted matrix),
#'                         adjustedMatrixT (transposed adjusted matrix)
#-------------------------------------------------------------------------------
run_deseq_viz <- function(dds,
                          cond_label,
                          groups = c("condition", "samples"),
                          primary_group = "condition",
                          output_prefix = "LocusSpecific") {
  # create output dirs
  fig_dir  <- file.path(output_prefix, "figures")
  data_dir <- file.path(output_prefix, "data")
  dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1. rlog transform
  rld <- rlog(dds, blind = FALSE)
  
  # 2. PCA plots (using default ntop)
  plt_full <- plotPCA(rld, intgroup = groups)
  plt_cond <- plotPCA(rld, intgroup = primary_group)
  
  file_full <- file.path(fig_dir,
                         paste0("PCA_", cond_label, "_",
                                paste(groups, collapse = "_"), ".pdf"))
  file_cond <- file.path(fig_dir,
                         paste0("PCA_", cond_label, "_by_",
                                primary_group, ".pdf"))
  
  ggsave(filename = file_full, plot = plt_full, width = 12, height = 12)
  ggsave(filename = file_cond, plot = plt_cond, width = 12, height = 12)
  
  # 3. extract rlog matrix & compute correlations
  rld_mat <- assay(rld)
  cor_mat <- cor(rld_mat)
  
  # 4. write out as TSV
  rld_df <- as.data.frame(rld_mat) %>%
            rownames_to_column("RNAelementID")
  out_tsv <- file.path(data_dir,
                       paste0("DESeq2_rlog_", cond_label, ".tsv"))
  fwrite(rld_df, file = out_tsv, sep = "\t")
  
  # 5. sample-distance heatmap
  samp_dist    <- dist(t(rld_mat))
  samp_dist_mx <- as.matrix(samp_dist)
  heat_file    <- file.path(fig_dir,
                            paste0("sampleDist_", cond_label, ".pdf"))
  
  pdf(heat_file, width = 10, height = 10)
  pheatmap(samp_dist_mx,
           clustering_distance_rows = samp_dist,
           clustering_distance_cols = samp_dist,
           color = colorRampPalette(
             rev(brewer.pal(9, "Blues")))(255)
  )
  dev.off()
  
  # 6. adjust matrix for heatmap & downstream use
  heatmap_matrixAdjusted <- cola::adjust_matrix(
    rld_mat,
    sd_quantile = 0.05,
    max_na = 0.25,
    verbose = TRUE
  )
  heatmap_matrixAdjustedT <- t(heatmap_matrixAdjusted)
  
  # return invisibly for downstream use
  invisible(list(
    rld                = rld,
    matrix             = rld_mat,
    correlation        = cor_mat,
    sampleDist         = samp_dist_mx,
    adjustedMatrix     = heatmap_matrixAdjusted,
    adjustedMatrixT    = heatmap_matrixAdjustedT
  ))
}




# head(rld_mat)


# now call the function on each
res_repeats  <- run_deseq_viz(dds_Repeats, cond_label = "Repeats")
# res_pCoding  <- run_deseq_viz(dds_pCoding, cond_label = "ProteinCoding")
res_L1  <- run_deseq_viz(dds_L1, cond_label = "LINE1")
res_ERV  <- run_deseq_viz(dds_ERVs, cond_label = "ERVs")
# names(res_repeats)
# head(res_repeats$matrix)


# sampleDists <- dist(t(res_repeats$matrix))
# sampleDistMatrix <- as.matrix(sampleDists)
# head(sampleDistMatrix)
library("pheatmap")
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# pdf("locusSpec_sampleDistMatrix.pdf")
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)
# dev.off()



# table(counts_Sq_pCoding_df_df$gene.type)
# uniq_genomicElements <- unique(counts_Sq_pCoding_df_df$gene.type)
all_types <- c("SINE",  "LTR","LINE", "DNA")

types_to_do <- all_types
dds_subsetTypes <- setNames(
  lapply(types_to_do, function(gt) {
    dds_main[rowData(dds_main)$gene.type == gt, ]
  }),
  types_to_do
)
names(dds_subsetTypes)

# 1. run the viz function on each subset
viz_results_by_type <- lapply(names(dds_subsetTypes), function(gt) {
  message("Processing gene.type: ", gt)
  run_deseq_viz(
    dds           = dds_subsetTypes[[gt]],
    cond_label    = gt,
    groups        = c("condition", "samples"),
    primary_group = "condition",
    output_prefix = "LocusSpecific"
  )
})

# 2. give it back the names for easy lookup
names(viz_results_by_type) <- names(dds_subsetTypes)


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

message("done with locuspecific line1 DESeq2 analysis")




##############################################################################################################################
##############################################clustering##########################################################
##############################################################################################################################
# metadata_short

# subSampleCondsDF <-  data.frame(conditions = subSampleConds)
# DrugsAnnot <- subSampleCondsDF %>% dplyr::filter(conditions %in% subSampleConds) %>% 

library(dplyr)
library(stringr)

DrugsAnnot <- metadata_short %>%
  # 1) replace hyphens with underscores
  mutate(condition = gsub("-", "_", condition)) %>%
  # 2) derive new columns
  mutate(
    # combo vs mono
    Combinations = if_else(str_detect(condition, "_"), "combo", "mono"),
    # BroadTargets: check the _combo_ patterns first
    BroadTargets = case_when(
      str_detect(condition, "QSTAT_CKi")   ~ "Chromatin+MEK",
      str_detect(condition, "QSTAT")       ~ "Chromatin",
      str_detect(condition, "SETDB1i_CKi") ~ "H3K9me+MEK",
      str_detect(condition, "SETDB1i_")    ~ "H3K9me",
      str_detect(condition, "AZA")         ~ "DNAme",
      str_detect(condition, "DMSO")         ~ "DMSO",
      str_detect(condition, "CKi")         ~ "MEK",
      TRUE                                 ~ NA_character_
    ),
    # Action: again, specific combos before general
    Action = case_when(
      str_detect(condition, "QSTAT_CKi")   ~ "HDACi+MEKi",
      str_detect(condition, "SETDB1i_CKi") ~ "H3K9mei+MEKi",
      str_detect(condition, "QSTAT")       ~ "HDACi",
      str_detect(condition, "SETDB1i_")    ~ "H3K9mei",
      str_detect(condition, "AZA")         ~ "DNMTi",
      str_detect(condition, "DMSO")         ~ "DMSO",
      str_detect(condition, "CKi")    ~ "MEKi",
      TRUE                                 ~ NA_character_
    )
  )



rownames(DrugsAnnot) <- DrugsAnnot$samples



#sanity checks
# names(res_repeats)
# heatmap_matrixAdjusted <- cola::adjust_matrix(head(res_repeats$matrix, 10), sd_quantile = 0.05, max_na = 0.25, verbose = TRUE)
# heatmap_matrixAdjusted_t = t(heatmap_matrixAdjusted)
# head(heatmap_matrixAdjusted_t)

if(!require('NbClust')) {
  install.packages('NbClust')
  library('NbClust')
}




## Function to run k-means clustering, visualize results, and summarize class proportions for any k
## Function to run k-means clustering, visualize results, and summarize class proportions for any k
cluster_heatmap_kmeans <- function(
  data_matrix,
  k,
  annotation_df,
  output_dir = ".",
  prefix = NULL,
  nstart = 25,
  seed = 983
) {
  # Load required packages
  require(ggplot2)
  require(factoextra)
  require(pheatmap)
  require(dplyr)
  require(tidyr)
  require(readr)
  require(stringr)
  require(tibble)
  require(scales)

  # Prepare output directory and filename prefix
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(prefix)) {
    prefix <- paste0("kmeans", k)
  }

  # Ensuring annotation_df has proper rownames for pheatmap
  annotation_df <- as.data.frame(annotation_df)
  if ("samples" %in% colnames(annotation_df)) {
    rownames(annotation_df) <- annotation_df$samples
    annotation_df$samples <- NULL
  }
  # Check matching
  if (!all(colnames(data_matrix) %in% rownames(annotation_df))) {
    stop("Row names of annotation_df must match column names of data_matrix.")
  }
  # Reorder annotation_df to match data_matrix
  annotation_df <- annotation_df[colnames(data_matrix), , drop = FALSE]

  # Run k-means
  set.seed(seed)
  km_res <- kmeans(data_matrix, centers = k, nstart = nstart)

  # Build clusterDF with RepID preserved
  clusterDF <- tibble(
    RepID   = names(km_res$cluster),
    Cluster = factor(km_res$cluster)
  ) %>%
    separate(
      col    = RepID,
      into   = c("chr", "start", "end", "RepName"),
      sep    = "\\|",
      remove = FALSE,
      extra  = "drop",
      fill   = "right"
    ) %>%
    mutate(
      GenomicCoords = paste0(chr, ":", start, "-", end),
      ClassFamily   = str_split(RepName, ":", simplify = TRUE)[,1]
    ) %>%
    separate(RepName, into = c("subfamily", "family", "class"), sep = ":")

  # 1) Scatter plot of clusters
  scatter_plot <- fviz_cluster(
    km_res,
    data         = data_matrix,
    geom         = "point",
    ellipse.type = "convex",
    palette      = "jco",
    ggtheme      = theme_minimal()
  )
  ggsave(
    filename = file.path(output_dir, paste0(prefix, "_clusters.png")),
    plot     = scatter_plot,
    width    = 6,
    height   = 4,
    dpi      = 300
  )

  # 2) Heatmap with annotations
  annot_row <- data.frame(Cluster = clusterDF$Cluster)
  rownames(annot_row) <- clusterDF$RepID
  pheatmap(
    mat            = data_matrix,
    annotation_row = annot_row,
    annotation_col = annotation_df,
    scale          = "row",
    show_rownames  = FALSE,
    cluster_rows   = TRUE,
    main           = paste0("k = ", k, " clusters"),
    filename       = file.path(output_dir, paste0(prefix, "_heatmap.pdf"))
  )

  # 3) Calculate class proportions
  prop_df <- clusterDF %>%
    group_by(Cluster, class) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Cluster) %>%
    mutate(prop = count / sum(count))

  # 4) Subset to major repeat classes
  major_classes <- c("SINE", "LTR", "LINE", "DNA")
  prop_major <- prop_df %>%
    filter(class %in% major_classes)

  # 5) Plot class proportions
  class_plot <- ggplot(prop_df, aes(x = Cluster, y = prop, fill = class)) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = percent) +
    labs(
      title = paste0("Repeat class proportions (k=", k, ")"),
      x     = "Cluster",
      y     = "Proportion",
      fill  = "Class"
    ) +
    theme_minimal()
  ggsave(
    filename = file.path(output_dir, paste0(prefix, "_class_prop.pdf")),
    plot     = class_plot,
    width    = 7,
    height   = 5
  )

  # 6) Plot major class proportions
  major_plot <- ggplot(prop_major, aes(x = Cluster, y = prop, fill = class)) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = percent) +
    labs(
      title = paste0("Major repeat class proportions (k=", k, ")"),
      x     = "Cluster",
      y     = "Proportion",
      fill  = "Class"
    ) +
    theme_minimal()
  ggsave(
    filename = file.path(output_dir, paste0(prefix, "_major_class_prop.pdf")),
    plot     = major_plot,
    width    = 7,
    height   = 5
  )

  # 7) Save results and return
  write_tsv(prop_df,    file = file.path(output_dir, paste0(prefix, "_prop.tsv")))
  write_tsv(prop_major, file = file.path(output_dir, paste0(prefix, "_major_prop.tsv")))
  save(
    data_matrix, clusterDF, annotation_df, k, km_res,
    file = file.path(output_dir, paste0(prefix, "_results.RData"))
  )

  invisible(list(
    kmeans    = km_res,
    clusterDF = clusterDF,
    plots     = list(
      scatter     = scatter_plot,
      class_prop  = class_plot,
      major_prop  = major_plot
    )
  ))
}

# [,-c(2,3) ]
message("running kmeans with k=3 clustering")
k3_clus_res <- cluster_heatmap_kmeans(data_matrix=res_repeats$adjustedMatrix, k=3, annotation_df = DrugsAnnot %>% select(-c(new_samples_name, condition_long, samples)), output_dir = "LocusSpecific/REP_kmeans_3", prefix = "REP_kmeans_3", nstart = 25, seed = 983)

# message("running kmeans clustering with k=OptimalK ")
# k3_clus_res <- cluster_heatmap_kmeans(data_matrix=res_repeats$adjustedMatrix, k=3, annotation_df = DrugsAnnot %>% select(-c(new_samples_name, condition_long, samples)), output_dir = "LocusSpecific/kmeans_3", prefix = "kmeans_3", nstart = 25, seed = 983)

message("running kmeans with k=2 clustering")
k2_clus_res <- cluster_heatmap_kmeans(data_matrix=res_repeats$adjustedMatrix, k=2, annotation_df = DrugsAnnot %>% select(-c(new_samples_name, condition_long, samples)), output_dir = "LocusSpecific/REP_kmeans_2", prefix = "REP_kmeans_2", nstart = 25, seed = 983)

# head(res_repeats$adjustedMatrix)
head(k3_clus_res$clusterDF)

message("running NbClust to find optimal number of clusters\n")
# what is the difference between adjustedMatrix(genes by samples)=group genes and adjustedMatrixT(samples by genes)=group samples?
# set.seed(123)  # seed for reproducibility
# clusterNum <- NbClust(res_repeats$adjustedMatrix, distance = "euclidean", min.nc = 2, max.nc = 12, method = "kmeans", index = "silhouette")
# message("optimal number of clusters: ", clusterNum$Best.nc[1], "; with value index: ", clusterNum$Best.nc[2])

# head(clusterNum$Best.partition)
# head(clusterNum$Best.nc)
# head(clusterNum$All.index)
# table(clusterNum$Best.nc[1])

##############################################################################################################################
########################################## End of Clustering ##################################################
##############################################################################################################################





##############################################################################################################################
########################################## Beginining of Annotations ##################################################
##############################################################################################################################

library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rtracklayer)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)  # GRanges of all known genes


annotate_repeat_clusters <- function(
  cluster_df,
  txdb,
  orgdb,
  annoDb       = "org.Mm.eg.db",
  class_filter = c("SINE","LTR","LINE","DNA"),
  tssRegion    = c(-3000, 3000),
  out_dir      = "repeat_cluster_results",
  prefix       = "cluster"
){
  # ─── Dependencies ───────────────────────────────────────────────────────────
  require(GenomicRanges)
  require(ChIPseeker)
  require(AnnotationDbi)
  require(dplyr)
  require(ggplot2)
  
  # ─── Prepare output folder ──────────────────────────────────────────────────
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
  
  # ─── 1) Clean types ─────────────────────────────────────────────────────────
  df <- cluster_df %>%
    mutate(
      chr   = as.character(chr),
      start = as.integer(start),
      end   = as.integer(end)
    )
  
  # ─── 2) Build GRanges ───────────────────────────────────────────────────────
  gr <- makeGRangesFromDataFrame(df,
                                 seqnames.field     = "chr",
                                 start.field        = "start",
                                 end.field          = "end",
                                 keep.extra.columns = TRUE)
  genes_gr <- GenomicFeatures::genes(txdb)
  
  # ─── 3) Intragenic vs Intergenic ────────────────────────────────────────────
  df$GenicStatus <- "intergenic"
  hits <- findOverlaps(gr, genes_gr)
  df$GenicStatus[queryHits(hits)] <- "intragenic"
  
  # ─── 4) Nearest gene & distance ─────────────────────────────────────────────
  nearest <- distanceToNearest(gr, genes_gr)
  matched <- queryHits(nearest)
  gene_ids <- genes_gr$gene_id[subjectHits(nearest)]
  syms     <- mapIds(orgdb,
                     keys     = gene_ids,
                     column   = "SYMBOL",
                     keytype  = "ENTREZID",
                     multiVals= "first")
  df$NearestGene    <- NA_character_
  df$DistanceToGene <- NA_integer_
  df$NearestGene[matched]    <- syms
  df$DistanceToGene[matched] <- mcols(nearest)$distance
  
  # ─── 5) Promoter/Intron/Exon annotation via ChIPseeker ──────────────────────
  peak_anno <- annotatePeak(
    gr,
    TxDb      = txdb,
    annoDb    = annoDb,
    tssRegion = tssRegion
  )
  anno_df <- as.data.frame(peak_anno) %>%
    mutate(
      GenomicCoords = paste0(seqnames, ":", start, "-", end),
      SYMBOL_unique = mapIds(orgdb,
                             keys      = geneId,
                             column    = "SYMBOL",
                             keytype   = "ENTREZID",
                             multiVals = "first")
    )
  
  df <- df %>%
    mutate(GenomicCoords = paste0(chr, ":", start, "-", end)) %>%
    left_join(
      anno_df %>%
        dplyr::select(GenomicCoords, annotation, geneId, SYMBOL_unique, distanceToTSS),
      by = "GenomicCoords"
    )
  
  # ─── 6) Compute proportions ─────────────────────────────────────────────────
  print(head(df))
  prop_df <- df %>%
    filter(class %in% class_filter) %>%
    group_by(Cluster, class, GenicStatus) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Cluster, class) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
  # ─── 7) Build & save bar‐chart ───────────────────────────────────────────────
  plt <- ggplot(prop_df, aes(x = class, y = proportion, fill = GenicStatus)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ Cluster, nrow = 1) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x     = "Repeat Class",
      y     = "Proportion",
      fill  = "Genic Status",
      title = "Intergenic vs Intragenic Proportions by Class"
    ) +
    theme_minimal() +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    )
  
  out_file <- file.path(out_dir, paste0(prefix, "_genic_proportions.pdf"))
  ggsave(filename = out_file, plot = plt, width = 10, height = 6)
  message("→ Saved proportions plot to: ", out_file)
  
  # ─── 8) Compute family‐level proportions and plot ───────────────────────────
  ## <=== Highlighted section: family‐level proportions and barplot ===
  family_df <- df %>%
    group_by(Cluster, family, GenicStatus) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Cluster, family) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()

  family_plot <- ggplot(
    family_df,
    aes(x = family, y = proportion, fill = GenicStatus)
  ) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Cluster, scales = "free_x") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x     = "Repeat Family",
      y     = "Proportion",
      fill  = "Genic Status",
      title = "Intergenic vs Intragenic Proportions by Family"
    ) +
    theme_minimal() +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    )
  fam_plot_file <- file.path(out_dir, paste0(prefix, "_family_proportions.pdf"))
  ggsave(fam_plot_file, plot = family_plot, width = 12, height = 6)
  message("→ Saved family proportions plot to: ", fam_plot_file)


   # ─── 7) Save results ────────────────────────────────────────────────────────
  # Annotated repeats
  annotated_file <- file.path(out_dir, paste0(prefix, "_annotated_repeats.csv"))
  write_csv(df, annotated_file)
  message("→ Saved annotated repeats to: ", annotated_file)

  # Proportions table
  prop_file <- file.path(out_dir, paste0(prefix, "_genic_proportions.csv"))
  write_csv(prop_df, prop_file)
  message("→ Saved proportions table to: ", prop_file)

  message("done annotating repeat clusters")
  # ─── 8) Return results ───────────────────────────────────────────────────────
# ─── 9) Plot distanceToTSS histogram by cluster ──────────────────────────────
  hist_plot <- ggplot(df, aes(x = distanceToTSS)) +
    geom_histogram(bins = 100) +
    facet_wrap(~ Cluster, scales = "free") +
    labs(x = "Distance to TSS (bp)", y = "Count",
         title = "Distance to TSS Distribution by Cluster") +
    theme_minimal()
  hist_file <- file.path(out_dir, paste0(prefix, "_distanceToTSS_histograms.pdf"))
  ggsave(hist_file, plot = hist_plot, width = 12, height = 6)
  message("→ Saved distanceToTSS histograms to: ", hist_file)

  # ─── 10) Return results ─────────────────────────────────────────────────────
  invisible(list(
    annotated_df   = df,
    prop_df        = prop_df,
    class_plot     = plt,
    family_prop    = family_df,
    family_plot    = family_plot,
    distance_hist  = hist_plot
  ))
}




##############################################################################
##########################.  run enrichemnt analysis   ##########################
##############################################################################
run_enrichment_pipeline <- function(
  annotated_df,
  orgdb,
  out_dir = "enrichment_results",
  prefix  = "cluster",
  categories2Show = 10
) {
  # ─── Dependencies ───────────────────────────────────────────────────────────
  require(clusterProfiler)
  require(DOSE)
  require(AnnotationDbi)
  require(dplyr)
  require(ggplot2)
  require(readr)
  
  # ─── Prepare output folder ──────────────────────────────────────────────────
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # ─── 1. Prepare gene lists per cluster ─────────────────────────────────────
  cluster_genes <- annotated_df %>%
    filter(!is.na(SYMBOL_unique)) %>%
    distinct(Cluster, SYMBOL_unique) %>%
    split(.$Cluster)
  names(cluster_genes) <- paste0(prefix, names(cluster_genes))
  
  # Extract symbol vectors
  cluster_symbols <- lapply(cluster_genes, function(df) df$SYMBOL_unique)
  
  message("Map SYMBOL to ENTREZID")
  # ─── 2. Map SYMBOL to ENTREZID ─────────────────────────────────────────────
  cluster_entrez <- lapply(cluster_symbols, function(symbols) {
    mapIds(orgdb,
          gene          = genes,
            # org.Mm.eg.db
          # annoDb       = "org.Mm.eg.db",
           keys     = symbols,
           column   = "ENTREZID",
           keytype  = "SYMBOL",
           multiVals= "first") %>%
      na.omit() %>% unique()
  })
  
  # ─── 3. Enrichment loops ────────────────────────────────────────────────────
  go_results <- list(); kegg_results <- list()
  for (cl in names(cluster_entrez)) {
    genes <- cluster_entrez[[cl]]
    go_results[[cl]] <- list(
      BP = enrichGO(genes, OrgDb=orgdb, keyType="ENTREZID", ont="BP",
                    pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.10),
      MF = enrichGO(genes, OrgDb=orgdb, keyType="ENTREZID", ont="MF", 
                    pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.10),
      CC = enrichGO(genes, OrgDb=orgdb, keyType="ENTREZID", ont="CC", 
                    pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.10)
    )
    kegg_results[[cl]] <- enrichKEGG(genes, organism="mmu", keyType="kegg",
                                     pAdjustMethod="BH", pvalueCutoff=0.05)
  }
  
  # ─── 4. Combine GO results ─────────────────────────────────────────────────
  combined_go <- bind_rows(lapply(names(go_results), function(cl) {
    bind_rows(lapply(c("BP","MF","CC"), function(ont) {
      res <- go_results[[cl]][[ont]]@result
      if (nrow(res)==0) return(tibble(Cluster=cl,Ontology=ont))
      df <- as_tibble(res); df$Cluster<-cl; df$Ontology<-ont; df
    }))
  }))
  write_csv(combined_go, file.path(out_dir, paste0(prefix, "_GO_combined.csv")))
  
  # ─── 5. Plot dotplots ──────────────────────────────────────────────────────
  for (ont in c("BP","MF","CC")) {
    pdf(file.path(out_dir, paste0(prefix, "_GO_", ont, ".pdf")), width=10, height=6)
    for (cl in names(go_results)) {
      res <- go_results[[cl]][[ont]]
      if (nrow(res@result)>0) print(dotplot(res, showCategory=categories2Show)+ ggtitle(paste(cl,ont)))
    }
    dev.off()
  }
  
  # ─── 6. KEGG dotplots ───────────────────────────────────────────────────────
  pdf(file.path(out_dir, paste0(prefix, "_KEGG.pdf")), width=10, height=6)
  for (cl in names(kegg_results)) {
    res <- kegg_results[[cl]]
    if (nrow(res@result)>0) print(dotplot(res, showCategory=categories2Show)+ ggtitle(cl))
  }
  dev.off()
  
  # ─── 7. Return lists ───────────────────────────────────────────────────────
  invisible(list(
    cluster_entrez = cluster_entrez,
    go_results     = go_results,
    kegg_results   = kegg_results,
    combined_go    = combined_go
  ))
}




res_clusterAnnotate <-  annotate_repeat_clusters(k3_clus_res$clusterDF, txdb   = TxDb.Mmusculus.UCSC.mm10.knownGene,
  orgdb  = org.Mm.eg.db,
  out_dir = "REP_Annot_results_k3",
  prefix  = "REP_k3")

# run_enrichment_pipeline(res_clusterAnnotate$annotated_df)
go_clusterAnnotate <- run_enrichment_pipeline(res_clusterAnnotate[["annotated_df"]], orgdb=org.Mm.eg.db, out_dir = "REP_enrichment_results_k3", prefix  = "REP_cluster_k3_")

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
library(UpSetR)

##QUE: do youy need different full length L1 datasets for this to run? No
##n
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/DE_FullLength.R --FullLengthL1 TRUE


##old way

####begin getting user inputs 
option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
        help="Print extra output [default]"),
    make_option(c("-q", "--quietly"), action="store_false",
        dest="verbose", help="Print little output"),
    make_option(c("-l", "--lfcCutoff"), type="numeric", default=0.5, help="Cutoff for adjusted p-value [default %default]"),
    make_option(c("-c", "--count"), type="integer", default=10,
        help="minimum RNA read counts [default %default]",
        metavar="number"),
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
MIN_ReadsCounts = arguments$count
smallestGroupSize <- 3
# set_fst_threads <- 4
counts_annot_path <- '/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/CT/counts_annot.tsv'
metadata_path <- "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_recoded.csv"
# path_map_locusL1_activeFullLen11 <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mapped_repeatMasker_AssignedTo_L1Base_mm10_mmflil1_8438.tsv"
path_map_locusL1_activeFullLen <- arguments$pathMapL1BaseRepMasker

sampleNames <- c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3","R.C.1", "R.C.2", "R.C.3", "R.Q.1", "R.Q.2", "R.Q.3","R.QC.1", "R.QC.2", "R.QC.3", "R.S.1", "R.S.2", "R.S.3","R.SC.1", "R.SC.2", "R.SC.3")
RefAlt_samples <- c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3")
DE_conditions <- data.frame(condition = c("DMSO;AZA", "DMSO;CKi", "DMSO;QSTAT", "DMSO;QSTAT-CKi"))

methodComputeMultipleLocusSpecL1ReadCounts = arguments$StatsMultipleL1
runFullLengthActiveL1 = arguments$FullLengthL1

##make ref and alt list to run a loop
REF <- "DMSO"
ALT <- "AZA"
REF_LIST <- c("DMSO", "DMSO", "DMSO",  "DMSO",   "DMSO",       "DMSO",        "SETDB1i",     "QSTAT")
ALT_LIST <- c("AZA",  "CKi",  "QSTAT", "SETDB1i", "QSTAT-CKi", "SETDB1i-CKi",  "SETDB1i-CKi", "QSTAT-CKi")

# REF_LIST <- c("DMSO", "DMSO", "DMSO")
# ALT_LIST <- c("AZA", "QSTAT", "CKi") 

MultiDE_conditions <- data.frame(REF = REF_LIST, ALT = ALT_LIST)
MultiDE_conditions$contrast <- paste0(MultiDE_conditions$ALT,"_",MultiDE_conditions$REF) #add the contrast column columns

# c("AZA_DMSO", "QSTAT_DMSO","SETDB1i_DMSO","CKi_DMSO", "QSTAT-CKi_DMSO", "SETDB1i-CKi_DMSO",
# "QSTAT-CKi_QSTAT","SETDB1i_QSTAT", 
# "QSTAT-CKi_CKi","SETDB1i-CKi_CKi",        
# "SETDB1i-CKi_SETDB1i","SETDB1i-CKi_QSTAT-CKi" )

# "QSTAT_CKi"

message("reading the metadata")
metadata_df <- read.csv(file = metadata_path, sep="," ,header = TRUE)




pathOverlaps <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/overlaps75SquireLINE_vs_L1_L1Base.bed"
# pathOverlaps <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/overlaps75SquireRep_vs_L1_L1Base.bed"
countsFile <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/SQuireRepeatsValidReadCounts_Aggregated_FWD_REV.bed"

dtOvs <- fread(pathOverlaps)

sampleNames <- c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3","R.C.1", "R.C.2", "R.C.3", "R.Q.1", "R.Q.2", "R.Q.3","R.QC.1", "R.QC.2", "R.QC.3", "R.S.1", "R.S.2", "R.S.3","R.SC.1", "R.SC.2", "R.SC.3")

setnames(dtOvs, c("V2","V3","V11", "V15"), c("start", "end","RNAelementID", "L1UID_name"))


message(paste0("\n\nrunning FullLengthActive? \n\n"))
#create folder for locus-specific lengths
dir.create("FullLengthActive/figures", recursive = TRUE)
dir.create("FullLengthActive/data", recursive = TRUE)




RepsCodingGenesDT <- fread("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/LocusSpecific/data/CountMatrix_squireRepeats_AggregatedFWDnREV_ProteinCoding.bed")

##get repeat ids
RepsOvCountsDT <- RepsCodingGenesDT[RNAelementID %in% dtOvs[,RNAelementID],]

colsInterested <- c("RNAelementID", sampleNames)
dtOvs <- dtOvs[RepsOvCountsDT[, ..colsInterested], on = c("RNAelementID")]


dtOvs[, RepeatWidth := (end - start) + 1]
# dtOvs[, idx := .I]
dtOvs[, isMaxWidth:=(RepeatWidth == max(RepeatWidth)), by = L1UID_name]
dtOvs[,TE_ID := RNAelementID]





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


if(methodComputeMultipleLocusSpecL1ReadCounts == "max") {
# To get other columns corresponding to these maximum rows, join back to the original data.table:
message("Selecting maximum width L1s")
dtOvs_selectReps <- dtOvs[isMaxWidth == TRUE,]
} else if (methodComputeMultipleLocusSpecL1ReadCounts == "geometericMean") {
    message("computing the geometric means of repeats that overlaps")
   #compute geometeric mean of read counts to combine multiple L1 that overlaps full length L1
dtOvs_selectReps <- dtOvs[, lapply(.SD, function(x){(as.integer(gm_mean(x)))}),by = L1UID_name, .SDcols = sampleNames]
}

dtOvs_selectReps[, gene.type := V9]

###sanity checks, which repeats are annotated LINE
RepStatsDF <- as.data.frame(table(dtOvs_selectReps$`V9`), stringsAsFactors = FALSE) %>% mutate(percentage = round(Freq/sum(Freq)*100, 2))
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

pltReps <- plotFunc(RepStatsDF, ylabIn = "Percentage", angles=90, textAngle=90, vjustIn =0.5)

ggsave(pltReps, filename = "FullLengthActive/figures/FullLengthsRepStatsDF.png", width = 10, height = 5)



##################################
##more data wrangling
#use this for master table 
dtOvs_selectReps[,RNAelementID := L1UID_name]
#### run 
ptCoding_allRepeats_activeL1 <- bind_rows(dtOvs_selectReps, RepsCodingGenesDT[gene.type == "protein_coding", ])
sum(duplicated(dtOvs_selectReps$L1UID_name))


sum(duplicated(RepsCodingGenesDT$RNAelementID))
sum(duplicated(ptCoding_allRepeats_activeL1$RNAelementID))
ptCoding_allRepeats_activeL1$RNAelementID[duplicated(ptCoding_allRepeats_activeL1$RNAelementID)]

###prepare data for active line 1 DE run
message("make matrix data for active L1")
colsDE <- c("RNAelementID", "gene.type", "L1UID_name",sampleNames) #use this for the final table
ptCoding_allRepeats_activeL1_DF <- as.data.frame(ptCoding_allRepeats_activeL1[, ..colsDE])
rownames(ptCoding_allRepeats_activeL1_DF) <- ptCoding_allRepeats_activeL1_DF$RNAelementID



######## function to run L1 full length
########################################
activeL1DE <- function(data_in, metadata_in, REF, ALT){
cond_collapsed <- paste0(c("condition",ALT,REF), collapse="_")
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

keep <- rowSums( counts(ddsInFunc) >= 10 ) >= 3
message(paste0("filtering out low counts\n",round(prop.table(table(keep)), 2)))
ddsInFunc <- ddsInFunc[keep,]

message("normalise repeats with non-repeats")
dds_main <- DESeq(ddsInFunc)

message("filter out repeats of interest")
dds_activeL1 <- dds_main[str_detect(rownames(dds_main),"UID-"),]
resActiveL1 <- results(dds_activeL1, contrast=c("condition",ALT, REF))
resActiveL1Df <- as.data.frame(resActiveL1)
resActiveL1Df$L1_ID <- row.names(resActiveL1Df)
message("summary of active L1 DE results")
message(summary(resActiveL1))
message("/n")

# shrink results for volcano
resultsActiveL1Shrink <- lfcShrink(dds_activeL1, contrast = c("condition",ALT, REF), res=resActiveL1, type = 'normal')

message("making volcano plots")
volActiveL1Shrink <- plotEVolcano(resultsActiveL1Shrink, title_in=cond_collapsed, subtitle_in="Differential Active L1 Expression", REF=REF, ALT=ALT)
volActiveL1NoShrink <- plotEVolcano(resActiveL1, title_in=cond_collapsed, subtitle_in="Differential Active L1 Expression", REF=REF, ALT=ALT)

ggsave(volActiveL1Shrink, filename = paste0("FullLengthActive/figures/volcanoShrink_ActiveLINE1_",cond_collapsed,".svg"))
ggsave(volActiveL1Shrink, filename = paste0("FullLengthActive/figures/volcanoShrink_ActiveLINE1_",cond_collapsed,".pdf"))
ggsave(volActiveL1NoShrink, filename = paste0("FullLengthActive/figures/volcanoNoShrink_ActiveLINE1_",cond_collapsed,".svg"))
ggsave(volActiveL1NoShrink, filename = paste0("FullLengthActive/figures/volcanoNoShrink_ActiveLINE1_",cond_collapsed,".pdf"))

# make CPM
message("making CPM")
ActiveL1_cpm_df <- makeCPM(DeSeq2_res=dds_activeL1, prior_cnts=2, keyColName="RNAelementID")

#save results
message("saving results to disk")
results2SaveActiveL1 <- list(ddsMain=dds_main, resultsObject=resActiveL1, resDF=resActiveL1Df, resShrinkedDF=resultsActiveL1Shrink, cpmDF=ActiveL1_cpm_df, volcano_ggplot = volActiveL1Shrink)
save(results2SaveActiveL1, file=paste0("FullLengthActive/data/Differential_ActiveL1ExpressionResults_",cond_collapsed,".RData"))  # Till here everything is ok. 
return(resActiveL1Df)
}



###### call main function
resList <- list()
##run
if(runFullLengthActiveL1 == TRUE){
  message("running locus specific l1")
resList <- lapply(1:nrow(MultiDE_conditions), function(x) {activeL1DE(data_in = ptCoding_allRepeats_activeL1_DF, metadata_in = metadata_df, REF = MultiDE_conditions[x,1], ALT = MultiDE_conditions[x,2])})
}else{
  message("skipping full length L1 analysis")
}

#### further analysis
names(resList) <- MultiDE_conditions$contrast
#resList[["SETDB1i_DMSO"]]
resAllSamplesDF <- rbindlist(resList, idcol = "condition", use.names = FALSE)

##Add categories 
resultsDTSample <- resAllSamplesDF %>% mutate(LFClabel = case_when(log2FoldChange > 0 & padj <= arguments$lfcCutoff ~ "up",
                                     log2FoldChange < 0 & padj <= arguments$lfcCutoff ~ "down",
                                     TRUE ~  "no_change")) 


#plot Up and Down counts
df_filtered <- resultsDTSample %>%
  filter(LFClabel %in% c("up", "down")) %>%
  group_by(condition, LFClabel) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(count = ifelse(LFClabel == "down", -count, count))  # Make "down" counts negative


df_ordered <- df_filtered %>%
  filter(LFClabel == "up") %>%
  arrange(desc(count)) %>%
  pull(condition) %>%
  unique()

# Apply this factor level ordering to the full dataset
df_filtered$condition <- factor(df_filtered$condition, levels = df_ordered)

# Plot with ggplot2
plt <- ggplot(df_filtered, aes(x = condition, y = count, fill = LFClabel)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "Up and Down Sig. Full Length Active Expression per condition", x = "condition", y = "Count", fill = "Label") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plt, file = paste0("FullLengthActive/figures/up_down_counts_per_Condition.png"), width = 10, height = 6)

# ggsave(pltReps, filename = "FullLengthActive/figures/FullLengthsRepStatsDF.png", width = 10, height = 5)


######Upset plot
df_binary <- resultsDTSample %>% filter(LFClabel %in% c("up", "down")) %>% 
  select(!c("lfcSE","stat","pvalue","padj","LFClabel", "baseMean", "log2FoldChange")) %>%  
  mutate(value = 1) %>% pivot_wider(names_from = condition, values_from = value, values_fill = 0L)



# Ensure L1_ID is not used as a set
df_matrix <- df_binary %>%
  column_to_rownames("L1_ID")  # Convert L1_ID to row names

# Convert to a standard dataframe (UpSetR does not like tibbles)
df_matrix <- as.data.frame(df_matrix)

png("FullLengthActive/figures/Upset_sigRepeats_condition.png", width = 13, height = 6, units = "in", res = 300)
# Generate and display the UpSet plot
upset(df_matrix, 
      sets = colnames(df_matrix), 
      nsets = length(colnames(df_matrix)), 
      nintersects = 40, 
      keep.order = TRUE)
dev.off()




#############HEATMAPS
##Heatmap
log2FoldChangeDF <- resultsDTSample %>% filter(LFClabel %in% c("up", "down")) %>% select(!c("lfcSE","stat","pvalue","padj","LFClabel", "baseMean")) %>% pivot_wider(names_from = condition, values_from = log2FoldChange) #%>% 
log2FoldChangeDF <- log2FoldChangeDF %>% column_to_rownames("L1_ID") 
# log2FoldChangeDF[is.na(log2FoldChangeDF)] = 0

pdf(file = paste0("FullLengthActive/figures/Heatmaplog2FoldChangeConditions.pdf"), width=12, height=9)
pheatmap::pheatmap(log2FoldChangeDF, main = paste0("abs(LFC) > 1 & padj < 0.05"), 
scale = "row", 
annotation_names_row=FALSE,
show_rownames = TRUE, 
cluster_rows = FALSE,
cluster_cols = FALSE,
na_col = "black")
dev.off()

log2FoldChangeZeroNADF <- log2FoldChangeDF
log2FoldChangeZeroNADF[is.na(log2FoldChangeZeroNADF)] = 0

# Define color palette: blue to white to red
my_palette <- colorRampPalette(c("blue", "white", "red"))(20)

# Set color breaks centered at 0
max_val <- max(abs(log2FoldChangeZeroNADF), na.rm = TRUE)
breaks <- seq(-max_val, max_val, length.out = 21)


pdf(file = paste0("FullLengthActive/figures/Heatmaplog2FoldChangeConditionsClustering.pdf"), width=12, height=9)
pheatmap::pheatmap(log2FoldChangeZeroNADF, main = paste0("abs(LFC) > 1 & padj < 0.05"), 
scale = "row", 
annotation_names_row=FALSE,
show_rownames = TRUE, 
cluster_rows = TRUE,
cluster_cols = TRUE,
color = my_palette,
breaks = breaks)
dev.off()


##write data to disk
fwrite(resultsDTSample, file = paste0("FullLengthActive/data/merged_FullLengthActiveL1_DESeq2Results.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

message("done\n")
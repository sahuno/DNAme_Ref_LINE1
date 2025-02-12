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
#read squire data
counts_df <- fread(squireCounts_path)
counts_annot <- read_tsv(paste0(counts_annot_path))
counts_annot <- counts_annot[counts_annot$`gene.type`=="protein_coding",]
# table(counts_annot$`gene.type`)


counts_repeats_protein_df <- rbind(counts_df, counts_annot, fill=TRUE)
counts_repeats_protein_df <- counts_repeats_protein_df %>% relocate(c(
    "R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3",
    "R.C.1", "R.C.2", "R.C.3", "R.Q.1", "R.Q.2", "R.Q.3",
    "R.QC.1", "R.QC.2", "R.QC.3", "R.S.1", "R.S.2", "R.S.3",
    "R.SC.1", "R.SC.2", "R.SC.3"), .after = last_col())

setDT(counts_repeats_protein_df)

#counts_repeats_protein_df[row]

#get rowsums that are above threshold
keepDF <- rowSums(counts_repeats_protein_df[,c(
    "R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3",
    "R.C.1", "R.C.2", "R.C.3", "R.Q.1", "R.Q.2", "R.Q.3",
    "R.QC.1", "R.QC.2", "R.QC.3", "R.S.1", "R.S.2", "R.S.3",
    "R.SC.1", "R.SC.2", "R.SC.3")] >= MIN_ReadsCounts) >= smallestGroupSize


counts_repeats_protein_df <- counts_repeats_protein_df[keepDF,]

#################
#create separate id for for repeats and genes
counts_repeats_protein_df[,genomicElementID:=(fcase(!is.na(te.id), te.id,
                                                    !is.na(gene.id), gene.id ) )]


dt_4DE <- counts_repeats_protein_df[, c("genomicElementID", "R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3")]
df_4DE <- as.data.frame(dt_4DE) 
rownames(df_4DE) <- df_4DE$genomicElementID
df_4DE <- df_4DE[, -1]

df_4DE <- df_4DE %>% mutate(across(c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3"), as.integer))


#read the metadata
metadata_df <- read.csv(file = metadata_path, sep="," ,header = TRUE)

#drop samples if necesssary
metadata_df <- metadata_df %>% filter(str_detect(samples, "^R.A|^R.0")) #%>% mutate(condition)
metadata_df <- metadata_df %>% dplyr::filter(condition %in% c("DMSO", "AZA"))

#################################
#################################
# run differenrial gene expression
dds <- DESeqDataSetFromMatrix(countData = data.matrix(df_4DE),
                              colData = metadata_df,
                              design = ~ condition)


dds <- DESeq(dds)

# rownames(dds)[str_detect(rownames(dds),"L1")]
### filter L1 only
dds_L1 <- dds[str_detect(rownames(dds),"L1"),]
ddsResults <- results(dds_L1)

#shrink results for volcano
results_L1Shrink <- lfcShrink(dds_L1, contrast = c("condition",'AZA','DMSO'), res=ddsResults, type = 'normal')


############################
#volcano plots of line 1
volcano_ggplot <- EnhancedVolcano(results_L1Shrink,
    lab = rownames(results_L1Shrink),
    x = 'log2FoldChange',
    y = 'padj',
    title = '5-AZA versus DMSO',
    subtitle = "Differential L1 expression",
    col=c('grey', 'grey', 'grey', 'blue'),
    ylab = bquote(~-Log[10]~ 'padjust'))
ggsave(volcano_ggplot, filename = "volcano_ggplot_LINE1.png")


results2Save <- list(dds, dds_L1, ddsResults, results_L1Shrink)
save(results2Save, file="DifferentialL1ExpressionResults.RData")  # Till here everything is ok. 

#save copies of results
ddsResultsDf <- as.data.frame(ddsResults)
summary(ddsResults)



#get expression matrix
dds_sf_Disp_L1_cpm <- edgeR::cpm(dds_L1, log = TRUE, prior.count = 2) #prior.count = 2 is added to avoid log(0) error
dim(dds_sf_Disp_L1_cpm)
dds_sf_Disp_L1_cpm_df <- as.data.frame(dds_sf_Disp_L1_cpm)
dds_sf_Disp_L1_cpm_df <- dds_sf_Disp_L1_cpm_df %>% rownames_to_column(var = "te.id")

######################################################################################
#################### create key; DNAme and RNA
###########################################################################################
# read full length line1 from L1Base
l1Base_entireLength <- fread(l1_path, col.names = c("chrom", "start", "end", "RepeatID", "score", "strand", "start2", "end2", "color"))



#filter out L1 from squire counts table 
counts_df_L1 <- counts_df %>% dplyr::filter(str_detect(te.id, "L1")) %>% dplyr::select(`te.id`,`te.name`) #!str_detect(te.id, "HLA")
counts_df_L1 <- counts_df_L1 %>% separate_wider_delim(te.name, "|" , names = c("chrom", "start", "end", "name","number","strand")) %>% separate_wider_delim(name, ":" ,names = c("consensus","clade", "class"))

# table(counts_df_L1$consensus)

counts_gr_L1 <- makeGRangesFromDataFrame(counts_df_L1, keep.extra.columns = TRUE)
#plot l1 width distribution
counts_gr_L1$SeqWidthSQuIRE <- width(counts_gr_L1)

#### find overlaps between l1base a 
l1Base_entireLength_gr <- makeGRangesFromDataFrame(l1Base_entireLength, keep.extra.columns = TRUE)
l1Base_entireLength_gr$SeqWidth_L1Base <- width(l1Base_entireLength_gr) #add width 

# filter standard chromosomes only
standard_chromosomes <- paste0("chr", c(1:19, "X", "Y"))  # Exclude "chrM" if not needed
counts_gr_L1 <- keepSeqlevels(counts_gr_L1, standard_chromosomes, pruning.mode = "coarse")



lsRNA_l1Base_keys <- gUtils::gr.findoverlaps(l1Base_entireLength_gr, subject=counts_gr_L1, first=FALSE,  qcol = c("RepeatID", "SeqWidth_L1Base"), scol=c("te.id", "consensus","clade", "class", "number", "SeqWidthSQuIRE"), return.type = "data.table")
##Note that start & End in `lsRNA_l1Base_keys` corresponds to that of Squire Cordinates

# lsRNA_l1Base_keys_minOv50 <- gUtils::gr.findoverlaps(l1Base_entireLength_gr, subject=counts_gr_L1, minoverlap = 10, first=FALSE,  qcol = c("RepeatID", "SeqWidth_L1Base"), scol=c("te.id", "consensus","clade", "class", "number", "SeqWidthSQuIRE"), return.type = "data.table")


##how many l1base cordiantes overlaps the rna counts data from squire; seems there are repeats `UID-1` mapping multiple loci
lsRNA_l1Base_keys %>% group_by(RepeatID) %>% summarise(n = n()) %>% dplyr::filter(n > 1)
dim(lsRNA_l1Base_keys)

# has_overlap <- gr.in(lsRNA_l1Base_keys, lsRNA_l1Base_keys)


#merged adjacent regions
# Convert data.table to GRanges
lsRNA_l1Base_keys_gr <- dt2gr(lsRNA_l1Base_keys)
# Apply gr.reduce() with grouping by "RepeatID"
lsRNA_l1Base_keys_reduced_gr <- gr.reduce(lsRNA_l1Base_keys_gr, by = "RepeatID", ignore.strand = TRUE, span = FALSE)
head(table(mcols(lsRNA_l1Base_keys_reduced_gr)$RepeatID), 20)


#########################
### convert to bedfile
# Convert GRanges to a data frame
# source('/DNAme_Ref_LINE1/scripts/R/gr2bed.R')


##pivort longer for ploting
lsRNA_l1Base_keys_dt <- gUtils::gr2dt(lsRNA_l1Base_keys)
lsRNA_l1Base_keys_dt_longer <- pivot_longer(lsRNA_l1Base_keys_dt, cols = c(SeqWidthSQuIRE, SeqWidth_L1Base))

l1baseSquireSeqLengths_histogram <- ggplot(lsRNA_l1Base_keys_dt_longer, aes(value) ) + geom_histogram() + facet_wrap(~name)
ggsave(l1baseSquireSeqLengths_histogram, filename = "overlaps_l1baseSquireSeqLengths_histogram.pdf")


# lsRNA_keys %>% filter(!str_detect(seqnames, "chrY|chrM|chrX"))
######################################################################################
######################################################################################

head(dds_sf_Disp_L1_cpm_df); message('\n'); head(lsRNA_l1Base_keys)
dds_sf_Disp_L1_cpm_df_withRepeatID <- dds_sf_Disp_L1_cpm_df %>% left_join(lsRNA_l1Base_keys %>% dplyr::select(c("RepeatID","te.id", "consensus","clade", "class")), by="te.id") %>% as_tibble()

dds_sf_Disp_L1_cpm_df_withRepeatID <- dds_sf_Disp_L1_cpm_df_withRepeatID %>% filter(!is.na(RepeatID))

dds_sf_Disp_L1_cpm_df_withRepeatID_piv <- dds_sf_Disp_L1_cpm_df_withRepeatID %>% pivot_longer(col = starts_with("R."), names_to="RNAsamples", values_to = "cpm") %>% mutate(samples = str_replace_all(RNAsamples,"\\.","-"), samples = str_replace_all(samples,"R","D"))
# head(dds_sf_Disp_L1_cpm_df_withRepeatID)
#sanity checks
# dds_sf_Disp_L1_cpm_df_withRepeatID_piv %>% filter(!is.na(RepeatID))



####################################################################################################
################## # see volcano of active line 1 ##################
###Que: among all the line 1 from squire, that overlapped l1base, do differential expression
dds_ActiveL1 <- dds[rownames(dds) %in% unique(lsRNA_l1Base_keys[,`te.id`]),]

#########################################################################
############################################################################


### read dna methyltaion rates data;  overlaps already been done

#read rates
DNAme_rates_dt <- read_fst(path_DNAmeRates, as.data.table = TRUE)
# DNAme_rates_dt[variable=="Freq_5mCpG.N" & value >= 3.0,]
DNAme_rates_geom_mean_dt_ <- DNAme_rates_dt[variable=="Freq_5mCpG.geom_mean",]


DNAme_rates_dt_statsWider <- dcast(DNAme_rates_dt, ... ~ variable, value.var = "value")[`Freq_5mCpG.N` >= 3,]
DNAme_rates_dt_statsWiderGmean <- DNAme_rates_dt_statsWider[,!c("Freq_5mCpG.N", "Freq_5mCpG.mean", "Freq_5mCpG.median", "Freq_5mCpG.entropy_in_log2")]
# summary($)

DNAme_rates_wide_dt <- dcast(DNAme_rates_dt_statsWiderGmean[,!c("id")], ... ~ samples, value.var = "Freq_5mCpG.geom_mean")

# DNAme_rates_wide_dt <- dcast(DNAme_rates_dt[variable=="Freq_5mCpG.geom_mean",][,!c("variable", "id")], ... ~ samples, value.var = "value")
# table(DNAme_rates_dt$variable)

# list_paths_overlaps <- list.files("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/results/OverlapsL1PromoterDNAme/", recursive = TRUE, full.names = TRUE, pattern=".bed")
# list_paths_overlaps <- list.files("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/results/OverlapsL1PromoterDNAme/", recursive = TRUE, full.names = TRUE, pattern=".bed")
# paths_fullLength_nonCpGIs <- list_paths_overlaps[str_detect(list_paths_overlaps,"fullLength/nonCGI")]
# paths_fullLength_nonCpGIs <- paths_fullLength_nonCpGIs[str_detect(paths_fullLength_nonCpGIs,"_minCov5")]

message("\nfilenames of under consideration\n")
# message(paths_fullLength_nonCpGIs)
message("end of file names")



# DNAmeRatesOverlaps <- lapply(paths_fullLength_nonCpGIs, function(x){fread(x, col.names = c("chrom", "start", "end", "RepeatID", "score", "strand", "start2", "end2", "color", "min", "max", "mean", "median", "count"))})
# names(DNAmeRatesOverlaps) <- basename(paths_fullLength_nonCpGIs)

# DNAmeRatesOverlaps <- rbindlist(DNAmeRatesOverlaps, idcol = "samples")

# DNAmeRatesOverlaps[, minCov:=(str_extract(samples, "minCov[0-9]*"))]
# OverlapsDNAmeRNA_df[, samples:=(str_extract(samples, ".*(?=_5mCpG)"))]
# DNAmeRatesOverlaps[, samples:=(str_extract(samples, "^[^_]+"))]

#convert to numeric to enable ploting of columns 
# DNAmeRatesOverlaps <- DNAmeRatesOverlaps %>% mutate(across(c("min",   "max",    "mean", "median" , "count"), as.numeric))


#sanity checks, how many line with valid results
# nRepeatswithValidDNAmeRateBed <- DNAmeRatesOverlaps %>% filter(count > 0) %>% group_by(samples) %>% summarise(n=n())
message("\nBegin valid overlaps DNAme")
# nRepeatswithValidDNAmeRateBed
message("done valid repeats DNAme\n")

# repeats with more than 3 cpgs
# DNAmeRatesOverlaps %>% filter(count > 3) %>% group_by(samples) %>% summarise(n=n())

##merge Differential repeats (CPM) with DNAmerates computed with bedtools, precedence is given to `dds_sf_Disp_L1_cpm_df_withRepeatID_piv` entries
DNAmeRatesOverlaps_joined <- dds_sf_Disp_L1_cpm_df_withRepeatID_piv %>% left_join(DNAme_rates_dt_statsWiderGmean, by = c("samples", "RepeatID"))

# DNAmeRatesOverlaps_joined <- DNAmeRatesOverlaps_joined %>% mutate(across(c("min",   "max",    "mean", "median" , "count"), as.numeric))

# DNAmeRatesOverlaps_joined %>% filter(is.na(count))
###keep only overlapping repeats
DNAmeRatesOverlaps_joined <- DNAmeRatesOverlaps_joined %>% filter(!is.na(RepeatID))

DNAmeRatesOverlaps_joined <- DNAmeRatesOverlaps_joined %>% left_join(metadata_df, by= c("RNAsamples" = "samples"))


ctsRepeatsPerSample <- DNAmeRatesOverlaps_joined %>% group_by(RepeatID) %>% summarize(cntsRepeats = n()) %>% print(n = 60) 
L1_multiMapping <- ctsRepeatsPerSample %>% filter(cntsRepeats > 6) %>% pull(RepeatID)
DNAmeRatesOverlaps_joined %>% filter(RepeatID %in% L1_multiMapping) %>% dplyr::select(`te.id`,RepeatID) %>% distinct() %>% arrange(RepeatID) %>% dplyr::select(RepeatID, `te.id`) %>% group_by(RepeatID) %>% summarise(countsRNALoci = n())
# group_by(RepeatID, `te.id`) %>% duplicated() #print(n = 100) 

DNAmeRatesOverlaps_joined %>% group_by(te.id) %>% summarize(cntsRepeats = n())
# DNAmeRatesOverlaps_joined %>% filter(is.na(mean)) %>% group_by(samples) 

message("repeats counts per sample")
ctsRepeatsPerSample
# DNAmeRatesOverlaps_joined
###############
# DNAmeRatesOverlaps_joined$mean <- as.numeric(DNAmeRatesOverlaps_joined$mean)
# DNAmeRatesOverlaps_joined$RepeatID <- paste0(DNAmeRatesOverlaps_joined$chrom, ":",DNAmeRatesOverlaps_joined$start,"_", DNAmeRatesOverlaps_joined$end, ":", DNAmeRatesOverlaps_joined$RepeatID, ":",DNAmeRatesOverlaps_joined$consensus)
DNAmeRatesOverlaps_joined$RepeatID <- paste0(DNAmeRatesOverlaps_joined$id, "_",DNAmeRatesOverlaps_joined$RepeatID, ":", DNAmeRatesOverlaps_joined$te.id)



##write table to disk
write_tsv(DNAmeRatesOverlaps_joined, file = "L1Base_DNAme_locusSpecificL1_RNA.tsv")

RepeatIDGroups <- split(DNAmeRatesOverlaps_joined, as.factor(DNAmeRatesOverlaps_joined$RepeatID))
# names(RepeatIDGroups)

# plotL1perCondition <- function(l1Name){
# plotOut <- ggplot(RepeatIDGroups[[l1Name]], aes(x=mean, y=cpm, color=condition)) + geom_point() + 
# geom_text_repel(aes(label=samples)) + labs(title = paste(l1Name), y = "cpm", x = "mean DNAme") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
# }


# plotL1SingleCondition <- function(l1Name){
# plotOut <- ggplot(RepeatIDGroups[[l1Name]] %>% dplyr::filter(condition=="DMSO"), aes(x=mean, y=cpm)) + geom_point() + 
# geom_text_repel(aes(label=samples)) + labs(title = paste(l1Name), y = "cpm", x = "mean DNAme") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
# }

# # dir.create("results/boxplotsRepPerCond/")
# pdf(paste0("scatter_DMSOAZA_lspecific_L1RNA_DNAme",".pdf"))
# lapply(names(RepeatIDGroups), plotL1perCondition)
# dev.off()

# pdf(paste0("scatter_DMSO_lspecific_L1RNA_DNAme",".pdf"))
# lapply(names(RepeatIDGroups), plotL1SingleCondition)
# dev.off()

# dmso color; #color="#00BFC4"

plotOutAll <- ggplot(DNAmeRatesOverlaps_joined %>% dplyr::filter(condition=="DMSO"), aes(x=`Freq_5mCpG.geom_mean`, y=cpm, color = samples)) + geom_point() + 
labs(title = "Active line1 from L1Base & locus Specific RNA", y = "RNA/cpm", x = "DNA methylation/mean 5mC") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

plotOutperSample <- ggplot(DNAmeRatesOverlaps_joined %>% dplyr::filter(condition=="DMSO"), aes(x=`Freq_5mCpG.geom_mean`, y=cpm, color = samples)) + geom_point() + facet_wrap(~samples, scale = "free")+
labs(title = "Active line1 from L1Base & locus Specific RNA", y = "RNA/cpm", x = "DNA methylation/mean 5mC") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))


plotOutAllConditions <- ggplot(DNAmeRatesOverlaps_joined , aes(x=`Freq_5mCpG.geom_mean`, y=cpm, color=condition)) + geom_point() + 
labs(title = "Active line1 from L1Base & locus Specific RNA", y = "RNA/cpm", x = "DNA methylation/mean 5mC") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

ggsave(plotOutAll, filename = "scatter_combined_DMSO_specific_L1RNA_DNAme.png")
ggsave(plotOutAllConditions, filename = "scatter_combined_DMSO_5Aza_specific_L1RNA_DNAme.png")
ggsave(plotOutperSample, filename = "scatter_combined_DMSO_specific_L1RNA_DNAme_perSample.png")


# ##plot variable genes 
# rld <- rlog(dds, blind=FALSE) #use, dseq2 object with all samples/conditions/contrasts, doeen't matter since it's experiment-wide investigation
# topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 10) #top variable genes across all samples
# rld_top  <- assay(rld)[topVarGenes, ]
# rld_top  <- rld_top - rowMeans(rld_top)
# anno <- as.data.frame(colData(rld)[, c("samples","condition")])




## make a unique id for repeats and genes
# length(keepDF)
# counts_repeats_protein_df %>% rowwise() %>% mutate(Rowsums = sum(across(starts_with("R."))))


# R.0.1, R.0.2, R.0.3, R.A.1, R.A.2, R.A.3, R.C.1, R.C.2, R.C.3, R.Q.1, R.Q.2, R.Q.3, R.QC.1, R.QC.2, R.QC.3, R.S.1, R.S.2, R.S.3, R.SC.1, R.SC.2, R.SC.3 



#############################
#############################
metadata_df <- metadata_df %>% mutate(new_samples_name = str_replace_all(new_samples_name,"R","D")) 

DNAmeRatesOverlaps_2plot <- DNAmeRatesOverlaps_joined %>% dplyr::filter(samples %in% metadata_df$new_samples_name) %>% left_join(metadata_df, by= c("samples"="new_samples_name"))
DNAmeRatesOverlaps_2plot_filtered <- DNAmeRatesOverlaps_2plot  %>% dplyr::filter(`condition.x`=="DMSO")

plotHist <- ggplot(data=DNAmeRatesOverlaps_2plot_filtered, aes(Freq_5mCpG.geom_mean)) + geom_histogram() + facet_wrap(~samples) + theme_minimal() +
labs(x = "mean DNAme", title = "DNA methylation distribution of murine active full length L1 - DMSO")
ggsave(plotHist, filename="histogram_DNAme.png", width = 9, height = 7) 

#############################
#############################

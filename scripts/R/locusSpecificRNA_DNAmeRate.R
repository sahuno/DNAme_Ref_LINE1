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

options("width"=200)


## To  run
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/locusSpecificRNA_DNAmeRate.R

#### program options
ref_variable = "DMSO"
MIN_ReadsCounts = 10
smallestGroupSize <- 3

metadata_path <- "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_recoded.csv"
l1_path <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed"
dir_overlaps <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/results/OverlapsL1PromoterDNAme/"

samplesOfChoice <- c("R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3")
sampleList <- c(
    "R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3",
    "R.C.1", "R.C.2", "R.C.3", "R.Q.1", "R.Q.2", "R.Q.3",
    "R.QC.1", "R.QC.2", "R.QC.3", "R.S.1", "R.S.2", "R.S.3",
    "R.SC.1", "R.SC.2", "R.SC.3")


#read squire data
counts_df <- fread("/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/CT/squire_te_fwd.tsv")
counts_annot <- read_tsv(paste0('/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/CT/counts_annot.tsv'))
dim(counts_df)
counts_repeats_protein_df <- rbind(counts_df, counts_annot, fill=TRUE)

counts_repeats_protein_df <- counts_repeats_protein_df %>% relocate(sampleList, .after = last_col())

setDT(counts_repeats_protein_df)

#counts_repeats_protein_df[row]
dim(counts_repeats_protein_df)
#get rowsums that are above minimum read thresholds
keepDF <- rowSums(counts_repeats_protein_df[, ..sampleList, with = FALSE] >= MIN_ReadsCounts) >= smallestGroupSize


# keepDF <- rowSums(counts_repeats_protein_df[,c(
#     "R.0.1", "R.0.2", "R.0.3", "R.A.1", "R.A.2", "R.A.3",
#     "R.C.1", "R.C.2", "R.C.3", "R.Q.1", "R.Q.2", "R.Q.3",
#     "R.QC.1", "R.QC.2", "R.QC.3", "R.S.1", "R.S.2", "R.S.3",
#     "R.SC.1", "R.SC.2", "R.SC.3")] >= MIN_ReadsCounts) >= smallestGroupSize


counts_repeats_protein_df <- counts_repeats_protein_df[keepDF,]

#create separate id for for repeats and genes
counts_repeats_protein_df[,genomicElementID:=(fcase(!is.na(te.id), te.id,
                                                    !is.na(gene.id), gene.id ) )]


# table(unique(counts_repeats_protein_df$`gene.type`))
cols2Filter <- c("genomicElementID", samplesOfChoice)
dt_4DE <- counts_repeats_protein_df[, ..cols2Filter]
df_4DE <- as.data.frame(dt_4DE) 
rownames(df_4DE) <- df_4DE$genomicElementID
df_4DE <- df_4DE[, -1]

df_4DE <- df_4DE %>% mutate(across(samplesOfChoice, as.integer))
dim(df_4DE)

#read the metadata
metadata_df <- read.csv(file = metadata_path, sep=",", header = TRUE)
#drop samples if necesssary
metadata_df <- metadata_df %>% filter(str_detect(samples, "^R.A|^R.0")) #%>% mutate(condition)
metadata_df <- metadata_df %>% dplyr::filter(condition %in% c("DMSO", "AZA"))

#################################
# run differenrial gene expression
#################################
#create DeSeq2 object
dds <- DESeqDataSetFromMatrix(countData = data.matrix(df_4DE),
                              colData = metadata_df,
                              design = ~ condition)


dds <- DESeq(dds)

# rownames(dds)[str_detect(rownames(dds),"L1")]
### filter L1 only
dds_L1 <- dds[str_detect(rownames(dds), "L1"), ]
dim(dds_L1)
ddsResults <- results(dds_L1)

ddsResultsDf <- as.data.frame(ddsResults)
summary(ddsResults)

dds_sf_Disp_L1_cpm <- edgeR::cpm(dds_L1, log = TRUE, prior.count = 2) #prior.count = 2 is added to avoid log(0) error
dim(dds_sf_Disp_L1_cpm)
dds_sf_Disp_L1_cpm_df <- as.data.frame(dds_sf_Disp_L1_cpm)
dds_sf_Disp_L1_cpm_df <- dds_sf_Disp_L1_cpm_df %>% rownames_to_column(var = "te.id")

######################################################################################
#################### create key; DNAme and RNA
###########################################################################################
# read full length line1 from L1Base
l1Base_entireLength <- fread(l1_path, col.names = c("chrom", "start", "end", "RepeatID", "score", "strand", "start2", "end2", "color"))

## the goal is to get cordinates of all locus-specific line1 irrespective of read counts
#filter out L1 from squire counts table 
counts_df_L1 <- counts_df %>% dplyr::filter(str_detect(te.id, "L1")) %>% dplyr::select(`te.id`,`te.name`) #!str_detect(te.id, "HLA")
counts_df_L1 <- counts_df_L1 %>% separate_wider_delim(te.name, "|" , names = c("chrom", "start", "end", "name","number","strand")) %>% separate_wider_delim(name, ":" ,names = c("consensus","clade", "class"))
dim(counts_df_L1)

counts_gr_L1 <- makeGRangesFromDataFrame(counts_df_L1, keep.extra.columns = TRUE)

###Notes 
#### find overlaps between l1base a 
gr_l1Base_entireLength <- makeGRangesFromDataFrame(l1Base_entireLength, keep.extra.columns = TRUE)
lsRNA_l1Base_keys <- gUtils::gr.findoverlaps(gr_l1Base_entireLength, subject=counts_gr_L1, first=FALSE,  qcol = c("RepeatID"), scol=c("te.id", "consensus","clade", "class", "number"), return.type = "data.table")

##how many l1base cordiantes overlaps the rna counts data from squire; seems there are repeats `UID-1` mapping multiple loci
lsRNA_l1Base_keys %>% group_by(RepeatID) %>% summarise(n())
dim(lsRNA_l1Base_keys)
# lsRNA_keys %>% filter(!str_detect(seqnames, "chrY|chrM|chrX"))
######################################################################################
######################################################################################

head(dds_sf_Disp_L1_cpm_df); message('\n');head(lsRNA_l1Base_keys)
dds_sf_Disp_L1_cpm_df_withRepeatID <- dds_sf_Disp_L1_cpm_df %>% left_join(lsRNA_l1Base_keys %>% dplyr::select(c("RepeatID","te.id", "consensus","clade", "class")), by="te.id") %>% as_tibble()

dds_sf_Disp_L1_cpm_df_withRepeatID <- dds_sf_Disp_L1_cpm_df_withRepeatID %>% filter(!is.na(RepeatID))

dds_sf_Disp_L1_cpm_df_withRepeatID_piv <- dds_sf_Disp_L1_cpm_df_withRepeatID %>% pivot_longer(col = starts_with("R."), names_to="RNAsamples", values_to = "cpm") %>% mutate(samples = str_replace_all(RNAsamples,"\\.","-"), samples = str_replace_all(samples,"R","D"))
# head(dds_sf_Disp_L1_cpm_df_withRepeatID)
#sanity checks
# dds_sf_Disp_L1_cpm_df_withRepeatID_piv %>% filter(!is.na(RepeatID))

########################################
### read overlaps data
list_paths_overlaps <- list.files(dir_overlaps, recursive = TRUE, full.names = TRUE, pattern=".bed")
paths_fullLength_nonCpGIs <- list_paths_overlaps[str_detect(list_paths_overlaps,"fullLength/nonCGI")]
paths_fullLength_nonCpGIs <- paths_fullLength_nonCpGIs[str_detect(paths_fullLength_nonCpGIs,"_minCov5")]

message("\nfilenames of under consideration\n")
message(paths_fullLength_nonCpGIs)
message("end of file names")


DNAmeRatesOverlaps <- lapply(paths_fullLength_nonCpGIs, function(x){fread(x, col.names = c("chrom", "start", "end", "RepeatID", "score", "strand", "start2", "end2", "color", "min", "max", "mean", "median", "count"))})
names(DNAmeRatesOverlaps) <- basename(paths_fullLength_nonCpGIs)

DNAmeRatesOverlaps <- rbindlist(DNAmeRatesOverlaps, idcol = "samples")

DNAmeRatesOverlaps[, minCov:=(str_extract(samples, "minCov[0-9]*"))]
# OverlapsDNAmeRNA_df[, samples:=(str_extract(samples, ".*(?=_5mCpG)"))]
DNAmeRatesOverlaps[, samples:=(str_extract(samples, "^[^_]+"))]

#convert to numeric to enable ploting of columns 
DNAmeRatesOverlaps <- DNAmeRatesOverlaps %>% mutate(across(c("min",   "max",    "mean", "median" , "count"), as.numeric))


#sanity checks, how many line with valid results
nRepeatswithValidDNAmeRateBed <- DNAmeRatesOverlaps %>% filter(count > 0) %>% group_by(samples) %>% summarise(n=n())
message("\nBegin valid overlaps DNAme")
nRepeatswithValidDNAmeRateBed
message("done valid repeats DNAme\n")

# repeats with more than 3 cpgs
DNAmeRatesOverlaps %>% filter(count > 3) %>% group_by(samples) %>% summarise(n=n())

##merge Differential repeats (CPM) with DNAmerates computed with bedtools, precedence is given to `dds_sf_Disp_L1_cpm_df_withRepeatID_piv` entries
DNAmeRatesOverlaps_joined <- dds_sf_Disp_L1_cpm_df_withRepeatID_piv %>% left_join(DNAmeRatesOverlaps, by = c("samples", "RepeatID"))

DNAmeRatesOverlaps_joined <- DNAmeRatesOverlaps_joined %>% mutate(across(c("min",   "max",    "mean", "median" , "count"), as.numeric))

DNAmeRatesOverlaps_joined %>% filter(is.na(count))
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
DNAmeRatesOverlaps_joined$RepeatID <- paste0(DNAmeRatesOverlaps_joined$chrom, ":",DNAmeRatesOverlaps_joined$start,"_", DNAmeRatesOverlaps_joined$end, ":", DNAmeRatesOverlaps_joined$RepeatID, ":",DNAmeRatesOverlaps_joined$consensus)



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

plotOutAll <- ggplot(DNAmeRatesOverlaps_joined %>% dplyr::filter(condition=="DMSO"), aes(x=mean, y=cpm, color = samples)) + geom_point() + 
labs(title = "Active line1 from L1Base & locus Specific RNA", y = "RNA/cpm", x = "DNA methylation/mean 5mC") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

plotOutAllConditions <- ggplot(DNAmeRatesOverlaps_joined , aes(x=mean, y=cpm, color=condition)) + geom_point() + 
labs(title = "Active line1 from L1Base & locus Specific RNA", y = "RNA/cpm", x = "DNA methylation/mean 5mC") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))


ggsave(plotOutAll, filename = "scatter_combined_DMSO_specific_L1RNA_DNAme.png")
ggsave(plotOutAllConditions, filename = "scatter_combined_DMSO_5Aza_specific_L1RNA_DNAme.png")


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

DNAmeRatesOverlaps_2plot <- DNAmeRatesOverlaps %>% dplyr::filter(samples %in% metadata_df$new_samples_name) %>% left_join(metadata_df, by= c("samples"="new_samples_name"))
DNAmeRatesOverlaps_2plot_filtered <- DNAmeRatesOverlaps_2plot %>% filter(count >= 3) %>% dplyr::filter(condition=="DMSO")

plotHist <- ggplot(data=DNAmeRatesOverlaps_2plot_filtered, aes(mean)) + geom_histogram() + facet_wrap(~samples) + theme_minimal() +
labs(x = "mean DNAme", title = "DNA methylation distribution of murine active full length L1 - DMSO")
ggsave(plotHist, filename="histogram_DNAme.png", width = 9, height = 7) 

#############################
#############################

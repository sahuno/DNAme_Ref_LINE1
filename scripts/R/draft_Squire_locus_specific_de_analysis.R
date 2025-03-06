#squire locus specific de analysis
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
options("width"=200)


###general approach is to run 
#1: `squire_te_fwd.tsv` this is the counts file you need
#2: normalized with coding genes
#remove coding genes
rna_path <- "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA"
counts_path <- paste0(rna_path, "/CT/squire_te_fwd.tsv")
# counts_df <- read_tsv(counts_path)
counts_df <- fread(counts_path)

metadata <- "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_DNArecoded.csv"
mouseTriEpi_metadata <- read_csv(metadata)


dim(counts_df)
rownames(counts_df) <- counts_df$te.id
head(counts_df)
counts_df_2 <- counts_df %>% dplyr::select(!c("te.id", "te.name"))


workflow_dir <- "/data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/"
metadata_path <- paste0("/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_recoded.csv")
blind_transform <- TRUE #should the rlog transformation be blind
# drop_samples <- "R.S.2" #samples to drop from analysis
drop_samples <- c("R.S.2", "R.C.3") #samples to drop from analysis
# drop_samples <- NULL #samples to drop from analysis
ref_variable = "DMSO"
MIN_ReadsCounts = 50
smallestGroupSize <- 3

source(paste0(workflow_dir,"scripts/helper_functions.R"))
metadata_df <- read.csv(file = metadata_path, sep="," ,header = TRUE)

metadata_df <- metadata_df[!metadata_df$samples %in% drop_samples,]
# counts_df <- counts_df[,..metadata_df$samples]
counts_df_2 <- counts_df_2 %>% dplyr::select(!all_of(drop_samples))


# dds <- DESeqDataSetFromMatrix(countData = data.matrix(counts_df_2),
#                               colData = metadata_df,
#                               design = ~ condition)

# sapply(data.matrix(counts_df_2), is.na)
# counts_df_2_rowSums <- rowSums(counts_df_2)
# counts_df_2_rowSums_less1 <- counts_df_2_rowSums[counts_df_2_rowSums < 1.0]
# counts_df_2_rowSums_less1[counts_df_2_rowSums_less1 > 0.00]

#check if there's is line 1
#filter out L1
countsL1_df <- counts_df %>% dplyr::filter(str_detect(te.id, "L1"))
countsL1_df_withoutHLA <- counts_df %>% dplyr::filter(str_detect(te.id, "L1")) #!str_detect(te.id, "HLA")

#check length of line1
HeadcountsL1_df <- head(countsL1_df, 10)[,c(1:5)]
# HeadcountsL1_df <- countsL1_df[,c(1:5)] #uncomment to run all

# HeadcountsL1_df
HeadcountsL1_df_sep <- HeadcountsL1_df %>% separate_wider_delim(te.name, "|" ,names = c("chr", "start", "end", "name","number","strand"))
HeadcountsL1_df_sep <- HeadcountsL1_df_sep %>% separate_wider_delim(name, ":" ,names = c("consensus","clade", "class"))
L1HS_filter <- HeadcountsL1_df_sep %>% dplyr::filter(str_detect(consensus, "L1Md")) #!str_detect(te.id, "HLA")
unique(L1HS_filter$consensus)

unique(HeadcountsL1_df_sep$consensus)
HeadcountsL1_gr_sep <- makeGRangesFromDataFrame(HeadcountsL1_df_sep, keep.extra.columns = TRUE)
width(HeadcountsL1_gr_sep)
# counts_dfp
# dim(counts_df_2)
#make a matrix out of squire counts
counts_df_te.id <- counts_df$`te.id`
counts_df_matrix <- data.matrix(counts_df[,!c(1,2)])
rownames(counts_df_matrix) <- counts_df_te.id
counts_df_matrix_noZeros <- counts_df_matrix[rowSums(counts_df_matrix) > 1,]
# dim(counts_df_noZeros)
counts_df_matrix_noZeros_int <- apply(counts_df_matrix_noZeros, 2, as.integer)
rownames(counts_df_matrix_noZeros_int) <- rownames(counts_df_matrix_noZeros) #put geneIDs back

#drop low variable samples
counts_df_matrix_noZeros_int_samplesRemov <- counts_df_matrix_noZeros_int[,-which(colnames(counts_df_matrix_noZeros_int) %in% drop_samples)]


keep_TEs <- rowSums(counts_df_matrix_noZeros_int_samplesRemov >= MIN_ReadsCounts) >= smallestGroupSize
# keep <- rowSums(counts(dds)) >= 10 #keep only counts where rowSums is above 10
counts_df_matrix_noZeros_int_samplesRemov <- counts_df_matrix_noZeros_int_samplesRemov[keep_TEs,] 




################################################################################################
#load and subset protein coding datasets from rna-seq results folder - sasha
################################################################################################
# counts_annot <- fread(paste0(rna_path,'/CT/counts_annot.tsv'))
counts_annot <- read_tsv(paste0(rna_path,'/CT/counts_annot.tsv'))
counts_annot_proteinCodingOnly <- counts_annot %>% dplyr::filter(gene.type == "protein_coding")
countsOnly_proteinCodingOnly_DF <- counts_annot_proteinCodingOnly[,11 : ncol(counts_annot_proteinCodingOnly)]
# unique(counts_annot$gene.type)
countsOnly_proteinCodingOnly_DF_removeOutlierSamples <- countsOnly_proteinCodingOnly_DF %>% dplyr::select(!all_of(drop_samples))
# dim(countsOnly_proteinCodingOnly_DF_removeOutlierSamples)

keep <- rowSums(countsOnly_proteinCodingOnly_DF_removeOutlierSamples >= MIN_ReadsCounts) >= smallestGroupSize
# keep <- rowSums(counts(dds)) >= 10 #keep only counts where rowSums is above 10
countsOnly_proteinCodingOnly_DF_removeOutlierSamples <- countsOnly_proteinCodingOnly_DF_removeOutlierSamples[keep,] #filter dds
countsOnly_proteinCodingOnly_DF_removeOutlierSamples <- as.data.frame(countsOnly_proteinCodingOnly_DF_removeOutlierSamples) #tibble doesn't allow rownames so first convert to data.frame
rownames(countsOnly_proteinCodingOnly_DF_removeOutlierSamples) <- counts_annot_proteinCodingOnly$gene.id[keep]

#cbind repeats and protein codinggenes
countsRepeats_codingGenes <- rbind(countsOnly_proteinCodingOnly_DF_removeOutlierSamples, counts_df_matrix_noZeros_int_samplesRemov)

##keep dmso & Aza only
metadata_df_AZA_DMSO <- metadata_df %>% dplyr::filter(condition %in% c("DMSO", "AZA"))

#make dseq2 objects
rownames(metadata_df) <- metadata_df$samples # this is important for the DESeqDataSetFromMatrix function
# metadata_df$condition <- factor(metadata_df$condition, ref)

countsRepeats_codingGenes_AZA_DMSO <- countsRepeats_codingGenes %>% dplyr::select(all_of(metadata_df_AZA_DMSO$samples))

metadata_df_AZA_DMSO$condition <- factor(metadata_df_AZA_DMSO$condition, levels = c("DMSO", "AZA"))



dds <- DESeqDataSetFromMatrix(countData = data.matrix(countsRepeats_codingGenes_AZA_DMSO),
                              colData = metadata_df_AZA_DMSO,
                              design = ~ condition)

message("setting CTRL as reference group")
# dds$condition <- relevel(dds$condition, ref = ref_variable)

#sanity checks
all(rownames(metadata_df) == colnames(countsRepeats_codingGenes))


# keep <- rowSums(counts(dds) >= MIN_ReadsCounts) >= smallestGroupSize
# # keep <- rowSums(counts(dds)) >= 10 #keep only counts where rowSums is above 10
# dds <- dds[keep,] #filter dds
# # head(counts(dds))
# names(assays(dds))
# names(mcols(dds))

#either i do this or run step by step
dds <- DESeq(dds)


#run DESeq step by step, normalizing with coding genes; 
dds_sf <- estimateSizeFactors(dds)
dds_sf_Disp <- estimateDispersions(dds_sf)
dds_sf_Disp_TEs <- dds_sf_Disp[!grepl("^ENSMUSG", rownames(dds_sf_Disp)), ] #now remove protein coding genes for the test
dds_sf_Disp_TEs <- nbinomWaldTest(dds_sf_Disp_TEs)




#do PCA to check sample for outliers
rlogTEs <- rlog(dds_sf_Disp_TEs, blind = blind_transform)
rlogTEs_withcoding <- rlog(dds, blind = blind_transform)

pltPCA_samples <- plotPCA(rlogTEs, intgroup=c("condition", "samples"))
pltPCA_ConditionsOnly <- plotPCA(rlogTEs, intgroup=c("condition"))
ggsave(pltPCA_samples, file = paste0("figures/PCA_by_condition_and_samples_locusSpecific_SquireRepeats_normalizedWithCoding.png"), width = 9, height = 7)
ggsave(pltPCA_ConditionsOnly, file = paste0("figures/PCA_by_conditionOnly_locusSpecific_SquireRepeats_normalizedWithCoding.png"), width = 9, height = 7)


pltPCA_samples2 <- plotPCA(rlogTEs_withcoding, intgroup=c("condition", "samples"))
pltPCA_ConditionsOnly2 <- plotPCA(rlogTEs_withcoding, intgroup=c("condition"))
ggsave(pltPCA_samples2, file = paste0("figures/PCA_by_condition_and_samples_locusSpecific_SquireRepeats_normalizedWithCoding_DESeq.png"), width = 9, height = 7)
ggsave(pltPCA_ConditionsOnly2, file = paste0("figures/PCA_by_conditionOnly_locusSpecific_SquireRepeats_normalizedWithCoding_DESeq.png"), width = 9, height = 7)


# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rlogTEs)
rld_cor <- cor(rld_mat)
rld_df <- rld_mat %>% as.data.frame() %>% rownames_to_column("gene.id")



## differential expression of line1
dds_sf_Disp_L1 <- dds_sf_Disp_TEs[str_detect(rownames(dds_sf_Disp_TEs), "L1"), ] #now get LINE1 elements
dds_sf_Disp_L1 <- DESeq(dds_sf_Disp_L1, test = "LRT", reduced = ~1)
dds_sf_Disp_L1_result <- results(dds_sf_Disp_L1)



#conver to cpm 
dds_sf_Disp_L1_cpm <- edgeR::cpm(dds_sf_Disp_L1, log = TRUE, prior.count = 2) #prior.count = 2 is added to avoid log(0) error
dds_sf_Disp_L1_cpm_df <- as.data.frame(dds_sf_Disp_L1_cpm)
dds_sf_Disp_L1_cpm_df <- dds_sf_Disp_L1_cpm_df %>% rownames_to_column(var = "te.id")

dds_sf_Disp_L1_cpm_df$te.id <- gsub("\\.\\.",".",dds_sf_Disp_L1_cpm_df$te.id)

dds_sf_Disp_L1_cpm_df2 <- dds_sf_Disp_L1_cpm_df %>% separate_wider_delim(te.id, "." ,names = c("consensus","clade", "class","chrom", "start", "end", "strand")) %>% mutate(strand = case_when(strand == "minus" ~ "-", strand == "plus" ~ "+"))
# dds_sf_Disp_L1_cpm_df2 %>% mutate(strand = case_when(strand == "minus" ~ "-", strand == "plus" ~ "+"))
fwrite(dds_sf_Disp_L1_cpm_df2, file="line1_RNA_DMSO_AZA.tsv")
# fwrite(dds_sf_Disp_L1, )

table(dds_sf_Disp_L1_cpm_df2$chrom)
dds_sf_Disp_L1_cpm_df2 <- dds_sf_Disp_L1_cpm_df2 %>% filter(!str_detect(chrom, "_") )

gr_rna <- makeGRangesFromDataFrame(dds_sf_Disp_L1_cpm_df2, keep.extra.columns = TRUE)
# gr_rna <- dt2gr(dds_sf_Disp_L1_cpm_df2)


##test with bed file
df_bed_path <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/results/OverlapsL1PromoterDNAme/D-A-1_4000/fullLength/nonCGI/D-A-1_4000_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov10.bed"
df_bed <- fread(df_bed_path)

names(df_bed) <- c("chrom", "start", "end", "RepeatID", "score", "strand", "start2", "end2", "color", "min", "max", "mean", "median", "count")
gr_dna <- dt2gr(df_bed)

overlapsActiveL1_DNARNA <- gUtils::gr.findoverlaps(gr_dna, subject=gr_rna, first=FALSE,  qcol = c("RepeatID", "mean", "count"), scol=c("consensus" ,      "clade"  ,     "class" ,    "R.0.1",   "R.0.2",  "R.0.3" ,    "R.A.1" ,    "R.A.2" ,    "R.A.3"))


## now do for all bed files
#non cgi files

####################################################################################################################################
####################################################################################################################################
list_paths_overlaps <- list.files("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/results/OverlapsL1PromoterDNAme/", recursive = TRUE, full.names = TRUE, pattern=".bed")

# nCgi <- 
paths_100bp5UTR_CGI <- list_paths_overlaps[str_detect(list_paths_overlaps,"100bp5UTR/CpGIs")]
paths_100bp5UTR_nonCGI <- list_paths_overlaps[str_detect(list_paths_overlaps,"100bp5UTR/nonCGI")]

paths_400600bp5UTR_CpGIs <- list_paths_overlaps[str_detect(list_paths_overlaps,"400600bp5UTR/CpGIs")]
paths_400600bp5UTR_nonCpGIs <- list_paths_overlaps[str_detect(list_paths_overlaps,"400600bp5UTR/nonCGI")]

paths_fullLength_nonCpGIs <- list_paths_overlaps[str_detect(list_paths_overlaps,"fullLength/nonCGI")]
paths_fullLength_CpGIs <- list_paths_overlaps[str_detect(list_paths_overlaps,"fullLength/CpGIs")]

l1_path <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed"

l1_entireLength <- fread(l1_path, col.names = c("chrom", "start", "end", "RepeatID", "score", "strand", "start2", "end2", "color"))
l1_keys <- gUtils::gr.findoverlaps(makeGRangesFromDataFrame(l1_entireLength, keep.extra.columns = TRUE), subject=gr_rna, first=FALSE,  qcol = c("RepeatID"), scol=c("consensus" ,      "clade"  ,     "class" ,    "R.0.1",   "R.0.2",  "R.0.3" ,    "R.A.1" ,    "R.A.2" ,    "R.A.3"), return.type = "data.table")
l1_keys <- l1_keys %>% pivot_longer(cols = starts_with("R."), names_to="RNAsamples", values_to = "cpm" ) %>% mutate(samples = str_replace_all(RNAsamples,"\\.","-"), samples = str_replace_all(samples,"R","D"))


# names(df_bed) <- 
DNAmeOverlaps <- lapply(paths_fullLength_nonCpGIs, function(x){fread(x, col.names = c("chrom", "start", "end", "RepeatID", "score", "strand", "start2", "end2", "color", "min", "max", "mean", "median", "count"))})
names(DNAmeOverlaps) <- basename(paths_fullLength_nonCpGIs)


# gUtils::gr.findoverlaps(makeGRangesFromDataFrame(DNAmeOverlaps[[1]], keep.extra.columns = TRUE), subject=gr_rna, first=FALSE,  qcol = c("RepeatID", "mean", "count"), scol=c("consensus" ,      "clade"  ,     "class" ,    "R.0.1",   "R.0.2",  "R.0.3" ,    "R.A.1" ,    "R.A.2" ,    "R.A.3"))

OverlapsDNAmeRNA <- lapply(DNAmeOverlaps, function(x){
    ovs <- gUtils::gr.findoverlaps(gUtils::dt2gr(x), subject=gr_rna, first=FALSE,  qcol = c("RepeatID", "mean", "count"), scol=c("consensus" ,      "clade"  ,     "class" ), return.type = "data.table")
})


# OverlapsDNAmeRNA <- lapply(DNAmeOverlaps, function(x){
#     ovs <- gUtils::gr.findoverlaps(gUtils::dt2gr(x), subject=gr_rna, first=FALSE,  qcol = c("RepeatID", "mean", "count"), scol=c("consensus" ,      "clade"  ,     "class" ,    "R.0.1",   "R.0.2",  "R.0.3" ,    "R.A.1" ,    "R.A.2" ,    "R.A.3"), return.type = "data.table")
# })



# lapply(OverlapsDNAmeRNA, length)
# DNAmeOverlaps_df <- rbindlist(DNAmeOverlaps, idcol = "sample")

OverlapsDNAmeRNA_df <- rbindlist(OverlapsDNAmeRNA, idcol = "samples")

OverlapsDNAmeRNA_df[, minCov:=(str_extract(samples, "minCov[0-9]*"))]
# OverlapsDNAmeRNA_df[, samples:=(str_extract(samples, ".*(?=_5mCpG)"))]
OverlapsDNAmeRNA_df[, samples:=(str_extract(samples, "^[^_]+"))]
OverlapsDNAmeRNA_df_joined <- OverlapsDNAmeRNA_df %>% left_join(l1_keys, by = c("samples", "RepeatID"))

OverlapsDNAmeRNA_df_joined <- OverlapsDNAmeRNA_df_joined %>% left_join(mouseTriEpi_metadata)
 
OverlapsDNAmeRNA_df_joined %>% filter(is.na(cpm))
#so now remove anything aside these samples [A-Z]-[0|A] - [1-3]
DMSO_AZaL1DNAmeCPM <- OverlapsDNAmeRNA_df_joined %>% filter(str_detect(samples, "D-0|D-A") & minCov == "minCov5")

DMSO_AZaL1DNAmeCPM_sorted <- DMSO_AZaL1DNAmeCPM %>% group_by(condition) %>% arrange(cpm, .by_group = TRUE) %>% ungroup()
# unique(DMSO_AZaL1DNAmeCPM_sorted$RepeatID)
# DMSO_L1DNAmeCPM <- DMSO_AZaL1DNAmeCPM %>% filter(str_detect(samples, "D-0") & minCov == "minCov5")

DMSO_AZaL1DNAmeCPM_sorted$mean <- as.numeric(DMSO_AZaL1DNAmeCPM_sorted$mean)
DMSO_AZaL1DNAmeCPM_sorted$RepeatID <- paste0(DMSO_AZaL1DNAmeCPM_sorted$`seqnames.x`, ":",DMSO_AZaL1DNAmeCPM_sorted$`start.x`,"_", DMSO_AZaL1DNAmeCPM_sorted$`end.x`, ":", DMSO_AZaL1DNAmeCPM_sorted$RepeatID, ":",DMSO_AZaL1DNAmeCPM_sorted$`consensus.x`)


RepeatIDGroups <- split(DMSO_AZaL1DNAmeCPM_sorted, as.factor(DMSO_AZaL1DNAmeCPM_sorted$RepeatID))
# names(RepeatIDGroups)

plotL1perCondition <- function(l1Name){
plotOut <- ggplot(RepeatIDGroups[[l1Name]], aes(x=mean, y=cpm, color=condition)) + geom_point() + 
geom_text_repel(aes(label=samples)) + labs(title = paste(l1Name), y = "cpm", x = "mean DNAme") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
}


plotL1SingleCondition <- function(l1Name){
plotOut <- ggplot(RepeatIDGroups[[l1Name]] %>% dplyr::filter(condition=="DMSO"), aes(x=mean, y=cpm)) + geom_point() + 
geom_text_repel(aes(label=samples)) + labs(title = paste(l1Name), y = "cpm", x = "mean DNAme") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
}

# dir.create("results/boxplotsRepPerCond/")
pdf(paste0("scatter_DMSOAZA_lspecific_L1RNA_DNAme",".pdf"))
lapply(names(RepeatIDGroups), plotL1perCondition)
dev.off()

pdf(paste0("scatter_DMSO_lspecific_L1RNA_DNAme",".pdf"))
lapply(names(RepeatIDGroups), plotL1SingleCondition)
dev.off()
# sum(is.na(OverlapsDNAmeRNA_df_joined$cpm))
# OverlapsDNAmeRNA_df[, samples:=(str_replace_all(samples, ".", ""))]


##########################
# ### get cordinates
# TEcounts <- fread("/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/mapping/R-0-1/squire/aln_TEcounts.txt.gz")
# # make a bed file of TE cordinates from the TEcounts 
# namesTEs <- c("TE_chr",   "TE_start",   "TE_stop" ,  "TE_name" ,   "TE_strand",  "milliDiv" ,  "tx_chr",      "tx_start",    "tx_stop"  ,   "TE_ID"   ,    "fpkm" ,    "tx_strand"  , "Sample"   ,   "alignedsize","uniq_counts", "tot_counts" , "tot_reads",   "score")
# TEcounts_ColnamesRearranged <- TEcounts[, ..namesTEs]

# setnames(TEcounts_ColnamesRearranged, c("TE_chr",   "TE_start",   "TE_stop" ,  "TE_name" ,   "TE_strand"), c("chr",   "start",   "end" ,  "TE_name" ,   "strand"))

# fwrite(TEcounts_ColnamesRearranged, "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/R01_TEcounts_DNA_RNA.bed", sep="\t", quote = FALSE, row.names = FALSE)
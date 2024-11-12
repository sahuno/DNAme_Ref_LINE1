




#get locus specific L1 rna count data
counts_df <- fread("/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/CT/squire_te_fwd.tsv")

# data.frame (gene.id=NA,ensembl.gene.id=NA, chr=NA,   gene.symbol=NA, gene.type=NA, length=NA description=NA entrez.gene.id uniprot.id go.gene.id)

#filter out L1
counts_df_L1 <- counts_df %>% dplyr::filter(str_detect(te.id, "L1")) #!str_detect(te.id, "HLA")

#check length of line1
# HeadcountsL1_df <- head(counts_df, 10)[,c(1:8)]
# HeadcountsL1_df <- countsL1_df[,c(1:5)] #uncomment to run all

counts_df_L1 <- counts_df_L1 %>% separate_wider_delim(te.name, "|" ,names = c("chrom", "start", "end", "name","number","strand"))
counts_df_L1 <- counts_df_L1 %>% separate_wider_delim(name, ":" ,names = c("consensus","clade", "class"))

#turn data.table to granges
l1RNA_gr <- makeGRangesFromDataFrame(HeadcountsL1_df_sep, keep.extra.columns = TRUE)

# read full length line1 from L1Base
l1_path <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed"
l1_entireLength <- fread(l1_path, col.names = c("chrom", "start", "end", "RepeatID", "score", "strand", "start2", "end2", "color"))


lsRNA_keys <- gUtils::gr.findoverlaps(makeGRangesFromDataFrame(l1_entireLength, keep.extra.columns = TRUE), subject=l1RNA_gr, first=FALSE,  qcol = c("RepeatID"), scol=c("te.id", "name",  "number",  sort(paste0(rep("R.", 6),rep(c("0.","A.")) ,c(1,2,3))) ), return.type = "data.table")

lsRNA_keys %>% filter(!str_detect(seqnames, "chrY|chrM|chrX"))
#filter locus specific l1 with enough reads
keep_TEs <- rowSums(lsRNA_keys[,c(11:16)] >= 10) >= 3

lsRNA_keys[keep_TEs, ]


### read overlaps data
list_paths_overlaps <- list.files("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/results/OverlapsL1PromoterDNAme/", recursive = TRUE, full.names = TRUE, pattern=".bed")
paths_fullLength_nonCpGIs <- list_paths_overlaps[str_detect(list_paths_overlaps,"fullLength/nonCGI")]
paths_fullLength_nonCpGIs <- paths_fullLength_nonCpGIs[str_detect(paths_fullLength_nonCpGIs,"_minCov5")]


DNAmeRatesOverlaps <- lapply(paths_fullLength_CpGIs, function(x){fread(x, col.names = c("chrom", "start", "end", "RepeatID", "score", "strand", "start2", "end2", "color", "min", "max", "mean", "median", "count"))})
names(DNAmeRatesOverlapsv) <- basename(paths_fullLength_CpGIs)
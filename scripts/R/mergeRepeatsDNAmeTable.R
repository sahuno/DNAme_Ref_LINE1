# merge all the dNAme Stats
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/mergeRepeatsDNAmeTable.R

# load libraries
library(data.table)

dir="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/RepeatsDNAmeStats"

files <- list.files(dir, pattern = "_repeatsDNAme.tsv", recursive=TRUE,full.names = TRUE)
#[1:3]

repeatDT_list <- lapply(files, function(x){
    repeatDT <- fread(x)
    sampleName <- basename(x)
    return(repeatDT)
})

names(repeatDT_list) <- gsub("_repeatsDNAme.tsv", "", basename(files))

# cols <- c(
#   "ID", "fracMethyl.N", "fracMethyl.mean", "fracMethyl.median", 
#   "fracMethyl.entropy_in_log2", "fracMethyl.geom_mean", "seqnames", "start", 
#   "end", "strand", "width", "name", "number", "consensus", "clade", 
#   "class", "validReadCounts", "TE_ID"
# )

cols <- c(
  "ID", "TE_ID","fracMethyl.N", "fracMethyl.mean", "fracMethyl.median", 
  "fracMethyl.entropy_in_log2", "fracMethyl.geom_mean"
)

# metadataCols <- c(
#   "ID", "seqnames", "start", 
#   "end", "strand", "width", "name", "number", "consensus", "clade", 
#   "class", "TE_ID", "validReadCounts")

metadataCols <- c(
  "ID", "chr", "start", 
  "end", "strand", "name", "number", "consensus", "clade", 
  "class", "TE_ID", "validReadCounts")


metaData <- repeatDT_list[[1]][,..metadataCols]

repeatStatsDT_list <- lapply(names(repeatDT_list), function(x){
    repeatDT_list[[x]][,..cols]
})

names(repeatStatsDT_list) <- gsub("_repeatsDNAme.tsv","",basename(files))


repeatStatsDT <- rbindlist(repeatStatsDT_list, idcol = "sample")

repeatStatsDT_plusMeta <- repeatStatsDT[unique(metaData), on = c("TE_ID")]

## write these to disk
fwrite(repeatStatsDT_plusMeta, sep="\t", file = paste0("RepeatsDNAmeStatsAllSamples_ValidnNonValidReadCounts_Aggregated_FWD_REV.bed"))

# sum(duplicated(unique(metaData)$TE_ID))

# unique(metaData)
# unique(metaData)[TE_ID=="chr1|3014751|3021011|L1Md_F2:L1:LINE|84|-",]

# grep -E "chr1|3014751|3021011|L1Md_F2:L1:LINE|84|-" /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/RepeatsDNAmeStats/D-0-3_4000_repeatsDNAme.tsv | less

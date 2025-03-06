#author: samuel ahuno
#date: march 3rd 2025
#purpose: assign name to full length L1

#Note:
# run this script first
# $ /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/mergeSquireRNA_cordinates_master.sh

library(data.table)
#genomics utils package from marcin imel.
library(gUtils)

#options
path_file <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/overlapsSquireL1.bed"
bedfile_names <- c("seqnames",    "start",      "end", "rm_name", "rm_score", "strand", "start", "end", "L1UID_itemRgb")
analysis_key = "mm10_mmflil1_8438" #change


df <- fread(path_file)
#df <- read.table(path_file)


#assign names to colmns
oldColumns <- c("V1","V2","V3","V4","V5","V6", "V7", "V8","V9", "V10","V11", "V12", "V13",  "V14", "V15")
#newColumns <- c("seqnames", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "rmSeqname", "rmStart", "rmEnd", "rmName", "rmScore", "rmStrand")
#newColumns <- c("L1UID_seqnames", "L1UID_start", "L1UID_end", "L1UID_name", "L1UID_score", "L1UID_strand", "L1UID_thickStart", "L1UID_thickEnd", "L1UID_itemRgb", "rmSeqname", "rmStart", "rmEnd", "rmName", "rmScore", "rmStrand")
newColumns <- c("L1UID_seqnames", "L1UID_start", "L1UID_end", "L1UID_name", "L1UID_score", "L1UID_strand", "L1UID_thickStart", "L1UID_thickEnd", "L1UID_itemRgb", "seqnames", "start", "end", "rm_name", "rm_score", "strand")

data.table::setnames(df, old=oldColumns, new=newColumns)

#data.table::setnames(df, old=c("V1","V2","V3","V4","V5","V6", "V7", "V8","V9", "V10","V11", "V12", "V13", "V14"), new=c("seqnames", "start", "end", "name", "color", "strand","score", "thickStart", "thickEnd", "itemRgb", "rmSeqname", "rmStart", "rmEnd", "rmName","rmcolor", "rmStrand"))
head(df)
gr <- dt2gr(df)
mcols(gr)$width_rmL1 <- width(gr)


#### randomly sample, write to bed to visualize in IGV
dt_2save <- gr2dt(gr)
dtUID100 <- dt_2save[L1UID_name == "UID-100", ..bedfile_names]
dtUID2 <- dt_2save[L1UID_name == "UID-2", ..bedfile_names]
dtUID2500 <- dt_2save[L1UID_name == "UID-2500", ..bedfile_names]

fwrite(dtUID100, file="mapped_repeatMasker_L1Base_UID100_mm10.bed", sep="\t", col.names = FALSE)
fwrite(dtUID2, file="mapped_repeatMasker_L1Base_UID2_mm10.bed", sep="\t",col.names = FALSE)
fwrite(dtUID2500, file="mapped_repeatMasker_L1Base_UID2500_mm10.bed", sep="\t", col.names = FALSE)

# chr1_148378685_148384680_UID100
# chr1_5816878_5827110_UID2_mm10
# chr7_17241390_17251787_UID2500_mm10
grl <- lapply(split(gr, gr$L1UID_name), function(x) {sort(x)})

#how many repeat masker l1 are in L1Base cordianates
grl_lens <- lapply(grl, length)
# grl[['UID-100']]
# sort(grl[['UID-100']])

#sum should be around 10kb
sum(grl[['UID-100']]$width_rmL1)

# data.frame(grl_lens)
# rbindlist(grl_lens)
#plot rmasker per uid

#compute pair-wise distance
grl_pairwiseDistance <- lapply(grl, function(x){gr.dist(x)})
grl_pairwiseDistance[['UID-100']]
# grl_pairwiseDistance[['UID-100']]

################################################
## which is longest repeat
################################################
# grl[['UID-100']]$width_rmL1[which.max(grl[['UID-100']]$width_rmL1)]
grl_dtl <- lapply(grl, function(x){gr2dt(x)[which.max(x$width_rmL1),]})
grl_dt <- rbindlist(grl_dtl)

message("writing data to disk")
newBedColumnNames <- c(bedfile_names, "width_rmL1", "L1UID_name", "L1UID_strand")
fwrite(grl_dt, file="mapped_repeatMasker_L1Base.tsv", sep="\t")
fwrite(grl_dt[,..newBedColumnNames], file=paste0("mapped_repeatMasker_AssignedTo_L1Base_",analysis_key,".tsv"), sep="\t") #for filtering purposes 


# analysis_key
# grlgrl_dt
#grl <- GenomicRanges::makeGRangesListFromDataFrame(df, split.field = "L1UID_name", keep.extra.columns = TRUE)

#width(grl)#
# x$width_rmL1 <-
#grl_lens <- lapply(grl, function(x) {width(x)})
#lapply(seq_along(grl), function(x) {mcols(grl[[x]])$width_rmL1 <- width(grl[[x]])})



#grList <- split(gr, gr$name)
#width(gr) <- 


# cat squireL1.bed | cut -f 1-6 > squireL1_cordinatesOnly.bed
# mergeSquireRNA_cordinates_master.sh
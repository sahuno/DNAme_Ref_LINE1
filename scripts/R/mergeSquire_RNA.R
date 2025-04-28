#merge Squire RNA results
#TO RUN:
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/mergeSquire_RNA.R

#load libraries
library(data.table)
library(fst)
library(S4Vectors)
library(GenomicRanges)
library(annotatr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
# library(gUtils)

minReadCounts <- 10

mergeSquireL1 <- "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/mapping/"
#mergeSquireL1 <- "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/mapping/R-0-2/squire/aln_TEcounts.txt.gz

pathsDT <- list.files(mergeSquireL1, recursive = TRUE, full.names = TRUE, pattern = "aln_TEcounts.txt.gz")
# pathsDT <- pathsDT[1:3]

dir.create("mergedSquireRepRNA")
dir.create("mergedSquireRepRNA/perSample", recursive = TRUE)

###################################
############### test
###################################
computeCols = c("tot_counts", "fpkm")
colsInterested <- c("TE_chr", "TE_start" , "TE_stop", "TE_ID") #"TE_name"
colsInterestedFinal <- c("TE_chr", "TE_start" , "TE_stop", "TE_ID", "tot_counts", "fpkm") 
colBed <- c("chr", "start", "end", "name", "number" , "strand", "consensus", "clade",  "class", "validReadCounts", "TE_ID")
colAnno <- c("chr", "start", "end", "name", "number" , "strand", "consensus", "clade",  "class","TE_ID")


# dt1 <- fread("/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/mapping/R-0-2/squire/aln_TEcounts.txt.gz")
# dt1_Stats <- dt1[, lapply(.SD, function(x){sum(x)}),by = TE_ID, .SDcols = computeCols]

# dt1_Stats_Annot <- dt1_Stats[unique(dt1[,..colsInterested]), on="TE_ID"]
# dt1_Stats_Annot <- dt1_Stats_Annot[, ..colsInterestedFinal]
# dt1_Stats_Annot <- dt1_Stats_Annot[order(TE_chr, TE_start)]
#########################################################################################################
#########################################################################################################


#sanity checks!
# dt1[TE_ID == "chr12|69159294|69159594|7SLRNA:srpRNA:srpRNA|20|+",]
# dt1[TE_ID == "chrY_JH584302_random|115914|118440|L1Md_T:L1:LINE|97|+",]
# dt1[grepl("chrY|90016996|90017986|L1Md_F2:L1:LINE|77", TE_ID),]

SQuireRepList <- lapply(pathsDT, function(x){
    RepeatsDT <- fread(x)
    sampleName <- stringr::str_split_i(x, "/", -3)
    RepeatsStatsDT <- RepeatsDT[, lapply(.SD, function(x){sum(x)}),by = TE_ID, .SDcols = computeCols]
    RepeatsAnnotDT <- RepeatsStatsDT[unique(RepeatsDT[,..colsInterested]), on="TE_ID"]
    RepeatsAnnotDT <- RepeatsAnnotDT[order(TE_chr, TE_start), ..colsInterestedFinal] #sort and order
    RepeatsAnnotDT[,c("chr", "start", "end", "name","number","strand") := tstrsplit(`TE_ID`, "|", fixed = TRUE)][, c("consensus","clade", "class"):=tstrsplit(name, ":", fixed = TRUE)] ##split the `TE_ID` column
    write_fst(RepeatsAnnotDT, path = paste0("mergedSquireRepRNA/perSample/",sampleName, "_squireRepeatsAnnot_Aggregated_FWD_REV.fst"))
})

#merge data for all samples
sampleNames <- stringr::str_split_i(pathsDT, "/", -3)
names(SQuireRepList) <- sampleNames

SQuireRepDT <- rbindlist(SQuireRepList, idcol= "sample")
SQuireRepDT <- SQuireRepDT[!grep("\\?",class),]
SQuireRepDT <- SQuireRepDT[!grep("\\?",clade),] #remove and not sure annotations


# table(SQuireRepDT$chr)
SQuireRepDT <- SQuireRepDT[!grepl("Un|_", chr)]
## needed for annotation
uniqueSQuireRepDT <- unique(SQuireRepDT[, ..colAnno])


# setnames(SQuireRepDT, "chr", "seqnames")
setnames(uniqueSQuireRepDT, "strand", "ignoredStrand")

# SQuireRepGr <- gUtils::dt2gr(SQuireRepDT)

uniqueSQuireRepDF <- as.data.frame(uniqueSQuireRepDT)
# dim(uniqueSQuireRepDF)
uniqueSQuireRepGr <- makeGRangesFromDataFrame(head(uniqueSQuireRepDF, 10000), keep.extra.columns = TRUE, seqnames.field="chr", start.field="start", end.field="end", starts.in.df.are.0based = TRUE, ignore.strand = TRUE)

head(uniqueSQuireRepDF)
dim(unique(head(uniqueSQuireRepDF, 10000)))

# 2. See which mm10 annotation sets are available
all_annots <- builtin_annotations()
mouse_annots <- all_annots[grepl("^mm10_", all_annots)]
print(mouse_annots)
# 3. Build a subset of annotations
annots_mm10 <- build_annotations(
  genome     = "mm10",
  annotations = c(
    # "mm10_genes_promoters",
    "mm10_genes_exons",
    "mm10_genes_introns",
    "mm10_genes_intergenic",
    "mm10_genes_5UTRs",
    "mm10_genes_3UTRs"
    )
)

message("annotations beginning")
annotated_mm10 <- annotate_regions(
  regions     = uniqueSQuireRepGr,
  annotations = annots_mm10,
  ignore.strand = TRUE
)

#convert granges to Data.table
annotated_mm10DT <- gUtils::gr2dt(annotated_mm10)#[, c("seqnames","start","end","annot.type")]
# annotated_mm10DT
setnames(annotated_mm10DT, "seqnames", "chr")
########################
#### discovery! not all regions are annotated
#annotated_mm10DT[TE_ID == "chr1|4606627|4607142|L1Md_T:L1:LINE|12|+",]
########################################
# unique(annotated_mm10DT)

annotated_mm10DT[
  ,
  `:=`(
    annot = fcase(
      grepl("intergenic", annot.type),                                 "intergenic",
      grepl("introns|promoters|exons|5UTRs|3UTRs", annot.type),        "intragenic",
      default = NA_character_
    ),
    geneID_n_Symbol = paste0(annot.gene_id, "|", annot.symbol)
  )
]
# annotated_mm10DT[TE_ID == "chr1|4606627|4607142|L1Md_T:L1:LINE|12|+",]

# ## bind the `annot.type` column to always refer back to the original
# annotated_mm10DT[
#   , annotTypeBinded := paste(annot.type, collapse = ";")
#   , by = TE_ID
# ]


# annotated_mm10DT[TE_ID == "chr1|58466514|58466621|U6:snRNA:snRNA|56|-",]

# annotated_mm10DT[str_detect("\\;", annotTypeBinded),]
# ## 
# annotated_mm10DT[
#   , annotTx_idBinded := paste(annot.tx_id, collapse = ";")
#   , by = TE_ID
# ]

# 3) drop the unwanted `annot.*` columns
drop_cols <- c(
  "annot.seqnames", "annot.start", "annot.end", "annot.width", "annot.strand",
  "annot.id",       "annot.tx_id",  "annot.gene_id", "annot.symbol", "annot.type", "strand", "width" 
)

final_classification <- annotated_mm10DT[
  , setdiff(names(annotated_mm10DT), drop_cols), with = FALSE
]

dim(unique(final_classification))

final_classificationUniq <- unique(final_classification)
# final_classificationUniq[, !c(annot, geneID_n_Symbol)]
# final_classificationUniq[, c("annot", "geneID_n_Symbol") := NULL]

dim(unique(final_classificationUniq))

# dim(final_classificationUniq[, !c("annot", "geneID_n_Symbol"), with = FALSE])

final_classificationUniq[
  , geneID_n_Symbol_collapsed := paste(unique(geneID_n_Symbol), collapse = ";")
  , by = TE_ID
][,geneID_n_Symbol := NULL]

final_classificationUniq2 <- unique(final_classificationUniq)

# head(uniqueSQuireRepDF, 1000)$TE_ID
# final_classificationUniq2$TE_ID

missing_TE_IDs <- setdiff(uniqueSQuireRepDF$TE_ID, final_classificationUniq2$TE_ID)
# final_classificationUniq2

missingDF <- uniqueSQuireRepDF[uniqueSQuireRepDF$TE_ID %in% missing_TE_IDs, ]
setDT(missingDF)
missingDF[ , `:=`(
  annot = NA_character_,
  geneID_n_Symbol_collapsed = NA_character_
)]



final_classificationWithFailedAnnotations <- rbind(missingDF, final_classificationUniq2)
dim(final_classificationWithFailedAnnotations)
# head(uniqueSQuireRepDF, 1000)[TE_ID %in% missing_TE_IDs,]

# print(missing_TE_IDs_dt)


##sanity checks
# final_classificationUniq2 %>% filter(TE_ID %in% c("chr1|58466514|58466621|U6:snRNA:snRNA|56|-", "chr1|4778130|4778258|RSINE1:B4:SINE|339|-", "chr1|86935678|86936896|L1Md_F2:L1:LINE|58|+")) 








final_classification <- annotated_mm10DT %>%
    mutate(annot= case_when(str_detect(annot.type, "intergenic") ~ "intergenic", 
                            str_detect(annot.type, "introns|promoters|exons|5UTRs|3UTRs") ~ "intragenic",
                            TRUE ~ NA_character_),
                            geneID_n_Symbol = paste0(`annot.gene_id`,"|" ,`annot.symbol`)) %>% dplyr::select(!c("annot.seqnames", annot.start, annot.end, annot.width, annot.strand, annot.id, annot.tx_id, annot.gene_id, annot.symbol,annot.type))

message("dim of unique classification")
dim(unique(final_classification))

final_classificationV2 <- unique(final_classification)
final_classificationV2
# sum(grepl("NA|NA", final_classificationV2$geneID_n_Symbol))
final_classificationV2 %>% group_by(TE_ID) %>% summarise(n = n()) %>% filter(n > 1) %>% arrange(desc(n))

final_classificationV2 %>% filter(duplicated(TE_ID) | duplicated(TE_ID, fromLast = TRUE))
final_classificationV2 %>% filter(TE_ID == "chr1|4778130|4778258|RSINE1:B4:SINE|339|-") 
final_classificationV2 %>% filter(TE_ID == "chr1|58466514|58466621|U6:snRNA:snRNA|56|-") 

testRowwise <- final_classificationV2 %>% filter(TE_ID %in% c("chr1|58466514|58466621|U6:snRNA:snRNA|56|-", "chr1|4778130|4778258|RSINE1:B4:SINE|339|-", "chr1|86935678|86936896|L1Md_F2:L1:LINE|58|+")) 


setDT(testRowwise)

testRowwise[
  , geneID_n_Symbol_collapsed := paste(unique(geneID_n_Symbol), collapse = ";")
  , by = TE_ID
][,geneID_n_Symbol := NULL]
unique(testRowwise)

testRowwise[,geneID_n_Symbol := NULL]

testRowwise %>% row_wise() %>% mutate(GeneAnnoFinal = paste0(TE_ID, "|", geneID_n_Symbol)) %>% select(GeneAnnoFinal) %>% ungroup() %>% distinct() %>% arrange(GeneAnnoFinal)

final_classificationV2 %>% summarise(
    n = n(),
    n_intergenic = sum(annot == "intergenic"),
    n_intragenic = sum(annot == "intragenic"),
    n_annotated = sum(!is.na(annot))
)

SQuireRepeatsAnnotationsDT <- unique(final_classification)
# head(as.data.frame(annotated_mm10)[, c("seqnames","start","end","annot.type")])
SQuireRepeatsAnnotationsDT[,geneID_n_Symbol:=NULL]
unique(final_classification[, -c("geneID_n_Symbol")])


# GenomicRanges::
# SQuireRepDF <- as.data.table(plyr::ldply(SQuireRepList, data.table, .id = "sample"))


#round to be used by DeSeq2
SQuireRepeatsAnnotationsDT[, tot_counts := round(tot_counts)]
# SQuireRepeatsAnnotationsDT[TE_ID == "chr1|4773954|4774028|MIRc:MIR:SINE|192|-",]
SQuireRepeatsAnnotationsDT[, validReadCounts:=(sum(tot_counts) > minReadCounts), by = c("TE_ID")]


setcolorder(SQuireRepeatsAnnotationsDT, c(setdiff(names(SQuireRepeatsAnnotationsDT), "sample"), "sample")) #need for bedFile

# SQuireRepeatsAnnotationsDT[validReadCounts == TRUE]
fwrite(SQuireRepeatsAnnotationsDT, sep="\t", file = paste0("mergedSquireRepRNA/squireRepeats_AggregatedFWDnREV_masterTable.bed"), col.names = FALSE)
write_fst(SQuireRepeatsAnnotationsDT, path = paste0("mergedSquireRepRNA/squireRepeats_AggregatedFWDnREV_masterTable.fst"), compress = 100)


# SQuireRepeatsAnnotationsDT <- read_fst(path = paste0("squireRepeats_AggregatedFWDnREV_masterTable.fst"), as.data.table = TRUE)
# SQuireRepeatsAnnotationsDT[TE_ID == "chr1|4773954|4774028|MIRc:MIR:SINE|192|-",]

# 
##need bedlike format 
SQuireRepeatsAnnotationsDT2Bed <- unique(SQuireRepeatsAnnotationsDT[validReadCounts == TRUE, ..colBed])

fpkm_totalCounts <- unique(SQuireRepeatsAnnotationsDT[validReadCounts == TRUE, c("TE_ID", "tot_counts",  "fpkm", "sample")])
tot_countsdt <- dcast(fpkm_totalCounts[,.SD, .SDcols=!c("fpkm")], ... ~ sample, value.var = "tot_counts", fill = 0)

#use this instead for line1Overlaps
SQuireRepeatsAnnotationsDT2BedCounts <- SQuireRepeatsAnnotationsDT2Bed[tot_countsdt, on="TE_ID"]
fwrite(SQuireRepeatsAnnotationsDT2BedCounts, sep="\t", file = paste0("mergedSquireRepRNA/SQuireRepeatsValidReadCounts_Aggregated_FWD_REV_totCounts.bed"), col.names = FALSE)
# fwrite(SQuireRepeatsAnnotationsDT2BedCounts[class == "LINE",], sep="\t", file = paste0("mergedSquireRepRNA/SQuireL1ValidReadCounts_Aggregated_FWD_REV_totCounts.bed"), col.names = FALSE)

# unique(SQuireRepeatsAnnotationsDT[,..colBed])[TE_ID == "chr1|4773954|4774028|MIRc:MIR:SINE|192|-",]
# fwrite(SQuireRepeatsAnnotationsDT2Bed[class == "LINE",], sep="\t", file = paste0("SQuireL1ValidReadCounts_Aggregated_FWD_REV.bed"), col.names = FALSE) #use this for overlaps


fwrite(SQuireRepeatsAnnotationsDT2Bed, sep="\t", file = paste0("mergedSquireRepRNA/SQuireRepeatsValidReadCounts_Aggregated_FWD_REV.bed"), col.names = FALSE)
fwrite(SQuireRepeatsAnnotationsDT2Bed[class == "LINE",], sep="\t", file = paste0("mergedSquireRepRNA/SQuireLINEsValidReadCounts_Aggregated_FWD_REV.bed"), col.names = FALSE) #use this for overlaps
fwrite(SQuireRepeatsAnnotationsDT2Bed[class == "LINE" & clade == "L1",], sep="\t", file = paste0("mergedSquireRepRNA/SQuireL1ValidReadCounts_Aggregated_FWD_REV.bed"), col.names = FALSE) #use this for overlaps
write_fst(SQuireRepeatsAnnotationsDT2Bed[class == "LINE",], path = paste0("mergedSquireRepRNA/SQuireLINEsValidReadCounts_Aggregated_FWD_REV.fst"))
write_fst(SQuireRepeatsAnnotationsDT2Bed[class == "LINE" & clade == "L1",], path = paste0("mergedSquireRepRNA/SQuireLINE1ValidReadCounts_Aggregated_FWD_REV.fst"))
write_fst(unique(SQuireRepeatsAnnotationsDT[class == "LINE" & clade == "L1",..colBed]), path = paste0("mergedSquireRepRNA/SQuireL1ValidnNonValidReadCounts_Aggregated_FWD_REV.fst")) #needed for RNA+MethylationSummary master table

fwrite(unique(SQuireRepeatsAnnotationsDT[class == "LINE" & clade == "L1",..colBed]), file = paste0("mergedSquireRepRNA/SQuireL1ValidnNonValidReadCounts_Aggregated_FWD_REV.bed"), col.names = FALSE, sep = "\t") #needed for RNA+MethylationSummary master table
fwrite(unique(SQuireRepeatsAnnotationsDT[class == "LINE" & clade == "L1" & consensus %in% c("L1Md_A", "L1Md_T", "L1Md_F", "L1Md_F2", "L1Md_F3", "L1Md_Gf"),..colBed]), file = paste0("mergedSquireRepRNA/SQuireL1MostActive_ValidnNonValidReadCounts_Aggregated_FWD_REV.bed"), col.names = FALSE, sep = "\t") #needed for minimal methlylation analysi





fwrite(unique(SQuireRepeatsAnnotationsDT[,..colBed]), sep="\t", file = paste0("mergedSquireRepRNA/SQuireRepeatsValidnNonValidReadCounts_Aggregated_FWD_REV.bed"), col.names = FALSE) #needed for rna+methylation master table
write_fst(unique(SQuireRepeatsAnnotationsDT[,..colBed]), path = paste0("mergedSquireRepRNA/SQuireRepeatsValidnNonValidReadCounts_Aggregated_FWD_REV.fst")) #needed for RNA+MethylationSummary master table


# unique(SQuireRepeatsAnnotationsDT[,..colBed][TE_ID == "chr1|4773954|4774028|MIRc:MIR:SINE|192|-",])


write_fst(SQuireRepeatsAnnotationsDT2Bed, path = paste0("mergedSquireRepRNA/SQuireRepeatsValidReadCounts_Aggregated_FWD_REV.fst"))

# SQuireRepeatsAnnotationsDT[,c("chr", "start", "end", "name","number","strand") := tstrsplit(`TE_ID`, "|", fixed = TRUE)][, c("consensus","clade", "class"):=tstrsplit(name, ":", fixed = TRUE)] ##split the te.name column



# dcast(SQuireRepeatsAnnotationsDT[,.SD, .SDcols=!c("fpkm")], TE_ID ~ tot_counts, value.var = "tot_counts")
wide_DT_tot_counts <- dcast(SQuireRepeatsAnnotationsDT[,.SD, .SDcols=!c("fpkm")], ... ~ sample, value.var = "tot_counts", fill = 0)
wide_DT_fpkm <- dcast(SQuireRepeatsAnnotationsDT[,.SD, .SDcols=!c("tot_counts")], TE_chr + TE_start + TE_stop + TE_ID + validReadCounts ~ sample, value.var = "fpkm", fill = 0)

write_fst(wide_DT_tot_counts, path = paste0("mergedSquireRepRNA/squireRepeatsAnnot_Aggregated_FWD_REV_totalCounts.fst"))
write_fst(wide_DT_fpkm, path = paste0("mergedSquireRepRNA/squireRepeatsAnnot_Aggregated_FWD_REV_fpkm.fst"))

# mergedSquireRepRNA/SQuireRepeatsValidReadCounts_Aggregated_FWD_REV_totCounts.bed
cmd = "bedtools intersect -a mergedSquireRepRNA/SQuireRepeatsValidReadCounts_Aggregated_FWD_REV.bed -b /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed -wa -wb -f 0.75 > mergedSquireRepRNA/overlaps75SquireRep_vs_L1_L1Base.bed"
cmdLINEOverlapsWaWb = "bedtools intersect -a mergedSquireRepRNA/SQuireL1ValidReadCounts_Aggregated_FWD_REV.bed -b /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed -wa -wb -f 0.75 > mergedSquireRepRNA/overlaps75SquireLINE_vs_L1_L1BaseWaWb.bed"
cmdLINEOverlapsNoWaWb = "bedtools intersect -a mergedSquireRepRNA/SQuireL1ValidReadCounts_Aggregated_FWD_REV.bed -b /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed -f 0.75 > mergedSquireRepRNA/overlaps75SquireLINE_vs_L1_L1BaseNoWaWb.bed"

##those with and without sufficient read counts
cmdLINEAllOverlapsNoWaWb = "bedtools intersect -a mergedSquireRepRNA/SQuireL1ValidnNonValidReadCounts_Aggregated_FWD_REV.bed -b /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed -f 0.75 > mergedSquireRepRNA/overlaps75SquireLINE1All_vs_L1_L1BaseNoWaWb.bed"
# cmd = "bedtools intersect -a mergedSquireRepRNA/SQuireRepeatsValidReadCounts_Aggregated_FWD_REV_totCounts.bed -b /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed -wa -wb -f 0.75 > overlapsSquireRep_vs_L1_L1Base.bed"




system(command = paste0(cmd))

#oevrlaps with L1
system(command = paste0(cmdLINEOverlapsWaWb))
system(command = paste0(cmdLINEOverlapsNoWaWb))
system(command = paste0(cmdLINEAllOverlapsNoWaWb))

## repeatsInterest <- c("SINE", "LINE")
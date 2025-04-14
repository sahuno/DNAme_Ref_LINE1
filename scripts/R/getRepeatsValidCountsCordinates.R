library(data.table)
pathsRepeats <-  "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DE_repeats/validReadCounts_ptCoding_LociSpecL1_ActiveL1_allOtherRepeats.tsv"


repProCodingDT <- fread(pathsRepeats)

repDT <- repProCodingDT[!gene.type == "protein_coding",]


repDT[repDT == "" | repDT == " "] <- NA

colsInterest <- c("te.name", "gene.type", "chr",   "start", "end", "name", "number", "strand", "consensus",  "clade",  "class", "te.name", "te.id", "genomicElementID")

colsBed <- c("chr",   "start", "end", "name", "number", "strand", "genomicElementID", "gene.type","consensus",  "clade",  "class", "te.name", "te.id")

repDT[,..colsInterest]

repBedDT <- repDT[,..colsBed]
# repDT[!is.na(L1UID_name), ]

# repDT[L1UID_name == "", ]
fwrite(repBedDT, sep = "\t", col.names = FALSE, file = "repeatswithValidCounts.bed")

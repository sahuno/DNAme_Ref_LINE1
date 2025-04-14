library(tidyverse)

pathsLocusSpecData <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DE_repeats/LocusSpecific/data"

### load RNAseq data/ full length active L1s
filesRDATA <- list.files(pathsLocusSpecData, pattern = ".RData", recursive = TRUE, full.names = TRUE) 

filesRDATA_STDBi <- filesRDATA[grepl("SETDB1i", filesRDATA)]

# filesRDATA <- filesRDATA[1:2]

loadedRData <- lapply(filesRDATA_STDBi, function(x) {
    e <- new.env()
    load(x, envir = e)
    as.list(e)
  })


  ##rename to list 
baseNames <- gsub(".*condition_", "", basename(filesRDATA_STDBi))
names(loadedRData) <- gsub(".RData", "", baseNames)



resultsCPMRepeats_list <- lapply(loadedRData, function(x) pivot_longer(x$results2SaveLspecL1$cpmDF, !te.id, names_to = "RNAsamples", values_to = "cpm"))

resultsCPMRepeats_DF <- plyr::ldply(resultsCPMRepeats_list, data.frame)
resultsCPMRepeats_DF <- mutate(resultsCPMRepeats_DF, sampleID = gsub("R.", "sN.", RNAsamples))
unique(resultsCPMRepeats_DF$sampleID)



# names(loadedRData[["DMSO_SETDB1i"]]$results2SaveLspecL1)
# loadedRData[["DMSO_SETDB1i"]]$results2SaveLspecL1[["resDF"]]
# loadedRData[["DMSO_SETDB1i"]]$results2SaveLspecL1[["cpmDF"]]


# dim(loadedRData[["DMSO_SETDB1i"]]$results2SaveLspecL1[["resDF"]])

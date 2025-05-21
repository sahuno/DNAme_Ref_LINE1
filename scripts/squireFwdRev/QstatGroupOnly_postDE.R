#
# merge DNA methylation rates with RNAseq counts
# mkdir -p DNAme_RNA_rln
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/DE_RepeatsHeatmaps.R


library(data.table)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(UpSetR)
library(optparse)
library(pvclust)
library(ggnewscale)
library(viridisLite)
library(cola)
library(knitr)
library(tidyplots)
library(factoextra)
library(clusterProfiler)
library(DOSE)
library(AnnotationDbi)
# Directory to save plots (create if needed)
library(ggplot2)
library(clusterProfiler)


if(!require('NbClust')) {
  install.packages('NbClust')
  library('NbClust')
}


cividis_palette <- viridis(100, option = "cividis")
options("width"=200)
# set_fst_threads <- 4
# threads_fst(set_fst_threads)

# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/DE_RepeatsHeatmaps.R \
# --lfcCutoff 0.05 \
# --l1Separately FALSE \
# --paths_rdata "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/LocusSpecific" \
# --FulllengthAnno "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/overlaps75SquireLINE_vs_L1_L1Base_wa.bed"

option_list <- list(
    make_option(c("-l", "--lfcCutoff"), type="numeric", default=0.05, help="Cutoff for adjusted p-value [default %default]"),
    make_option(c("-s", "--l1Separately"), type="logical", default=FALSE, help="run locus specific line-1 separately [default %default]"),
    make_option(c("-d", "--paths_rdata"), type="character", default=NULL,help="Path to Rdata folder"),
    make_option(c("-f", "--FulllengthAnno"), type="character", default=NULL, help="Path to full-length LINE-1 annotation file")
    )

parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

opt$paths_rdata <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/LocusSpecific"
opt$FulllengthAnno <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/overlaps75SquireLINE_vs_L1_L1Base_wa.bed"

# After parsing, split comma-separated values into vectors
paths_rdata <- opt$paths_rdata


## Therapy list
combinedTherapy <- c("SETDB1i-CKi_SETDB1i" ,"SETDB1i-CKi_CKi","QSTAT-CKi_QSTAT" ,"QSTAT-CKi_CKi","SETDB1i-CKi_QSTAT-CKi")
monTherapy <- c("AZA_DMSO", "CKi_DMSO","QSTAT_DMSO", "QSTAT-CKi_DMSO", "SETDB1i_DMSO", "SETDB1i-CKi_DMSO")

#customize the therapy list
subSampleConds <-  c("CKi_DMSO", "QSTAT_DMSO", "QSTAT-CKi_DMSO")


pat <- paste(subSampleConds, collapse="|")

pat2 <- paste0("_condition_(",
               paste(subSampleConds, collapse="|"),
               ")\\.RData$")



if(opt$l1Separately){
  message("running locus specific LINE-1 separately")
} else {
  message("skipping locus specific LINE-1 separately")
}






###############################################################
###### Repeat analysis for all repeats
###############################################################
filesRDATAReps <- list.files(opt$paths_rdata, pattern = "Differential_locusSpecific_RepeatsExpressionResults", recursive = TRUE, full.names = TRUE) 

filesRDATAReps <- str_subset(filesRDATAReps, pat2)

loadedRDataReps <- lapply(filesRDATAReps, function(x) {
    e <- new.env()
    load(x, envir = e)
    as.list(e)
  })

names(loadedRDataReps[[1]][[1]])
# loadedRDataReps[[]]
##rename to list 
resultsDT <- lapply(loadedRDataReps, function(x){
  res <- x[["results2SaveLspecRepeats"]]$resDF
  res <- as.data.frame(res)
  res$RepID <- row.names(res)
  return(res)
})
names(resultsDT) <- gsub(".RData", "", gsub(".*condition_", "", filesRDATAReps))

# table(resultsDT$condition)

resultsDTSample <- rbindlist(resultsDT, idcol = "condition")
# table(resultsDTSample$condition)


# Check if 'fullLengthpath' is a column in resultsDTSample
if ("isfullLengthL1" %in% colnames(resultsDTSample)) {
  message("'isfullLengthL1' column exists in resultsDTSample.")
} else {
  message("'isfullLengthL1' column does not exist in resultsDTSample.")
  message("adding isfullLengthL1 column to resultsDTSample.")

if (file.exists(opt$FulllengthAnno) && !file.info(opt$FulllengthAnno)$isdir) {
  message(paste(opt$FulllengthAnno, "exists and is a file."))
colInterest <- "V11" # this will cause an error
fullDT <- fread(opt$FulllengthAnno, select = colInterest)
# if any entry in fullDT[[colInterest]] matched resultsDTSample$RepID create a new column `isfullLengthL1` in resultsDTSample with TRUE/FALSE
resultsDTSample[, isfullLengthL1 := RepID %in% fullDT[[colInterest]]]
} else {
  message(paste(opt$FulllengthAnno, "does not exist or is not a file."))
}
}

# fullLengthpath
# resultsDTSample[isfullLengthL1==TRUE,]

#Assign labels based on log2FoldChange and adjusted p-value
# resultsDTSample <- resultsDTSample %>% mutate(LFClabel = case_when(log2FoldChange > 1 & padj <= opt$lfcCutoff ~ "up",
#                                      log2FoldChange < -1 & padj <= opt$lfcCutoff ~ "down",
#                                      TRUE ~  "no_change")) 

resultsDTSample <- resultsDTSample %>% dplyr::mutate(LFClabel = case_when(
    padj < opt$lfcCutoff & log2FoldChange < -1 ~ "down",
    padj < opt$lfcCutoff & log2FoldChange > 1 ~ "up",
    padj > opt$lfcCutoff & log2FoldChange < -1 ~ "no_change",
    padj > opt$lfcCutoff & log2FoldChange > 1 ~ "no_change",
    padj < opt$lfcCutoff & abs(log2FoldChange) < 1 ~ "no_change",
    padj > opt$lfcCutoff & abs(log2FoldChange) < 1 ~ "no_change",
    TRUE ~ NA))





table(resultsDTSample$condition)

fwrite(resultsDTSample, file = file.path("resultsAllRepeatsRNA_DE_QSTATiOnlywithMetadata.txt"), sep = "\t")

# awk 'NR==1 || $10 == "TRUE"' resultsAllRepeatsRNA_DE_withMetadata.txt > resultsFromAllRepeatsRNA_DE_FullLengthL1_withMetadata.txt

################################################################################################
######## Repeats up and down counts ###########################################################
################################################################################################
resultsDTSample <- resultsDTSample %>% filter(LFClabel %in% c("up", "down")) %>% 
mutate(ClassFamily = str_split_i(RepID, "\\|", 4)) %>%
separate(ClassFamily, into = c("subfamily", "family", "class"), sep = ":") %>%
filter(class %in% c("SINE","LTR","LINE","DNA"))

# table(resultsDTSample$class)

head(resultsDTSample)
table(resultsDTSample$LFClabel)
resultsDTSampleUpDown <- resultsDTSample %>% filter(LFClabel %in% c("up", "down"))



# resultsDTSampleUpDown <- resultsDTSampleUpDown %>% 
# mutate(cordinate = str_split_i(RepID, "\\|", 1-3)) %>%
# separate(ClassFamily, into = c("chr", "family", "class"), sep = "|")
##separate_wider_delim()


resultsDTSampleUpDown <- resultsDTSampleUpDown %>%
  separate(RepID,
           into = c("chr","start","end","RepName"),
           sep  = "\\|",
           extra = "drop",
           fill  = "right", remove = FALSE) %>% 
           mutate(GenomicCoords = paste0(chr, ":", start, "-", end)) #%>% select(!`junk`)

library(readr)
#library(tidyverse)
#library(dplyr)
resultsDTSampleUpDownForBed <- resultsDTSampleUpDown %>% dplyr::select( c(RepID,  chr ,   start   ,   end  , RepName, isfullLengthL1, LFClabel,  subfamily,    family, class ,  LFClabel , GenomicCoords)) %>% distinct()

write_tsv(resultsDTSampleUpDownForBed %>% dplyr::select(chr, start, end, RepName), col_names = FALSE, file = "repeatsBed/resultsDTSampleUpDown.bed")
write_tsv(resultsDTSampleUpDownForBed %>% dplyr::select(GenomicCoords, RepName), col_names = FALSE, file = "repeatsBed/igvr_resultsDTSampleUpDown.txt")
# write_tsv(resultsDTSampleUpDown %>% select(GenomicCoords, RepName), col_names = FALSE, file = "igvr_resultsDTSampleUpDown.txt")

resultsDTSampleUpDownForBed %>% 
  group_by(class) %>% 
  group_walk(~ {
    write_tsv(
      .x %>% dplyr::select(GenomicCoords, RepName),
      file = paste0("repeatsBed/igvr_resultsDTSampleUpDown_", .y$class, ".txt"),
      col_names = FALSE
    )
  })


resultsDTSampleUpDownForBed %>% 
  group_by(class) %>% 
  group_walk(~ {
    write_tsv(
      .x %>% 
      arrange(chr, start, end) %>%
      select(chr, start, end, RepName),
      file = paste0("repeatsBed/Repeats_DE_UpDown_", .y$class, ".bed"),
      col_names = FALSE
    )
  })




# gr_short
# for each class‐specific bed…
bash_cmd <- '
for f in repeatsBed/Repeats_DE_UpDown_*.bed; do
  sort -k1,1 -k2,2n "$f" > "${f%.bed}.sorted.bed"
  mv "${f%.bed}.sorted.bed" "$f"
done
'

# system2 will invoke bash -c "<our script>"
system2(
  "bash",
  args = c("-c", shQuote(bash_cmd)),
  stdout = "",   # print stdout
  stderr = ""    # print stderr
)
head(resultsDTSampleUpDown)





#####################################################################
#####################################################################
### create promters for LINE1
# df_short <- head(resultsDTSampleUpDownForBed %>% filter(family == "L1")  %>% select(chr, start, end, RepName))

# if ("package:org.Mm.eg.db" %in% search()) {
#   detach("package:org.Mm.eg.db", unload = TRUE, character.only = TRUE)
# }
# if ("package:AnnotationDbi" %in% search()) {
#   detach("package:AnnotationDbi", unload = TRUE, character.only = TRUE)
# }
# if ("package:GenomicFeatures" %in% search()) {
#   detach("package:GenomicFeatures", unload = TRUE, character.only = TRUE)
# }

head(resultsDTSampleUpDownForBed)
LINE1Only <- resultsDTSampleUpDownForBed %>% filter(family == "L1")  %>% mutate(RepID = paste0(RepID,"|",isfullLengthL1)) %>% select(chr, start, end, RepID) 
LINE1Only_gr <- makeGRangesFromDataFrame(LINE1Only, keep.extra.columns = TRUE, ignore.strand = TRUE)
head(LINE1Only_gr)
head(LINE1Only)

resize_length <- c(100, 300, 500)

lapply(resize_length, function(x) {
  # 1) resize to the first x bp (clipped at the end if shorter)
  gr_resized <- resize(LINE1Only_gr,
                       width = x,
                       fix   = "start")
  
  names(gr_resized) <- gr_resized$RepID
  
  # 2) build a filename, e.g. "first_100bp.bed"
  outfile <- paste0("LINE1_promoters/first_", x, "bp.bed")

  # 3) export as BED
  export(gr_resized, con = outfile, format = "bed")
  
  # return the path invisibly if you like
  invisible(outfile)
})


#sort and order bed files
# for each class‐specific bed…
bash_cmdSortL1 <- '
for f in LINE1_promoters/first_*.bed; do
  sort -k1,1 -k2,2n "$f" > "${f%.bed}.sorted.bed"
  mv "${f%.bed}.sorted.bed" "$f"
done
'

# system2 will invoke bash -c "<our script>"
system2(
  "bash",
  args = c("-c", shQuote(bash_cmdSortL1)),
  stdout = "",   # print stdout
  stderr = ""    # print stderr
)
# head(resultsDTSampleUpDown)
# resultsDTSampleUpDown %>% 


write_tsv(resultsDTSampleUpDown, file = "resultsDTSampleUpDown_with_chromCord.txt", col_names = TRUE)

write_tsv(resultsDTSampleUpDownForBed, file = "uniqueRepeats_no_conditions_column_resultsDTSampleUpDown_with_chromCord.txt", col_names = TRUE)

#####################################################################
#####################################################################








resultsDTSample %>% filter(str_detect(condition, "SETDB1i")) %>% head() #(condition, LFClabel) %>% summarise(count = n(), .groups = "drop") %>% arrange(-count)
resultsDTSampleUpDown %>% filter(str_detect(condition, "SETDB1i")) %>% head() #(condition, LFClabel) %>% summarise(count = n(), .groups = "drop") %>% arrange(-count)

majorRepeatsClasses <- c("SINE","LTR","LINE","DNA")

# pltUPDown <- ggplot(resultsDTSampleUpDown, aes(x = condition, y = log2FoldChange)) +
#   geom_boxplot() +
#   geom_hline(yintercept = 0, color = "black") +
#   labs(title = "Up and Down Sig. Repeats per Condition",
#        x = "Condition",
#        y = "log2FoldChange",
#        color = "Label") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~ class, scales = "free_y")    # Facet by repeat class
# ggsave(pltUPDown, file = "RepeatsUp_down_LFC_per_Condition.png", width = 10, height = 6)


pltUPDown <- resultsDTSampleUpDown %>% tidyplot(x = condition, y = log2FoldChange) |> 
    add_boxplot(alpha = 0.4) |> 
  adjust_x_axis(rotate_labels = 45)|>
  split_plot(by = class, nrow = 3) |>
      save_plot(filename = "RepeatsAllClassesUp_down_LFC_QSTATi_GroupSamePage.pdf", view_plot = FALSE, multiple_files = FALSE)

#   save_plot(filename = "RepeatsAllClassesUp_down_LFC_QSTATi_Group.pdf")

resultsDTSampleUpDown %>%
  filter(class %in% c("SINE","LTR","LINE","DNA")) %>% filter(is.na(log2FoldChange))

DFmajorClasses <- resultsDTSampleUpDown %>%
  filter(class %in% c("SINE","LTR","LINE","DNA")) %>%
  count(class, condition)

counts_df <- resultsDTSampleUpDown %>%
  filter(class %in% c("SINE","LTR","LINE","DNA")) %>%
  group_by(class, condition) %>%
  summarise(
    n       = n(),
    y_top   = max(log2FoldChange),
    .groups = "drop"
  )



  pltUPDownMajorClasses <- resultsDTSampleUpDown %>% filter(class %in% c("SINE","LTR","LINE","DNA")) %>% tidyplot(x = condition, y = log2FoldChange) |> 
    add_boxplot(alpha = 0.4) |> 
    add_test_asterisks(hide.ns = FALSE, method = "wilcox_test", hide_info = TRUE) |>
  adjust_x_axis(rotate_labels = 45)|>
#   split_plot(by = class, ncol = 2, nrow = 2) |>
#   save_plot(filename = "RepeatsMajorClassesUp_down_LFC_QSTATi_Group.pdf", view_plot = FALSE) |>
  split_plot(by = class, nrow = 3) |>
    save_plot(filename = "RepeatsMajorClassesUp_down_LFC_QSTATi_GroupSamePage.pdf", view_plot = FALSE, multiple_files = FALSE)


  pltVolMajorClassesQCKI <- resultsDTSampleUpDown %>% 
  filter(class %in% c("SINE","LTR","LINE","DNA"), condition %in% c("QSTAT-CKi_DMSO")) %>% mutate(`-log10(padj)` = -log10(padj)) %>%
tidyplot(x = log2FoldChange, y = `-log10(padj)`, color = LFClabel) |> 
    add_data_points(alpha = 0.4) |> 
    # add_test_asterisks(hide.ns = FALSE, method = "wilcox_test", hide_info = TRUE) |>
  adjust_x_axis(rotate_labels = 45)|>
      add_reference_lines(x = c(-1, 1), y = 1.3,linetype = "dashed",linewidth = 0.25) |>
   add_title("QSTAT-CKi_DMSO") |>
  split_plot(by = c(class)) |>
#   add_title("QSTAT-CKi_DMSO") |>
  save_plot(filename = "RepeatsMajorClassesVolcano_LFC_QSTAT_CKIi_Group.pdf")


  pltVolMajorClassesCKI <- resultsDTSampleUpDown %>% 
  filter(class %in% c("SINE","LTR","LINE","DNA"), condition %in% c("CKi_DMSO")) %>% mutate(`-log10(padj)` = -log10(padj)) %>%
tidyplot(x = log2FoldChange, y = `-log10(padj)`, color = LFClabel) |> 
    add_data_points(alpha = 0.4) |> 
    # add_test_asterisks(hide.ns = FALSE, method = "wilcox_test", hide_info = TRUE) |>
  adjust_x_axis(rotate_labels = 45)|>
      add_reference_lines(x = c(-1, 1), y = 1.3,linetype = "dashed",linewidth = 0.25) |>
  add_title("CKi_DMSO") |>
#   rename_y_axis_labels("padj (-log10)") |>
  split_plot(by = c(class)) |>
  save_plot(filename = "RepeatsMajorClassesVolcano_LFC_CKIi_Group.pdf")


  pltVolMajorClassesQSTATi <- resultsDTSampleUpDown %>% 
  filter(class %in% c("SINE","LTR","LINE","DNA"), condition %in% c("QSTAT_DMSO")) %>% mutate(`-log10(padj)` = -log10(padj)) %>%
tidyplot(x = log2FoldChange, y = `-log10(padj)`, color = LFClabel) |> 
    add_data_points(alpha = 0.4) |> 
    add_reference_lines(x = c(-1, 1), y = 1.3,linetype = "dashed",linewidth = 0.25) |>
    # add_test_asterisks(hide.ns = FALSE, method = "wilcox_test", hide_info = TRUE) |>
  adjust_x_axis(rotate_labels = 45)|>
  add_title("QSTAT_DMSO") |>
  split_plot(by = c(class)) |>
  save_plot(filename = "RepeatsMajorClassesVolcano_LFC_QSTATi_Group.pdf")


head(resultsDTSampleUpDown)
  pltVolSINESQSTATi <- resultsDTSampleUpDown %>% 
  filter(class %in% c("SINE","LTR","LINE","DNA"), condition %in% c("QSTAT-CKi_DMSO"), class == "SINE") %>% mutate(`-log10(padj)` = -log10(padj)) %>%
tidyplot(x = log2FoldChange, y = `-log10(padj)`, color = LFClabel) |> 
    add_data_points(alpha = 0.4) |> 
    # add_test_asterisks(hide.ns = FALSE, method = "wilcox_test", hide_info = TRUE) |>
  adjust_x_axis(rotate_labels = 45)|>
      add_reference_lines(x = c(-1, 1), y = 1.3,linetype = "dashed",linewidth = 0.25) |>
   add_title("QSTAT-CKi_DMSO") |>
  split_plot(by = c(family)) |>
#   add_title("QSTAT-CKi_DMSO") |>
  save_plot(filename = "RepeatsSINESVolcano_LFC_QSTAT_CKIi_Group.pdf", view_plot = FALSE)



  pltVolLTRQSTATi <- resultsDTSampleUpDown %>% 
  filter(class %in% c("SINE","LTR","LINE","DNA"), condition %in% c("QSTAT-CKi_DMSO"), class == "LTR") %>% mutate(`-log10(padj)` = -log10(padj)) %>%
tidyplot(x = log2FoldChange, y = `-log10(padj)`, color = LFClabel) |> 
    add_data_points(alpha = 0.4) |> 
    # add_test_asterisks(hide.ns = FALSE, method = "wilcox_test", hide_info = TRUE) |>
  adjust_x_axis(rotate_labels = 45)|>
      add_reference_lines(x = c(-1, 1), y = 1.3,linetype = "dashed",linewidth = 0.25) |>
   add_title("QSTAT-CKi_DMSO") |>
  split_plot(by = c(family)) |>
#   add_title("QSTAT-CKi_DMSO") |>
  save_plot(filename = "RepeatsLTRVolcano_LFC_QSTAT_CKIi_Group.pdf")



# majorRepeatsClasses

df_filteredCondClass <- resultsDTSample %>% select(condition, RepID, class, family, log2FoldChange, padj, LFClabel) %>%
  arrange(-log2FoldChange) %>%
  group_by(condition, class, LFClabel) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(logCounts = log10(count + 1)) %>%
  mutate(logCounts = ifelse(LFClabel == "down", -logCounts, logCounts)) %>%
  mutate(count = ifelse(LFClabel == "down", -count, count))

df_filteredCondMajorClass <- resultsDTSample %>% filter(class %in% majorRepeatsClasses) %>% select(condition, RepID, class, family, log2FoldChange, padj, LFClabel) %>%
  arrange(-log2FoldChange) %>%
  group_by(condition, class, LFClabel) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(logCounts = log10(count + 1)) %>%
  mutate(logCounts = ifelse(LFClabel == "down", -logCounts, logCounts)) %>%
  mutate(count = ifelse(LFClabel == "down", -count, count))


pltClass <- ggplot(df_filteredCondClass, aes(x = condition, y = logCounts, fill = LFClabel)) +
  geom_bar(stat = "identity", position = "identity") +
  facet_wrap(~ class, scales = "free_y") +   # Facet by repeat class
  scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "Up and Down Sig. Repeats per Condition by Class",
       x = "Condition",
       y = "logCounts",
       fill = "Label") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(pltClass, file = "RepeatsUp_down_CountsPerGroup_per_Condition_Class.png", width = 10, height = 6)

pltClassUnlog <- ggplot(df_filteredCondClass, aes(x = condition, y = count, fill = LFClabel)) +
  geom_bar(stat = "identity", position = "identity") +
  facet_wrap(~ class, scales = "free_y") +   # Facet by repeat class
  scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "Up and Down Sig. Repeats per Condition by Class",
       x = "Condition",
       y = "Counts",
       fill = "Label") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(pltClassUnlog, file = "RepeatsUp_down_LinearScaleCountsPerGroup_per_Condition_Class.pdf", width = 10, height = 6)


df_filteredCondMajorClassCountsRaw <- resultsDTSample %>% filter(class %in% majorRepeatsClasses) %>% select(condition, RepID, class, family, log2FoldChange, padj, LFClabel) %>%
  arrange(-log2FoldChange) %>%
  group_by(condition, class, LFClabel) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(count = ifelse(LFClabel == "down", -count, count))


pltMajorClass <- ggplot(df_filteredCondMajorClass, aes(x = condition, y = logCounts, fill = LFClabel)) +
  geom_bar(stat = "identity", position = "identity") +
  facet_wrap(~ class, scales = "free_y") +   # Facet by repeat class
  scale_fill_manual(values = c("up" = "#D55E00", "down" = "#0072B2")) +
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "Up and Down Sig. Repeats per Condition by Class",
       x = "Condition",
       y = "logCounts",
       fill = "Label") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(pltMajorClass, file = "RepeatsUp_down_CountsPerGroup_per_Condition_MajorClasses.png", width = 10, height = 6)


pltMajorClassUnlog <- ggplot(df_filteredCondMajorClass, aes(x = condition, y = count, fill = LFClabel)) +
  geom_bar(stat = "identity", position = "identity") +
  facet_wrap(~ class, scales = "free_y") +   # Facet by repeat class
  scale_fill_manual(values = c("up" = "#D55E00", "down" = "#0072B2")) +
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "Up and Down Sig. Repeats per Condition by Class",
       x = "Condition",
       y = "Counts",
       fill = "Label") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(pltMajorClassUnlog, file = "RepeatsUp_down_CountsLinearPerGroup_per_Condition_MajorClasses.pdf", width = 10, height = 6)




df_filteredCondMajorClassFamily <- resultsDTSample %>% filter(class %in% majorRepeatsClasses) %>% select(condition, RepID, class, family, log2FoldChange, padj, LFClabel) %>%
  arrange(-log2FoldChange) %>%
  group_by(condition, class, family, LFClabel) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(logCounts = log10(count + 1)) %>%
  mutate(logCounts = ifelse(LFClabel == "down", -logCounts, logCounts)) %>%
  mutate(count = ifelse(LFClabel == "down", -count, count))



########################################
#### Here is the most important figure! 
########################################

pltMajorClassFamily <- ggplot(df_filteredCondMajorClassFamily, aes(x = reorder(family, -count), y = count, fill = condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ class, scales = "free", ncol = 2) +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_brewer(palette = "Set2") + # or scale_fill_manual(...)
  labs(
    title = "Counts per Family × Condition, by Class",
    x     = "Family",
    y     = "counts (neg if down)",
    fill  = "Condition"
  ) +
  # scale_y_continuous(breaks = y_breaks, limits = c(y_min, y_max)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text   = element_text(face = "bold")
  )

ggsave(pltMajorClassFamily, file = "RepeatsUp_down_CountsPerGroup_per_Condition_MajorClassesFamily.pdf", width = 10, height = 6)



y_rng   <- range(df_filteredCondMajorClassFamily$logCounts, na.rm=TRUE)
# round to the nearest 0.5
y_min   <- floor(y_rng[1] * 2) / 2
y_max   <- ceiling(y_rng[2] * 2) / 2
y_breaks <- seq(y_min, y_max, by = 0.5)

pltMajorClassFamilyLog <- ggplot(df_filteredCondMajorClassFamily, aes(x = family, y = logCounts, fill = condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ class, scales = "free", ncol = 2) +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_brewer(palette = "Set2") + # or scale_fill_manual(...)
  labs(
    title = "Counts per Family × Condition, by Class",
    x     = "Family",
    y     = "logCounts (neg if down)",
    fill  = "Condition"
  ) +
  # scale_y_continuous(breaks = y_breaks, limits = c(y_min, y_max)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text   = element_text(face = "bold")
  )

ggsave(pltMajorClassFamilyLog, file = "RepeatsUp_down_logCountsPerGroup_per_Condition_MajorClassesFamily.pdf", width = 10, height = 6)



resultsDTSample %>% filter(isfullLengthL1 == TRUE) %>% group_by(condition, LFClabel) %>% summarise(count = n(), .groups = "drop") %>% arrange(-count)
# resultsDTSample$log2FoldChange[is.na(resultsDTSample$log2FoldChange)] <- 0

# Use tidyr to reshape data: pivot so rows are LINE1 elements and columns are samples
heatmap_data <- resultsDTSample %>%
  select(condition, RepID, log2FoldChange) %>%
  pivot_wider(names_from = condition, values_from = log2FoldChange, values_fill = 0)

table(resultsDTSample$class)
# Setting the row names to repeats (RepID) and remove the RepID column from the data frame
heatmap_matrix <- as.data.frame(heatmap_data) %>% 
  column_to_rownames(var = "RepID") %>%
  as.matrix()

heatmap_matrixAdjusted <- cola::adjust_matrix(heatmap_matrix, sd_quantile = 0.05, max_na = 0.25, verbose = TRUE)

heatmap_matrixAdjusted_t = t(heatmap_matrixAdjusted)

# res <- NbClust(heatmap_matrixAdjusted, distance = "euclidean", min.nc=2, max.nc=8, 
#             method = "complete", index = "ch")



subSampleCondsDF <-  data.frame(conditions = subSampleConds)
DrugsAnnot <- subSampleCondsDF %>% dplyr::filter(conditions %in% subSampleConds) %>% 
mutate(Combinations = case_when(str_detect(conditions, "-") ~ "combo",TRUE ~ "mono"), 
BroadTargets = case_when(str_detect(conditions, "QSTAT_|SETDB1i_") ~ "Chromatin",
                          str_detect(conditions, "-CKi_") ~ "Chromatin+MEK",
                          str_detect(conditions, "AZA_") ~ "DNAme",
                          str_detect(conditions, "CKi_DMSO") ~ "MEK", TRUE ~ NA),
Action = case_when(str_detect(conditions, "QSTAT_") ~ "HDACi", 
                    str_detect(conditions, "SETDB1i_") ~ "H3K9mei",
                    str_detect(conditions, "QSTAT-CKi_DMSO") ~ "HDACi+MEKi",
                    str_detect(conditions, "SETDB1i-CKi_DMSO") ~ "H3K9mei+MEKi",
                    str_detect(conditions, "AZA_") ~ "DNMTi",
                    str_detect(conditions, "CKi_DMSO") ~ "MEKi", TRUE ~ NA))

DrugsAnnot <- DrugsAnnot %>% mutate(MEKi = case_when(str_detect(Action, "MEKi") ~ "Yes", TRUE ~ "No"))
# )
rownames(DrugsAnnot) <- DrugsAnnot$conditions

set.seed(540)
Koptimal=2
k <- pheatmap::pheatmap(heatmap_matrixAdjusted, kmeans_k = Koptimal, annotation = DrugsAnnot[,c("Combinations", "Action", "BroadTargets", "MEKi")], 
# colorRampPalette(c("#5D3A9B", "white", "#E66100"))(100),
main = paste0("abs(LFC) > 1 & padj < 0.05; repeats [NA=0]"), 
scale = "row", 
annotation_names_row=FALSE, 
show_rownames = F, 
filename = paste0("Kmeans_clustering_sig_repeats_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.png"))

# names(k$kmeans)
clusterDF <- as.data.frame(factor(k$kmeans$cluster))
colnames(clusterDF) <- "Cluster"
OrderByCluster <- heatmap_matrixAdjusted[order(clusterDF$Cluster), ]

set.seed(123)  # seed for reproducibility
clusterNum <- NbClust(heatmap_matrixAdjusted, distance = "euclidean", min.nc = 2, max.nc = 12, method = "kmeans", index = "silhouette")
message("optimal number of clusters: ", clusterNum$Best.nc[1], "; with value index: ", clusterNum$Best.nc[2])

head(clusterNum$Best.partition)
head(clusterNum$Best.nc)
head(clusterNum$All.index)
table(clusterNum$Best.nc[1])

# sort and order the clusters
orderedCluster <- sort(clusterNum$Best.partition)
heatmap_matrixAdjusted_sorted_NBClust <- heatmap_matrixAdjusted[match(names(orderedCluster), rownames(heatmap_matrixAdjusted)), ]


library(cluster)
set.seed(762)
pltclust <- fviz_nbclust(heatmap_matrixAdjusted_sorted_NBClust, nboot=1000, kmeans, method = "silhouette", nstart = 25)
ggsave("silhouette_optimal_clusters_repeatsV2.png", plot = pltclust, width = 6, height = 4, dpi = 300)

set.seed(542)
pltclustWSS <- fviz_nbclust(heatmap_matrixAdjusted_sorted_NBClust, nboot=1000, kmeans, method = "wss", nstart = 25)
ggsave("wss_optimal_clusters_repeats.png", plot = pltclustWSS, width = 6, height = 4, dpi = 300)

set.seed(321)
k <- clusterNum$Best.nc[1] 
km_res <- kmeans(heatmap_matrixAdjusted_sorted_NBClust, centers = k, nstart = 25)

# library(plotly)
set.seed(981)
plt_clusters <- fviz_cluster(km_res, data = heatmap_matrixAdjusted_sorted_NBClust, geom = "point",
             ellipse.type = "convex", palette = "jco",
             ggtheme = theme_minimal())
ggsave("kmeans_optimal_clusters_repeats.png", plot = plt_clusters, width = 6, height = 4, dpi = 300)

set.seed(981)
km_res2 <- kmeans(heatmap_matrixAdjusted_sorted_NBClust, centers = 5, nstart = 25)
plt_clustersk5 <- fviz_cluster(km_res2, data = heatmap_matrixAdjusted_sorted_NBClust, geom = "point",
             ellipse.type = "convex", palette = "jco",
             ggtheme = theme_minimal())
ggsave("kmeans5_clusters_repeats.png", plot = plt_clustersk5, width = 6, height = 4, dpi = 300)



set.seed(983)
km_resk4 <- kmeans(heatmap_matrixAdjusted_sorted_NBClust, centers = 4, nstart = 25)
plt_clustersk4 <- fviz_cluster(km_resk4, data = heatmap_matrixAdjusted_sorted_NBClust, geom = "point",
             ellipse.type = "convex", palette = "jco",
             ggtheme = theme_minimal())
ggsave("kmeans4_clusters_repeats.png", plot = plt_clustersk4, width = 6, height = 4, dpi = 300)




################################################################################################################
##kmeans 3 visualization #######################################################################################
################################################################################################################

set.seed(983)
km_resk3 <- kmeans(heatmap_matrixAdjusted_sorted_NBClust, centers = 3, nstart = 25)

plt_clustersk3 <- fviz_cluster(km_resk3, data = heatmap_matrixAdjusted_sorted_NBClust, geom = "point",
             ellipse.type = "convex", palette = "jco",
             ggtheme = theme_minimal())
ggsave("kmeans3_clusters_repeats.png", plot = plt_clustersk3, width = 6, height = 4, dpi = 300)

clusterDF_K3 <- as.data.frame(factor(km_resk3$cluster))
colnames(clusterDF_K3) <- "Cluster"

pheatmap::pheatmap(heatmap_matrixAdjusted_sorted_NBClust, 
annotation_row = clusterDF_K3,
scale = "row", show_rownames = FALSE,
    cluster_rows = FALSE,
    annotation = DrugsAnnot[,c("Combinations", "Action", "BroadTargets", "MEKi")], 
    main = paste0("abs(LFC) > 1 & padj < 0.05; repeats [NA=0]"),
    filename = paste0("kmeans3_after_pca_visualization_sig_Repeats_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.pdf"))

save(heatmap_matrixAdjusted_sorted_NBClust, orderedCluster, DrugsAnnot, clusterNum, km_resk4,file = "data.RData")


head(clusterDF_K3_df)

clusterDF_K3_df <- clusterDF_K3 %>%  rownames_to_column(var="RepID")  %>% 
mutate(ClassFamily = str_split_i(RepID, "\\|", 4)) %>%
separate(ClassFamily, into = c("subfamily", "family", "class"), sep = ":")

clusterDF_K3_df <- clusterDF_K3_df %>% 
  separate(RepID,
           into = c("chr","start","end","RepName"),
           sep  = "\\|",
           extra = "drop",
           fill  = "right", remove = FALSE) %>% 
           mutate(GenomicCoords = paste0(chr, ":", start, "-", end)) #%>% head()


head(clusterDF_K3_df)
# Calculate proportions
clusterDF_K3_class_prop <- clusterDF_K3_df %>% 
  group_by(Cluster, class) %>% 
  summarise(n = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(prop = n / sum(n))

clusterDF_K3_MajorClass_prop <- clusterDF_K3_df %>% filter(class %in% c("SINE","LTR","LINE","DNA"))  %>% 
  group_by(Cluster, class) %>% 
  summarise(n = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(prop = n / sum(n))



table(clusterDF_K3_class_prop$Cluster)

# Plot
plt_class_propk3 <- ggplot(clusterDF_K3_class_prop, aes(x = factor(Cluster), y = prop, fill = class)) +
  geom_col(position = "fill") +
  labs(
    title = "Proportion of Repeat Classes in Each Cluster",
    x = "Cluster",
    y = "Proportion",
    fill = "Class"
  ) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()

ggsave(plt_class_propk3, file = "proportion_repeat_class_per_clusterk3.pdf", width = 7, height = 5)


message("Proportion of Major Classes in Each Cluster")
clusterDF_K3_MajorClass_prop
write_tsv(clusterDF_K3_MajorClass_prop, file = "proportion_majorRepeat_class_per_clusterk3.tsv")

plt_majoClass_propk3 <- ggplot(clusterDF_K3_MajorClass_prop, aes(x = factor(Cluster), y = prop, fill = class)) +
  geom_col(position = "fill") +
  labs(
    title = "Proportion of Repeat Classes in Each Cluster",
    x = "Cluster",
    y = "Proportion",
    fill = "Class"
  ) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()

ggsave(plt_majoClass_propk3, file = "proportion_MajorRepeat_class_per_clusterk3.pdf", width = 7, height = 5)
####################################################################################################################################
######################## K means 3##################################################################################################
#####################################################################################################################################






################################################################################################################
###### Generic functions for different kmeans ##################################################################
################################################################################################################

## Function to run k-means clustering, visualize results, and summarize class proportions for any k
## Function to run k-means clustering, visualize results, and summarize class proportions for any k
cluster_heatmap_kmeans <- function(
  data_matrix,
  k,
  annotation_df,
  output_dir = ".",
  prefix = NULL,
  nstart = 25,
  seed = 983
) {
  # Load required packages
  require(ggplot2)
  require(factoextra)
  require(pheatmap)
  require(dplyr)
  require(tidyr)
  require(readr)
  require(stringr)
  require(tibble)
  require(scales)

  # Prepare output directory and filename prefix
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(prefix)) {
    prefix <- paste0("kmeans", k)
  }

  # Ensure annotation_df has proper rownames for pheatmap
  annotation_df <- as.data.frame(annotation_df)
  if ("conditions" %in% colnames(annotation_df)) {
    rownames(annotation_df) <- annotation_df$conditions
    annotation_df$conditions <- NULL
  }
  # Check matching
  if (!all(colnames(data_matrix) %in% rownames(annotation_df))) {
    stop("Row names of annotation_df must match column names of data_matrix.")
  }
  # Reorder annotation_df to match data_matrix
  annotation_df <- annotation_df[colnames(data_matrix), , drop = FALSE]

  # Run k-means
  set.seed(seed)
  km_res <- kmeans(data_matrix, centers = k, nstart = nstart)

  # Build clusterDF with RepID preserved
  clusterDF <- tibble(
    RepID   = names(km_res$cluster),
    Cluster = factor(km_res$cluster)
  ) %>%
    separate(
      col    = RepID,
      into   = c("chr", "start", "end", "RepName"),
      sep    = "\\|",
      remove = FALSE,
      extra  = "drop",
      fill   = "right"
    ) %>%
    mutate(
      GenomicCoords = paste0(chr, ":", start, "-", end),
      ClassFamily   = str_split(RepName, ":", simplify = TRUE)[,1]
    ) %>%
    separate(RepName, into = c("subfamily", "family", "class"), sep = ":")

  # 1) Scatter plot of clusters
  scatter_plot <- fviz_cluster(
    km_res,
    data         = data_matrix,
    geom         = "point",
    ellipse.type = "convex",
    palette      = "jco",
    ggtheme      = theme_minimal()
  )
  ggsave(
    filename = file.path(output_dir, paste0(prefix, "_clusters.png")),
    plot     = scatter_plot,
    width    = 6,
    height   = 4,
    dpi      = 300
  )

  # 2) Heatmap with annotations
  annot_row <- data.frame(Cluster = clusterDF$Cluster)
  rownames(annot_row) <- clusterDF$RepID
  pheatmap(
    mat            = data_matrix,
    annotation_row = annot_row,
    annotation_col = annotation_df,
    scale          = "row",
    show_rownames  = FALSE,
    cluster_rows   = FALSE,
    main           = paste0("k = ", k, " clusters"),
    filename       = file.path(output_dir, paste0(prefix, "_heatmap.pdf"))
  )

  # 3) Calculate class proportions
  prop_df <- clusterDF %>%
    group_by(Cluster, class) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Cluster) %>%
    mutate(prop = count / sum(count))

  # 4) Subset to major repeat classes
  major_classes <- c("SINE", "LTR", "LINE", "DNA")
  prop_major <- prop_df %>%
    filter(class %in% major_classes)

  # 5) Plot class proportions
  class_plot <- ggplot(prop_df, aes(x = Cluster, y = prop, fill = class)) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = percent) +
    labs(
      title = paste0("Repeat class proportions (k=", k, ")"),
      x     = "Cluster",
      y     = "Proportion",
      fill  = "Class"
    ) +
    theme_minimal()
  ggsave(
    filename = file.path(output_dir, paste0(prefix, "_class_prop.pdf")),
    plot     = class_plot,
    width    = 7,
    height   = 5
  )

  # 6) Plot major class proportions
  major_plot <- ggplot(prop_major, aes(x = Cluster, y = prop, fill = class)) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = percent) +
    labs(
      title = paste0("Major repeat class proportions (k=", k, ")"),
      x     = "Cluster",
      y     = "Proportion",
      fill  = "Class"
    ) +
    theme_minimal()
  ggsave(
    filename = file.path(output_dir, paste0(prefix, "_major_class_prop.pdf")),
    plot     = major_plot,
    width    = 7,
    height   = 5
  )

  # 7) Save results and return
  write_tsv(prop_df,    file = file.path(output_dir, paste0(prefix, "_prop.tsv")))
  write_tsv(prop_major, file = file.path(output_dir, paste0(prefix, "_major_prop.tsv")))
  save(
    data_matrix, clusterDF, annotation_df, k, km_res,
    file = file.path(output_dir, paste0(prefix, "_results.RData"))
  )

  invisible(list(
    kmeans    = km_res,
    clusterDF = clusterDF,
    plots     = list(
      scatter     = scatter_plot,
      class_prop  = class_plot,
      major_prop  = major_plot
    )
  ))
}

# Example usage
# cluster_heatmap_kmeans(
#   data_matrix    = heatmap_matrixAdjusted_sorted_NBClust,
#   k              = 3,
#   annotation_df  = DrugsAnnot[, c("conditions","Combinations","BroadTargets","Action","MEKi")],
#   output_dir     = "results"
# )


# run cluster for k=2:
cluster_heatmap_kmeans(
  data_matrix    = heatmap_matrixAdjusted_sorted_NBClust,
  k              = 2,
  annotation_df  = DrugsAnnot[, c("Combinations", "Action", "BroadTargets", "MEKi")],
  output_dir     = "results"
)

k3_clustering_res <- cluster_heatmap_kmeans(
  data_matrix    = heatmap_matrixAdjusted_sorted_NBClust,
  k              = 3,
  annotation_df  = DrugsAnnot[, c("Combinations", "Action", "BroadTargets", "MEKi")],
  output_dir     = "results_k3"
)

################################################################################################################
################################End of Generic function ##############################################
################################################################################################################



# colnames(matrixMonoCombo_DMSO_ref_cleaned)
################################################################################################################
######################### heatmaps with optimal clustering #####################################################
################################################################################################################
pheatmap::pheatmap(heatmap_matrixAdjusted_sorted_NBClust, 
annotation_row = clusterDF,
scale = "row", show_rownames = FALSE,
    cluster_rows = FALSE,
    annotation = DrugsAnnot[,c("Combinations", "Action", "BroadTargets", "MEKi")], 
    main = paste0("abs(LFC) > 1 & padj < 0.05; repeats [NA=0]"),
    filename = paste0("kmeans_withOptimumClustering_sig_Repeats_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.pdf"))


# Plot proportions of class in each orderedCluster
# library(ggplot2)
# library(dplyr)
################################################################
################################################################
################################################################

df_orderedCluster <- data.frame(orderedCluster)
df_orderedCluster <- df_orderedCluster %>% rownames_to_column(var = "RepID")
head(df_orderedCluster)

df_orderedCluster <- df_orderedCluster %>%  
mutate(ClassFamily = str_split_i(RepID, "\\|", 4)) %>%
separate(ClassFamily, into = c("subfamily", "family", "class"), sep = ":")




# Calculate proportions
df_class_prop <- df_orderedCluster %>% 
  group_by(orderedCluster, class) %>% 
  summarise(n = n(), .groups = "drop") %>%
  group_by(orderedCluster) %>%
  mutate(prop = n / sum(n))

table(df_orderedCluster$orderedCluster)

# Plot
plt_class_prop <- ggplot(df_class_prop, aes(x = factor(orderedCluster), y = prop, fill = class)) +
  geom_col(position = "fill") +
  labs(
    title = "Proportion of Repeat Classes in Each Cluster",
    x = "Cluster",
    y = "Proportion",
    fill = "Class"
  ) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()

ggsave(plt_class_prop, file = "proportion_repeat_class_per_optimal_cluster.pdf", width = 7, height = 5)



df_Major_class_prop <- df_orderedCluster %>% 
 filter(class %in% c("SINE","LTR","LINE","DNA"))%>%
  group_by(orderedCluster, class) %>% 
  summarise(n = n(), .groups = "drop") %>%
  group_by(orderedCluster) %>%
  mutate(prop = n / sum(n))


plt_MajorClass_prop <- ggplot(df_Major_class_prop, aes(x = factor(orderedCluster), y = prop, fill = class)) +
  geom_col(position = "fill") +
  labs(
    title = "Proportion of Repeat Classes in Each Cluster",
    x = "Cluster",
    y = "Proportion",
    fill = "Class"
  ) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()

ggsave(plt_MajorClass_prop, file = "proportion_MajorRepeat_class_per_optimal_cluster.pdf", width = 7, height = 5)




##################################################################################################################
###################################### End of optimal clustering #################################################
##################################################################################################################

##################################################################################################################
############### Annotation #######################################################################################
##################################################################################################################
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rtracklayer)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)  # GRanges of all known genes





# names(mergedRNA_Annot) <- names(splitRNA)
# head(mergedRNA_Annot[[1]])
# dim(mergedRNA_Annot[[1]])
# table(mergedRNA_Annot[[1]]$condition)
# dim(clusterDF_K3_df)




#############################################################################################
##### Gene set enrichment analysis ##########################################################
#############################################################################################


# ─── 0. Install & load packages ────────────────────────────────────────────────
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# pkgs <- c("clusterProfiler", "org.Mm.eg.db", "DOSE", "dplyr", "AnnotationDbi")
# for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p)



##############################################################################
##########################.                         ##########################
##############################################################################
refactorCluster <- function(clusterDat){
df_orderedCluster <- data.frame(clusterDat)
names(df_orderedCluster) <- c("Cluster")
df_orderedCluster <- df_orderedCluster %>% rownames_to_column(var = "RepID")
message("head of cluster data")
head(df_orderedCluster)
df_Cluster <- df_orderedCluster %>%  
mutate(ClassFamily = str_split_i(RepID, "\\|", 4))  %>%
    separate(ClassFamily, into = c("subfamily", "family", "class"), sep = ":") %>% 
    separate(RepID,
           into = c("chr","start","end","RepName"),
           sep  = "\\|",
           extra = "drop",
           fill  = "right", remove = FALSE) %>% 
           mutate(GenomicCoords = paste0(chr, ":", start, "-", end))
return(df_Cluster)
}


#test code
# clusterDat <- km_res$cluster
# ttt <- refactorCluster(km_res$cluster)
# head(ttt) #clusterDat

annotate_repeat_clusters <- function(
  cluster_df,
  txdb,
  orgdb,
  annoDb       = "org.Mm.eg.db",
  class_filter = c("SINE","LTR","LINE","DNA"),
  tssRegion    = c(-3000, 3000),
  out_dir      = "repeat_cluster_results",
  prefix       = "cluster"
){
  # ─── Dependencies ───────────────────────────────────────────────────────────
  require(GenomicRanges)
  require(ChIPseeker)
  require(AnnotationDbi)
  require(dplyr)
  require(ggplot2)
  
  # ─── Prepare output folder ──────────────────────────────────────────────────
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
  
  # ─── 1) Clean types ─────────────────────────────────────────────────────────
  df <- cluster_df %>%
    mutate(
      chr   = as.character(chr),
      start = as.integer(start),
      end   = as.integer(end)
    )
  
  # ─── 2) Build GRanges ───────────────────────────────────────────────────────
  gr <- makeGRangesFromDataFrame(df,
                                 seqnames.field     = "chr",
                                 start.field        = "start",
                                 end.field          = "end",
                                 keep.extra.columns = TRUE)
  genes_gr <- GenomicFeatures::genes(txdb)
  
  # ─── 3) Intragenic vs Intergenic ────────────────────────────────────────────
  df$GenicStatus <- "intergenic"
  hits <- findOverlaps(gr, genes_gr)
  df$GenicStatus[queryHits(hits)] <- "intragenic"
  
  # ─── 4) Nearest gene & distance ─────────────────────────────────────────────
  nearest <- distanceToNearest(gr, genes_gr)
  matched <- queryHits(nearest)
  gene_ids <- genes_gr$gene_id[subjectHits(nearest)]
  syms     <- mapIds(orgdb,
                     keys     = gene_ids,
                     column   = "SYMBOL",
                     keytype  = "ENTREZID",
                     multiVals= "first")
  df$NearestGene    <- NA_character_
  df$DistanceToGene <- NA_integer_
  df$NearestGene[matched]    <- syms
  df$DistanceToGene[matched] <- mcols(nearest)$distance
  
  # ─── 5) Promoter/Intron/Exon annotation via ChIPseeker ──────────────────────
  peak_anno <- annotatePeak(
    gr,
    TxDb      = txdb,
    annoDb    = annoDb,
    tssRegion = tssRegion
  )
  anno_df <- as.data.frame(peak_anno) %>%
    mutate(
      GenomicCoords = paste0(seqnames, ":", start, "-", end),
      SYMBOL_unique = mapIds(orgdb,
                             keys      = geneId,
                             column    = "SYMBOL",
                             keytype   = "ENTREZID",
                             multiVals = "first")
    )
  
  df <- df %>%
    mutate(GenomicCoords = paste0(chr, ":", start, "-", end)) %>%
    left_join(
      anno_df %>%
        select(GenomicCoords, annotation, geneId, SYMBOL_unique, distanceToTSS),
      by = "GenomicCoords"
    )
  
  # ─── 6) Compute proportions ─────────────────────────────────────────────────
  print(head(df))
  prop_df <- df %>%
    filter(class %in% class_filter) %>%
    group_by(Cluster, class, GenicStatus) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Cluster, class) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
  # ─── 7) Build & save bar‐chart ───────────────────────────────────────────────
  plt <- ggplot(prop_df, aes(x = class, y = proportion, fill = GenicStatus)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ Cluster, nrow = 1) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x     = "Repeat Class",
      y     = "Proportion",
      fill  = "Genic Status",
      title = "Intergenic vs Intragenic Proportions by Class"
    ) +
    theme_minimal() +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    )
  
  out_file <- file.path(out_dir, paste0(prefix, "_genic_proportions.pdf"))
  ggsave(filename = out_file, plot = plt, width = 10, height = 6)
  message("→ Saved proportions plot to: ", out_file)
  
  # ─── 8) Compute family‐level proportions and plot ───────────────────────────
  ## <=== Highlighted section: family‐level proportions and barplot ===
  family_df <- df %>%
    group_by(Cluster, family, GenicStatus) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Cluster, family) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()

  family_plot <- ggplot(
    family_df,
    aes(x = family, y = proportion, fill = GenicStatus)
  ) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Cluster, scales = "free_x") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x     = "Repeat Family",
      y     = "Proportion",
      fill  = "Genic Status",
      title = "Intergenic vs Intragenic Proportions by Family"
    ) +
    theme_minimal() +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    )
  fam_plot_file <- file.path(out_dir, paste0(prefix, "_family_proportions.pdf"))
  ggsave(fam_plot_file, plot = family_plot, width = 12, height = 6)
  message("→ Saved family proportions plot to: ", fam_plot_file)


   # ─── 7) Save results ────────────────────────────────────────────────────────
  # Annotated repeats
  annotated_file <- file.path(out_dir, paste0(prefix, "_annotated_repeats.csv"))
  write_csv(df, annotated_file)
  message("→ Saved annotated repeats to: ", annotated_file)

  # Proportions table
  prop_file <- file.path(out_dir, paste0(prefix, "_genic_proportions.csv"))
  write_csv(prop_df, prop_file)
  message("→ Saved proportions table to: ", prop_file)

  message("done annotating repeat clusters")
  # ─── 8) Return results ───────────────────────────────────────────────────────
# ─── 9) Plot distanceToTSS histogram by cluster ──────────────────────────────
  hist_plot <- ggplot(df, aes(x = distanceToTSS)) +
    geom_histogram(bins = 100) +
    facet_wrap(~ Cluster, scales = "free") +
    labs(x = "Distance to TSS (bp)", y = "Count",
         title = "Distance to TSS Distribution by Cluster") +
    theme_minimal()
  hist_file <- file.path(out_dir, paste0(prefix, "_distanceToTSS_histograms.pdf"))
  ggsave(hist_file, plot = hist_plot, width = 12, height = 6)
  message("→ Saved distanceToTSS histograms to: ", hist_file)

  # ─── 10) Return results ─────────────────────────────────────────────────────
  invisible(list(
    annotated_df   = df,
    prop_df        = prop_df,
    class_plot     = plt,
    family_prop    = family_df,
    family_plot    = family_plot,
    distance_hist  = hist_plot
  ))
}


# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# library(org.Mm.eg.db)
# testE <- refactorCluster(km_res$cluster)


##############################################################################
##########################.  run enrichemnt analysis   ##########################
##############################################################################
run_enrichment_pipeline <- function(
  annotated_df,
  orgdb,
  out_dir = "enrichment_results",
  prefix  = "cluster",
  categories2Show = 10
) {
  # ─── Dependencies ───────────────────────────────────────────────────────────
  require(clusterProfiler)
  require(DOSE)
  require(AnnotationDbi)
  require(dplyr)
  require(ggplot2)
  require(readr)
  
  # ─── Prepare output folder ──────────────────────────────────────────────────
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # ─── 1. Prepare gene lists per cluster ─────────────────────────────────────
  cluster_genes <- annotated_df %>%
    filter(!is.na(SYMBOL_unique)) %>%
    distinct(Cluster, SYMBOL_unique) %>%
    split(.$Cluster)
  names(cluster_genes) <- paste0(prefix, names(cluster_genes))
  
  # Extract symbol vectors
  cluster_symbols <- lapply(cluster_genes, function(df) df$SYMBOL_unique)
  
  message("Map SYMBOL to ENTREZID")
  # ─── 2. Map SYMBOL to ENTREZID ─────────────────────────────────────────────
  cluster_entrez <- lapply(cluster_symbols, function(symbols) {
    mapIds(orgdb,
          gene          = genes,
            # org.Mm.eg.db
          # annoDb       = "org.Mm.eg.db",
           keys     = symbols,
           column   = "ENTREZID",
           keytype  = "SYMBOL",
           multiVals= "first") %>%
      na.omit() %>% unique()
  })
  
  # ─── 3. Enrichment loops ────────────────────────────────────────────────────
  go_results <- list(); kegg_results <- list()
  for (cl in names(cluster_entrez)) {
    genes <- cluster_entrez[[cl]]
    go_results[[cl]] <- list(
      BP = enrichGO(genes, OrgDb=orgdb, keyType="ENTREZID", ont="BP",
                    pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.10),
      MF = enrichGO(genes, OrgDb=orgdb, keyType="ENTREZID", ont="MF", 
                    pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.10),
      CC = enrichGO(genes, OrgDb=orgdb, keyType="ENTREZID", ont="CC", 
                    pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.10)
    )
    kegg_results[[cl]] <- enrichKEGG(genes, organism="mmu", keyType="kegg",
                                     pAdjustMethod="BH", pvalueCutoff=0.05)
  }
  
  # ─── 4. Combine GO results ─────────────────────────────────────────────────
  combined_go <- bind_rows(lapply(names(go_results), function(cl) {
    bind_rows(lapply(c("BP","MF","CC"), function(ont) {
      res <- go_results[[cl]][[ont]]@result
      if (nrow(res)==0) return(tibble(Cluster=cl,Ontology=ont))
      df <- as_tibble(res); df$Cluster<-cl; df$Ontology<-ont; df
    }))
  }))
  write_csv(combined_go, file.path(out_dir, paste0(prefix, "_GO_combined.csv")))
  
  # ─── 5. Plot dotplots ──────────────────────────────────────────────────────
  for (ont in c("BP","MF","CC")) {
    pdf(file.path(out_dir, paste0(prefix, "_GO_", ont, ".pdf")), width=10, height=6)
    for (cl in names(go_results)) {
      res <- go_results[[cl]][[ont]]
      if (nrow(res@result)>0) print(dotplot(res, showCategory=categories2Show)+ ggtitle(paste(cl,ont)))
    }
    dev.off()
  }
  
  # ─── 6. KEGG dotplots ───────────────────────────────────────────────────────
  pdf(file.path(out_dir, paste0(prefix, "_KEGG.pdf")), width=10, height=6)
  for (cl in names(kegg_results)) {
    res <- kegg_results[[cl]]
    if (nrow(res@result)>0) print(dotplot(res, showCategory=categories2Show)+ ggtitle(cl))
  }
  dev.off()
  
  # ─── 7. Return lists ───────────────────────────────────────────────────────
  invisible(list(
    cluster_entrez = cluster_entrez,
    go_results     = go_results,
    kegg_results   = kegg_results,
    combined_go    = combined_go
  ))
}






# head(refactorCluster(km_res$cluster))
res_kOptim <- annotate_repeat_clusters(
  refactorCluster(km_res$cluster),
  txdb   = TxDb.Mmusculus.UCSC.mm10.knownGene,
  orgdb  = org.Mm.eg.db,
  out_dir = "results_kOptim",
  prefix  = "kOptim"
)
# length(res_k2)
head(res_kOptim$annotated_df)


# k3_clustering_res$clusterDF
head(k3_clustering_res$kmeans$cluster)

res_k3 <- annotate_repeat_clusters(
  refactorCluster(k3_clustering_res$kmeans$cluster),
  txdb   = TxDb.Mmusculus.UCSC.mm10.knownGene,
  orgdb  = org.Mm.eg.db,
  out_dir = "results_k3",
  prefix  = "k3"
)



# dim(resultsDTSample)
splitRNA <- split(resultsDTSample, resultsDTSample$condition) 
mergedRNA_Annot <- lapply(names(splitRNA), function(x){
  merge(splitRNA[[x]], res_k2$annotated_df, by = "RepID", all.y = TRUE)
})


go_kOPtim <- run_enrichment_pipeline(res_kOptim[["annotated_df"]], orgdb=org.Mm.eg.db, out_dir = "enrichment_results_kOptim", prefix  = "cluster_kOptim_")

go_k3 <- run_enrichment_pipeline(res_k3[["annotated_df"]], orgdb=org.Mm.eg.db, out_dir = "enrichment_resultsK3", prefix  = "clusterK3_")

names(go_k3)
















































































































































































































































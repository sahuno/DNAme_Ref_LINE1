#!/usr/bin/env Rscript

# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/upsetPlots_DE.R
#!/usr/bin/env Rscript
Sys.setenv(TZ = "America/New_York")

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ComplexUpset)
  library(tidyverse)
})

#––– OPTPARSE SETUP –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
option_list <- list(
  make_option(c("-r", "--rna"),    type="character", help="Path to RNA TSV (required)", metavar="file"),
  make_option(c("-d", "--dna"),    type="character", help="Path to DNA TSV (required)", metavar="file"),
  make_option(c("-f", "--focus"),  type="character", default=NULL,
              help="Comma‐separated list of conditions to keep [default: all]"),
  # make_option(c("-l", "--fullLengthpath"),  type="character", default=NULL,
  #             help="Comma‐separated list of conditions to keep [default: all]")
  make_option(c("-c", "--classes"), type="character", default="LINE",
              help="Comma‐separated classLabels to process [default: %default]"),
  make_option(c("-m", "--minsets"), type="integer", default=2,
              help="Minimum number of sets a feature must appear in [default: %default]"),
  make_option(c("-o", "--prefix"),  type="character", default="upsetPlots",
              help="Output file prefix [default: %default]")
)

opt_parser <- OptionParser(
  usage = "Usage: %prog -r RNA_TSV -d DNA_TSV [options]",
  option_list = option_list
)
opt <- parse_args(opt_parser)

# enforce required args
if (is.null(opt$rna) || is.null(opt$dna)) {
  print_help(opt_parser)
  quit(status = 1)
}

focusConds  <- if (!is.null(opt$focus))  strsplit(opt$focus,  ",")[[1]] else NULL
classLabels <- strsplit(opt$classes, ",")[[1]]

#––– FUNCTION DEFINITIONS ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# 1) Read inputs and optionally filter by condition
read_inputs <- function(rna_path, dna_path, focusConds=NULL) {
  df_rna <- fread(rna_path, sep="\t", header=TRUE, data.table=FALSE)
  if (!is.null(focusConds)) {
    df_rna <- df_rna %>% filter(condition %in% focusConds)
  }
  df_dna <- fread(dna_path, data.table=FALSE)
  list(rna = df_rna, dna = df_dna)
}

# 2) Extract subfamily/family/class from RepID
preprocess_rna <- function(df_rna) {
  df_rna %>%
    select(-any_of(c("UpDown","sig"))) %>%
    mutate(
      ClassFamily = vapply(strsplit(RepID, "\\|"), `[`, character(1), 4)
    ) %>%
    separate(ClassFamily,
             into = c("subfamily","family","class"),
             sep  = ":",
             extra = "merge"
    )
}

# 3) Filter by up/down and classLabel
get_stats <- function(df, label = c("up","down"), classLabel) {
  label   <- match.arg(label)
  df_filt <- df %>%
    filter(LFClabel == label, grepl(classLabel, class))
  stats   <- df_filt %>%
    count(condition, name = "n") %>%
    arrange(desc(n))
  list(df = df_filt, stats = stats)
}

# 4) Build full + filtered membership matrices
build_membership <- function(df_filt, min_sets = 2L) {
  mat <- df_filt %>%
    distinct(RepID, condition) %>%
    mutate(present = 1L) %>%
    pivot_wider(
      names_from  = condition,
      values_from = present,
      values_fill = list(present = 0)
    ) %>%
    column_to_rownames("RepID")
  mat_filt <- mat[rowSums(mat) > min_sets, , drop = FALSE]
  list(full = mat, filtered = mat_filt)
}

# 5) Save a membership matrix to TSV
save_membership <- function(mat, path) {
  mat %>%
    rownames_to_column("RepID") %>%
    write_tsv(path)
}

# 6) Plot an UpSet using ComplexUpset (ggplot2) and save
plot_upset <- function(mat, path) {
  if (is.matrix(mat)) mat <- as.data.frame(mat)
  plt <- ComplexUpset::upset(
    mat,
    intersect       = colnames(mat),
    n_intersections = 20,
    name            = "Condition"
  )
  ggsave(path, plt, width = 10, height = 8)
}

#––– MAIN ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

main <- function(rna_path, dna_path, focusConds,
                 classLabels, min_sets, prefix) {
  inputs <- read_inputs(rna_path, dna_path, focusConds)
  df_rna <- preprocess_rna(inputs$rna)

  for (cls in classLabels) {
    for (lbl in c("up","down")) {
      tag    <- paste0(cls, "_", lbl)
      info   <- get_stats(df_rna, lbl, classLabel = cls)

      cat("\n--", toupper(cls), "/", toupper(lbl), "regulated counts:\n")
      print(info$stats)

      mem    <- build_membership(info$df, min_sets)
      full   <- mem$full
      filt   <- mem$filtered

      # Filenames embed class + direction
      mat_all_file  <- sprintf("%s_membership_%s.tsv",         prefix, tag)
      mat_filt_file <- sprintf("%s_membership_filtered_%s.tsv", prefix, tag)
      plot_file     <- sprintf("%s_%s_upset.pdf",             prefix, tag)

      save_membership(full,  mat_all_file)
      save_membership(filt,  mat_filt_file)

      cat("[", cls, "/", lbl, "] plotting dims:", dim(filt), "\n")
      plot_upset(filt, plot_file)
    }
  }
}

if (!interactive()) {
  main(
    rna_path    = opt$rna,
    dna_path    = opt$dna,
    focusConds  = focusConds,
    classLabels = classLabels,
    min_sets    = opt$minsets,
    prefix      = opt$prefix
  )
}

## Run for only L1
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/upsetPlots_DE.R \
#   --rna /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecific/resultsL1LocusSpecificRNA_DE_withMetadata.txt \
#   --dna /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/promoterMethyl/results/merged_repeatsDNAme_RNA_table_WithConds.tsv \
#   --focus CKi_DMSO,SETDB1i-CKi_DMSO,QSTAT-CKi_DMSO,QSTAT_DMSO,SETDB1i_DMSO \
#   --minsets 2 \
#   --prefix minsets2 \
#   --classes LINE


### Run for all repeats
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/upsetPlots_DE.R \
#   --rna /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecific/resultsAllRepeatsRNA_DE_withMetadata.txt \
#   --dna /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/promoterMethyl/results/merged_repeatsDNAme_RNA_table_WithConds.tsv \
#   --focus CKi_DMSO,SETDB1i-CKi_DMSO,QSTAT-CKi_DMSO,QSTAT_DMSO,SETDB1i_DMSO \
#   --minsets 2 \
#   --prefix minsets2 \
#   --classes SINE,LTR,LINE,tRNA,scRNA,DNA,Other,RC,snRNA,RNA

#rRNA,srpRNA
# "SINE"   "LTR"    "LINE"   "tRNA"   "scRNA"  "DNA"    "Other"  "RC"     "snRNA"   "srpRNA" "RNA"  

#   --rna /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecific/resultsAllRepeatsRNA_DE_withMetadata.txt  \


# pathFulllength <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/overlaps75SquireLINE_vs_L1_L1Base_wa.bed"

# colInterest <- "V11"
# fullDT <- fread(pathFulllength, select = colInterest)

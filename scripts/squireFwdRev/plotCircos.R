#!/usr/bin/env Rscript
# Usage: Rscript circos_plot.R --input <path/to/data.tsv> --outdir <output/directory> --keyword <plot_keyword>

# install.packages(c("optparse","circlize","dplyr","stringr"))  # if needed
# install.packages(c("optparse"))  # if needed
# batch 
# Usage: bash /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/batch_circos_runner.sh --input-dir /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/upsetPlots/separateHistoneModifications/QSTATi --outdir /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/circosPlots --script /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/plotCircos.R
# Usage: Rscript circos_plot.R --input <path/to/data.tsv> --outdir <output/directory> --keyword <plot_keyword>

# install.packages(c("optparse","circlize","dplyr","stringr"))  # if needed
# Usage: Rscript circos_plot.R --input <path/to/data.tsv> --outdir <output/directory> --keyword <plot_keyword> [--shared]

# install.packages(c("optparse","circlize","dplyr","stringr"))  # if needed

library(optparse)
library(dplyr)
library(stringr)
library(circlize)

# 1) Parse command line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to input TSV file", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help = "Directory to save output files", metavar = "character"),
  make_option(c("-k", "--keyword"), type = "character", default = "",
              help = "Plotting keyword to append to filename and title", metavar = "character"),
  make_option(c("-s", "--shared"), action = "store_true", default = FALSE,
              help = "If set, only keep rows where the sum across numeric columns equals 3")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate input
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Error: --input must be specified", call. = FALSE)
}

# Create output directory if it doesn't exist
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

# Construct output filename and title
base_name <- paste0("repeat_counts", if (nzchar(opt$keyword)) paste0("_", opt$keyword) else "", if (opt$shared) "_shared" else "")
out_file   <- file.path(opt$outdir, paste0(base_name, ".pdf"))
plot_title <- paste0("Repeat Element Counts per Chromosome", 
                     if (nzchar(opt$keyword)) paste0(" - ", opt$keyword) else "", 
                     if (opt$shared) " (shared sum=3)" else "")

# 2) Read & filter TSV (remove unplaced and random contigs)
message("Reading data from: ", opt$input)
df <- read.delim(opt$input, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(!str_detect(RepID, "chrUn_"), !str_detect(RepID, "_random")) %>%
  mutate(chrom = sub("\\|.*", "", RepID)) %>%
  filter(!str_detect(chrom, "random"))

# 4) Clean up column names and summarize by chromosome
names(df) <- gsub("\\.|-", "_", names(df))

# 3) Optionally keep only rows summing to 3
if (opt$shared) {
  message("Filtering rows where sum of numeric columns == 3")
  df <- df %>%
    rowwise() %>%
    mutate(row_sum = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    ungroup() %>%
    filter(row_sum == 3) %>%
    select(-row_sum)
  message(sprintf("Selected %d rows with total count == 3", nrow(df)))
  df <- df %>% 
  select(RepID, shared = 2, chrom)

  counts <- df %>%
  group_by(chrom) %>%
  summarize(across(where(is.numeric), sum), .groups = "drop")

# Check for data validity
conds <- names(counts)[-1]
if (nrow(counts) == 0) {
  stop(sprintf("No data to plot after filtering for '%s'.", opt$keyword))
}
if (length(conds) == 0) {
  stop(sprintf("No numeric columns to plot for '%s'.", opt$keyword))
}

# 5) Prepare plotting parameters
t <- counts[, conds, drop = FALSE]
max_val <- max(as.numeric(unlist(t)), na.rm = TRUE)
track_cols_master <- c("skyblue", "salmon", "lightgreen")
track_cols <- rep(track_cols_master, length.out = length(conds))

# 6) Plot to PDF
message("Saving plot to: ", out_file)
pdf(file = out_file, width = 10, height = 10)
circos.clear()
circos.par(cell.padding = c(0,0,0,0), start.degree = 90, gap.degree = 2)

# Initialize genome ideogram
circos.initializeWithIdeogram(species = "mm10")

for (i in seq_along(conds)) {
  cond <- conds[i]
  circos.trackPlotRegion(
    factors      = counts$chrom,
    y            = counts[[cond]],
    ylim         = c(0, max_val),
    track.height = 0.12,
    bg.border    = NA,
    panel.fun    = function(x, y) {
      this_chr <- get.cell.meta.data("sector.index")
      lims     <- get.cell.meta.data("xlim")
      h        <- counts[[cond]][counts$chrom == this_chr]
      circos.rect(xleft = lims[1], ybottom = 0,
                  xright = lims[2], ytop = h,
                  col = track_cols[i], border = NA)
    }
  )
}

# 7) Add main title
title(main = plot_title, cex.main = 1.2, line = -1)

# 8) Draw legend
par(xpd = TRUE)
legend("bottomleft", legend = conds, fill = track_cols,
       border = NA, bty = "n", inset = c(0.02, 0.02), title = "Condition")

# 9) Clean up
circos.clear()
dev.off()
message("Plot complete.")
}else {
  message("No filtering applied")
  counts <- df %>%
  group_by(chrom) %>%
  summarize(across(where(is.numeric), sum), .groups = "drop")

# Check for data validity
conds <- names(counts)[-1]
if (nrow(counts) == 0) {
  stop(sprintf("No data to plot after filtering for '%s'.", opt$keyword))
}
if (length(conds) == 0) {
  stop(sprintf("No numeric columns to plot for '%s'.", opt$keyword))
}

# 5) Prepare plotting parameters
t <- counts[, conds, drop = FALSE]
max_val <- max(as.numeric(unlist(t)), na.rm = TRUE)
track_cols_master <- c("skyblue", "salmon", "lightgreen")
track_cols <- rep(track_cols_master, length.out = length(conds))

# 6) Plot to PDF
message("Saving plot to: ", out_file)
pdf(file = out_file, width = 10, height = 10)
circos.clear()
circos.par(cell.padding = c(0,0,0,0), start.degree = 90, gap.degree = 2)

# Initialize genome ideogram
circos.initializeWithIdeogram(species = "mm10")

for (i in seq_along(conds)) {
  cond <- conds[i]
  circos.trackPlotRegion(
    factors      = counts$chrom,
    y            = counts[[cond]],
    ylim         = c(0, max_val),
    track.height = 0.12,
    bg.border    = NA,
    panel.fun    = function(x, y) {
      this_chr <- get.cell.meta.data("sector.index")
      lims     <- get.cell.meta.data("xlim")
      h        <- counts[[cond]][counts$chrom == this_chr]
      circos.rect(xleft = lims[1], ybottom = 0,
                  xright = lims[2], ytop = h,
                  col = track_cols[i], border = NA)
    }
  )
}

# 7) Add main title
title(main = plot_title, cex.main = 1.2, line = -1)

# 8) Draw legend
par(xpd = TRUE)
legend("bottomleft", legend = conds, fill = track_cols,
       border = NA, bty = "n", inset = c(0.02, 0.02), title = "Condition")

# 9) Clean up
circos.clear()
dev.off()
message("Plot complete.")
}







# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/plotCircos.R \
# --input /Users/ahunos/upsetPlots/separateHistoneModifications/QSTATi/minsets1FullLength_membership_LINE_up.tsv \
# --outdir /Users/ahunos/Downloads/ \
# --keyword LINEs_QSTATi

# opt$input <- "/Users/ahunos/upsetPlots/separateHistoneModifications/QSTATi/minsets1FullLength_membership_LINE_up.tsv"
# opt$outdir <- "/Users/ahunos/Downloads/"
# ../upsetPlots/separateHistoneModifications/QSTATi/
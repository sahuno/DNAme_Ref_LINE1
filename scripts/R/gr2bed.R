# convert to bedfile

#########################
### convert to bedfile
# Convert GRanges to a data frame
# sorted_gr <- sort(l1Base_entireLength_gr)
counts_gr_L1 <- sort(counts_gr_L1)

bed_df <- data.frame(
  chrom = seqnames(counts_gr_L1),       # Chromosome
  start = start(counts_gr_L1) - 1,      # BED format uses 0-based start
  end = end(counts_gr_L1),              # BED format uses 1-based end
  name = counts_gr_L1$te.id,            # Name column (can be any metadata)
  score = ".",                          # Placeholder for BED score
  strand = strand(counts_gr_L1)         # Strand information
)
# Optional: Save extra metadata columns in a separate BED-like file
metadata_cols <- data.frame(mcols(counts_gr_L1))

# Combine metadata if needed
bed_df <- cbind(bed_df, metadata_cols)

bed_df$start <- format(bed_df$start, scientific = FALSE, trim = TRUE)
bed_df$end <- format(bed_df$end, scientific = FALSE, trim = TRUE)

write.table(
  bed_df[, 1:6],  # Standard BED format uses only first 6 columns
  file = "SquireL1.bed",
  quote = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = FALSE
)



###################################
##### step 2 #########################
###################################
# Convert GRanges to a data frame in BED format (BED6+ with extra metadata)
l1Base_entireLength_gr <- sort(l1Base_entireLength_gr)

bed_df <- data.frame(
  chrom = seqnames(l1Base_entireLength_gr),  # Chromosome
  start = start(l1Base_entireLength_gr) - 1, # BED format is 0-based start
  end = end(l1Base_entireLength_gr),        # BED format is 1-based end
  name = l1Base_entireLength_gr$RepeatID,   # Repeat ID as name
  score = l1Base_entireLength_gr$score,     # Score (integer)
  strand = strand(l1Base_entireLength_gr),  # Strand
  thickStart = l1Base_entireLength_gr$start2 - 1, # Thick start (optional, for UCSC tracks)
  thickEnd = l1Base_entireLength_gr$end2,        # Thick end (optional)
  color = l1Base_entireLength_gr$color,     # RGB color (for UCSC tracks)
  width = l1Base_entireLength_gr$SeqWidth_L1Base # Additional metadata
)


bed_df$start <- format(bed_df$start, scientific = FALSE, trim = TRUE)
bed_df$end <- format(bed_df$end, scientific = FALSE, trim = TRUE)


write.table(
  bed_df,
  file = "l1Base_entireLength.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

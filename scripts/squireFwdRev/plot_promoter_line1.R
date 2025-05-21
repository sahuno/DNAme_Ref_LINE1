# fread("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/results_LINE1_promoters/results/merged_all_samples.bed")

library(tidyverse)
library(data.table)
library(fst)


library(ggplot2)
library(dplyr)
library(ggrepel)    # for geom_text_repel()

# grep -E "chr15|63776172|63776324|Lx8:L1:LINE|285|[+]" /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/results_LINE1_promoters/results/first_300bp/merged_first_300bp_samples.bed

DNAme_first100 <- fread("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/results_LINE1_promoters/results/first_300bp/merged_first_300bp_samples.bed")
rna_DF <- fread("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/resultsDTSampleUpDown_with_chromCord.txt")
fpkm_path <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/squireRepeatsAnnot_Aggregated_FWD_REV_fpkm.fst"


DNAme_first100 %>% filter(str_detect(V1, "D-Q-")) %>% filter(str_detect(V5,"chr15|63776172|63776324|Lx8:L1:LINE|285|[+]"))





dim(DNAme_first100)

rna_DF_L1 <- rna_DF %>% filter(LFClabel %in% c("up", "down")) %>% filter(family == "L1")
rna_DF_L1_ids <- rna_DF_L1 %>% select(c(RepID, LFClabel, condition))

# rna_DF_L1  %>% filter(str_detect(condition, "QSTAT-CKi_DMSO"))%>% group_by()  %>% summarise(n=n()) %>% filter(n > 1) %>% arrange(desc(n))

DNAme_first100_QSTATonly <- DNAme_first100 %>% filter(str_detect(V1, "D-QC|D-0|D-C|D-Q")) #%>% filter(V5 %in% uniqueRepeatsL1)

table(DNAme_first100_QSTATonly$V1)
#DMSO only
DNAme_first100_QSTATonly <- DNAme_first100_QSTATonly %>% filter(!V1 %in% c("D-0-1_5000", "D-0-1_4000", "D-0-2_5000", "D-0-2_4000", "D-Q-1_5000", "D-Q-1_4000"))

drop_lowQualSample <- c("D-C-3")

# QC_DF <- rna_DF %>% 
# filter(condition=="QSTAT-CKi_DMSO", family == "L1")   %>%  mutate(RepID = paste0(RepID,"|",isfullLengthL1)) 

# uniqueRepeatsL1 <- unique(QC_DF$RepID)
1

# splitted_list <- split(
#   as.data.frame(DNAme_first100_QSTATonly),
#   DNAme_first100_QSTATonly$V1
# )

# DNAme_first100_QSTATonly %>% filter(V5 == "chr1|7180618|7185932|L1_Mus1:L1:LINE|93|+|FALSE")

# QC_DNA_RNA_DF <- QC_DF %>% left_join(DNAme_first100_QSTATonly , by = c("RepID" = "V5"))


# library(dplyr)

# QC_DNA_RNA_DF <- QC_DNA_RNA_DF %>%
#   mutate(across(c(V8, V9),
#                 ~ na_if(.x, ".") %>% as.numeric()))

# # check
# QC_DNA_RNA_DF %>% summarise(across(c(V8, V9), ~ class(.)))



#   pltDMSOQSTAT<- QC_DNA_RNA_DF %>% tidyplot(x = V8, y = log2FoldChange, color = LFClabel) |> 
#     add_data_points(alpha = 0.4) |> 
#     # add_test_asterisks(hide.ns = FALSE, method = "wilcox_test", hide_info = TRUE) |>
#   adjust_x_axis(rotate_labels = 45)|>
# #   split_plot(by = class, ncol = 2, nrow = 2) |>
# #   save_plot(filename = "RepeatsMajorClassesUp_down_LFC_QSTATi_Group.pdf", view_plot = FALSE) |>
#   split_plot(by = V1) |>
#     save_plot(filename = "L1_LFC_DNAme_expression.pdf", view_plot = FALSE, multiple_files = FALSE)



#   pltDMSOQSTAT_minCpGs <- QC_DNA_RNA_DF  %>% filter(V10 >= 3) %>% tidyplot(x = V8, y = log2FoldChange, color = LFClabel) |> 
#     add_data_points(alpha = 0.4, mapping = aes(shape = isfullLengthL1)) |> 
#     # add_test_asterisks(hide.ns = FALSE, method = "wilcox_test", hide_info = TRUE) |>
#   adjust_x_axis(rotate_labels = 45)|>
# #   split_plot(by = class, ncol = 2, nrow = 2) |>
# #   save_plot(filename = "RepeatsMajorClassesUp_down_LFC_QSTATi_Group.pdf", view_plot = FALSE) |>
#   split_plot(by = V1) |>
#     save_plot(filename = "L1_LFC_DNAme_expressionMinCpGs_3.pdf", view_plot = FALSE, multiple_files = FALSE)






# library(ggplot2)

# plt <- ggplot(
#   QC_DNA_RNA_DF %>% filter(V10 >= 3, !V1 %in% c("D-0-1_5000", "D-0-1_4000", "D-0-2_5000", "D-0-2_4000")),
#   aes(
#     x     = V8,
#     y     = log2FoldChange,
#     color = LFClabel,
#     shape = isfullLengthL1    # map shape here
#   )
# ) +
#   geom_point(alpha = 0.7, size = 3) +
#   facet_wrap(~ V1, scales = "free_x", ncol = 2) +
#   scale_shape_manual(
#     values = c(
#       "FALSE" = 16,  # solid circle
#       "TRUE"  = 17   # solid triangle
#     )
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.text  = element_text(face = "bold")
#   ) +
#   labs(
#     x     = "mean methylation level",
#     y     = "log2 Fold Change",
#     shape = "Full-Length L1?",
#     color = "LFC Label"
#   )
# ggsave(
#   filename = "L1_LFC_DNAme_expressionMinCpGs_3_ggplot2.pdf",
#   plot     = plt,
#   width    = 10,
#   height   =12
# )


# pltHistogram <- ggplot(
#   QC_DNA_RNA_DF %>% filter(V10 >= 3, !V1 %in% c("D-0-1_5000", "D-0-1_4000", "D-0-2_5000", "D-0-2_4000")),
#   aes(
#     x     = V8,
#     fill = LFClabel,
#     # shape = isfullLengthL1    # map shape here
#   )
# ) + geom_density() +
#   facet_wrap(~ V1, scales = "free_x", ncol = 2) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.text  = element_text(face = "bold")
#   ) +
#   labs(
#     x     = "mean methylation level",
#     y     = "Density",
#     shape = "Full-Length L1?",
#     color = "LFC Label"
#   )
# ggsave(
#   filename = "L1_LFC_DNAme_expressionMinCpGs_3_histogram.pdf",
#   plot     = pltHistogram,
#   width    = 10,
#   height   =12
# )
#   ungroup() -> DNAme_first100_QSTATonly
# split(V1) -> DNAme_first100_QSTATonly


# genomeChromSizes=/data1/greenbab/database/mm10/mm10.sorted.chrom.sizes
# genomeFasta=/data1/greenbab/database/mm10/mm10.fa
# LINE1_promoters=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LINE1_promoters/first_100bp.bed




##########################################

colNames <- c("TE_chr", "TE_start", "TE_stop","TE_ID","validReadCounts", "R-0-1", "R-0-2", "R-0-3", "R-C-1", "R-C-2", "R-C-3", "R-Q-1",  "R-Q-2", "R-Q-3", "R-QC-1", "R-QC-2", "R-QC-3")
dt_fpkm <- read_fst(fpkm_path, columns = colNames)
head(dt_fpkm)

# table(dt_fpkm$)
# DNAme_first100_QSTATonly %>% separate_wider_delim(V5)

DNAme_first100_QSTATonly <- DNAme_first100_QSTATonly %>%
  separate(V5, into = c("TE_ID", "isFullLength"), sep = "\\|(?=[^|]+$)")

table(DNAme_first100_QSTATonly$V1)

#get unique line 1s
TE_ID_DNAme <- unique(DNAme_first100_QSTATonly$TE_ID)
 
head(DNAme_first100_QSTATonly)

# table(dt_fpkm$TE_ID)
dt_fpkm_uniqueTEs <- dt_fpkm %>% filter(TE_ID %in% TE_ID_DNAme)
head(dt_fpkm_uniqueTEs)

DNAme_first100_QSTATonly <- DNAme_first100_QSTATonly %>% mutate(samples = gsub("_.*","", V1)) #%>% head()
head(DNAme_first100_QSTATonly)
table(DNAme_first100_QSTATonly$samples)

colnames(dt_fpkm_uniqueTEs) <- sub("^R-", "D-", colnames(dt_fpkm_uniqueTEs))
head(dt_fpkm_uniqueTEs)

head(dt_fpkm_uniqueTEs)
dt_fpkm_long <- dt_fpkm_uniqueTEs %>%
  pivot_longer(
    cols = matches("^D-"),
    names_to = "sample",
    values_to = "value"
  )

table(dt_fpkm_long$sample)
table(DNAme_first100_QSTATonly$samples)
head(DNAme_first100_QSTATonly)
head(dt_fpkm_long)

dt_fpkm_long %>% filter(str_detect(sample, "D-Q-")) %>% filter(TE_ID == "chr15|63776172|63776324|Lx8:L1:LINE|285|+")
DNAme_first100_QSTATonly %>% filter(str_detect(samples, "D-Q-")) %>% filter(TE_ID == "chr15|63776172|63776324|Lx8:L1:LINE|285|+")


dt_fpkm_long %>% left_join(DNAme_first100_QSTATonly, by = c("TE_ID" = "TE_ID", "sample" = "samples")) %>% 
  filter(!is.na(value)) %>% 
  select(TE_ID, sample, isFullLength, value, V8, V9) %>% 
  mutate(across(c(V8, V9), ~ na_if(.x, ".") %>% as.numeric())) -> dt_fpkm_long_DNAme

head(dt_fpkm_long_DNAme)


## drop low quality RNA samples
dt_fpkm_long_DNAme <- dt_fpkm_long_DNAme %>% filter(!str_detect(sample, drop_lowQualSample))

table(dt_fpkm_long_DNAme$sample)
# dt_fpkm_long_DNAme %>% filter(str_detect(sample, "D-QC-"))



############## sanity checks ################
rna_DF_L1_ids %>% filter(condition %in% "CKi_DMSO")
rna_DF_L1_ids %>% filter(condition %in% "QSTAT_DMSO")
rna_DF_L1_ids %>% filter(condition %in% "QSTAT-CKi_DMSO")
######################## sanity checks ################


table(rna_DF_L1_ids$condition)
cki_rna_dname <- dt_fpkm_long_DNAme %>% filter(str_detect(sample, "D-C-")) %>% left_join(rna_DF_L1_ids %>% filter(condition %in% "CKi_DMSO"), by = c("TE_ID" = "RepID"))
qc_rna_dname <- dt_fpkm_long_DNAme %>% filter(str_detect(sample, "D-QC-")) %>% left_join(rna_DF_L1_ids %>% filter(condition %in% "QSTAT-CKi_DMSO"), by = c("TE_ID" = "RepID"))
q_rna_dname <- dt_fpkm_long_DNAme %>% filter(str_detect(sample, "D-Q-")) %>% left_join(rna_DF_L1_ids %>% filter(condition %in% "QSTAT_DMSO"), by = c("TE_ID" = "RepID"))

dt_fpkm_long_DNAme %>% filter(str_detect(sample, "D-Q-")) %>% filter(TE_ID == "chr15|63776172|63776324|Lx8:L1:LINE|285|+")

dmso_rna_dname <- dt_fpkm_long_DNAme %>% filter(str_detect(sample, "D-0-")) %>% mutate(LFClabel=NA ,condition=NA)



merged_rna_dname <- rbind(
  cki_rna_dname,
  qc_rna_dname,
  q_rna_dname,
  dmso_rna_dname
)
dim(merged_rna_dname)

# dt_fpkm_long_DNAme %>% filter(str_detect(sample, "D-Q-")) %>% group_by(sample) %>% summarise(mean = mean(value, na.rm = TRUE)) %>% arrange(desc(mean)) 
dt_fpkm_long_DNAme %>% filter(str_detect(sample, "D-Q-"))%>% group_by()  %>% summarise(n=n()) %>% filter(n > 1) %>% arrange(desc(n))
dim(dt_fpkm_long_DNAme)

rna_DF_L1_ids %>% filter(str_detect(condition, "QSTAT-CKi_DMSO")) %>% pull(RepID) %>% head()
dt_fpkm_long_DNAme %>% filter(str_detect(sample, "D-QC-")) %>% pull(TE_ID) %>% head()
# dt_fpkm_long_DNAme %>% filter(str_detect(sample, "D-Q-")) %>% group_by(sample) %>% summarise(mean = mean(value, na.rm = TRUE)) %>% arrange(desc(mean)) 

# %>% filter(V10 >= 3)




merged_rna_dname %>% filter(TE_ID == "chr15|63776172|63776324|Lx8:L1:LINE|285|+") #%>% table()


df2 <- merged_rna_dname %>%
  filter(!is.na(condition), !is.na(LFClabel))

pltBox <- ggplot(df2, aes(x = sample, y = V8, fill = LFClabel)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  facet_wrap(~ condition, scales = "free_x", nrow = 1) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x      = "Sample",
    y      = "Mean DNAme",
    fill   = "LFC Label",
    title  = "Distribution of L1 mean DNAme"
  )

ggsave(
  filename = "L1_promoter_DNAme_boxplot.pdf",
  plot     = pltBox,
  width    = 10,
  height   = 12
)



#################### ploting 
plt <- ggplot(
  merged_rna_dname,
  aes(
    x     = V8,
    y     = value,
    color = LFClabel,
    shape = isFullLength    # map shape here
  )
) +
  geom_point(alpha = 0.7, size = 3) +
  facet_wrap(~ sample, scales = "free_x", ncol = 3) +
  # scale_shape_manual(
  #   values = c(
  #     "FALSE" = 16,  # solid circle
  #     "TRUE"  = 17   # solid triangle
  #   )) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text  = element_text(face = "bold")
  ) +
  labs(
    x     = "mean methylation level",
    y     = "FPKM",
    shape = "isFullLength",
    color = "LFC Label",
    title = "methylation level vs FPKM of LINE1"
  )
ggsave(
  filename = "L1_fpkm_DNAme_expressionMinCpGs_3_ggplot2_first300bp_highQualSamples.pdf",
  plot     = plt,
  width    = 10,
  height   =12
)




# 1) find top 5 per sample
top5 <- merged_rna_dname %>%
  group_by(sample) %>%
  slice_max(order_by = value, n = 3, with_ties = FALSE)

# 2) build your plot, adding geom_text_repel for just those points
plt_top <- ggplot(
  merged_rna_dname,
  aes(
    x     = V8,
    y     = value,
    color = LFClabel,
    shape = isFullLength    # map shape here
  )
) +
  geom_point(alpha = 0.7, size = 3) +
  geom_text_repel(
    data = top5,
    aes(label = TE_ID),
    size = 3,
    max.overlaps = 10
  ) +
  facet_wrap(~ sample, scales = "free_x", ncol = 3) +
  # scale_shape_manual(
  #   values = c("FALSE" = 16, "TRUE" = 17)
  # ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text  = element_text(face = "bold")
  ) +
  labs(
    x     = "mean methylation level",
    y     = "FPKM",
    shape = "Full-Length L1?",
    color = "LFC Label",
    title = "methylation level vs FPKM of LINE1"
  )

# 3) save
ggsave(
  "L1_fpkm_DNAme_expressionMinCpGs_3_ggplot2_top5labels_first300_highQual.pdf",
  plot   = plt_top,
  width  = 10,
  height = 12
)

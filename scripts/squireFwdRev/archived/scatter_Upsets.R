library(UpSetR)
library(dplyr)
library(tidyr)
library(data.table)
library(tidyverse)
library(ComplexUpset)
##
pathsDNAme <- "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/promoterMethyl/results/merged_repeatsDNAme_RNA_table_WithConds.tsv"
# 1) Read the data

df <- read.delim("resultsRNA_DE_withMetadata.txt", 
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dnarna <- fread(pathsDNAme) 

# stringr::str_split_i(head(df$RepID), "\\|", 4)


# names(df)
dfSelect <- df %>% select(!c(UpDown, sig)) #%>% head()

dfSelect <- dfSelect %>% mutate(ClassFamily = str_split_i(RepID, "\\|", 4)) %>%
separate(ClassFamily, into = c("subfamily", "family", "class"), sep = ":")

table(dfSelect$class)


slitDE <- split(dfSelect, dfSelect$condition)

#dnarna %>% filter(str_detect(condition, gsub("_","|","AZA_DMSO"))) %>% 

dfSelect %>% mutate(ClassFamily = str_split_i(RepID, "\\|", 4)) %>%
separate(ClassFamily, into = c("subfamily", "family", "class"), sep = ":")




# 2) Filter for upregulated LINE elements
df_up <- df %>% 
  filter(LFClabel == "up", grepl(":LINE", RepID))

df_upStats <- df_up %>% group_by(condition) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) # Check the number of upregulated LINE elements per condition

didDF <- df_up %>%
  select(RepID, condition)

# 3) Build a binary membership matrix: rows are RepID, columns are conditions
membership <- df_up %>%
  select(RepID, condition) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from  = condition,
    values_from = present,
    values_fill = list(present = 0)
  ) %>%
  column_to_rownames("RepID")

rowSums(membership) # check the number of conditions per RepID


membership_with_rownames <- membership %>%
  rownames_to_column(var = "RepID")
# Save the membership matrix to a file
write.table(membership_with_rownames, "membership_matrix_upReg.tsv")

# Check if the sum of conditions per RepID is greater than 5
membership_filtered <- membership[rowSums(membership) > 3, ]
rowSums(membership_filtered) # Display the filtered sums

# Save the filtered membership matrix to a file
membership_filtered_with_rownames <- membership_filtered %>%
  rownames_to_column(var = "RepID")
write_tsv(membership_filtered_with_rownames, "membership_filtered_matrix_upReg.tsv")

# 4) Draw the UpSet plot
#    - nsets = number of conditions to show
#    - nintersects = number of intersections to plot (NA = all)
png("upRegulatedLINE_upsetPlot.png", width = 900, height = 600, res = 100)
upset(
  membership_filtered,
  nsets       = ncol(membership_filtered),
  nintersects = NA,
  order.by    = "freq"
)
dev.off()



png("DownRegulatedLINE_upsetPlot.png", width = 900, height = 600, res = 100)
upset(
  filt,
  nsets       = ncol(filt),
  nintersects = NA,
  order.by    = "freq"
)
dev.off()


filt_filtered <- filt[rowSums(filt) > 3, ]

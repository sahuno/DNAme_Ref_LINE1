#compile somatic inertions for line1
library(tidyverse)
library(data.table)

options(width=250)

#read tldr outputs
dir_tldr_table <- '/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/LINE1_insertions_tldr/outputs'
paths_tldr_table <- list.files(dir_tldr_table, pattern=".table.txt", recursive = TRUE, full.names=TRUE)

#load mouse metadata
mouseTriEpi_metadata <- read_csv("/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_recoded.csv")
mouseTriEpi_metadata$samples <- gsub("R-","D-", mouseTriEpi_metadata$new_samples_name)
mouseTriEpi_metadata$condition <- gsub("AZA","5-AZA", mouseTriEpi_metadata$condition)

removeSamplesOut <- c("D-0-1_4000.table.txt", "D-0-1_5000.table.txt","D-0-2_4000.table.txt", "D-0-2_5000.table.txt", "D-Q-1_5000.table.txt", "D-Q-1_4000.table.txt", "D-S-1_5000.table.txt", "D-S-1_4000.table.txt")
paths_tldr_table <- paths_tldr_table[!grepl(paste(removeSamplesOut, collapse = "|"),paths_tldr_table)]

sampleNames_tldr <- gsub(".table.txt","", basename(paths_tldr_table))
sampleNames_tldr <- gsub("_.*","", sampleNames_tldr)

# sort(gsub("_.*","", sampleNames_tldr))


list_tldr_table <- lapply(paths_tldr_table, function(x){fread(x)})
names(list_tldr_table) <- sampleNames_tldr
lapply(list_tldr_table, dim)

dt_tldr <- rbindlist(list_tldr_table, idcol = "samples")
names(dt_tldr)
table(dt_tldr$samples)
#how many L1 insertions per sample
dt_tldr %>% group_by(samples) %>% filter(Family == "L1" &Filter == "PASS") %>% summarise(n())
Insertions_All_stats_df <- dt_tldr %>% group_by(samples, Family) %>% filter(Filter == "PASS" & NumSamples == 1) %>% summarise(n())
Insertions_All_stats_df <- Insertions_All_stats_df %>% left_join(mouseTriEpi_metadata)

Insertions_all <- ggplot(Insertions_All_stats_df, aes(x = samples, y = `n()`, fill = Family)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Insertions per sample & conditions",
       x = "samples",
       y = "Count",
       fill = "Family") + facet_wrap(~condition, scales = "free_x") #+
  #theme_minimal()#
ggsave(Insertions_all, filename="figures/nInsertionsTRI_Epi.png")


#there are 8 sub-families
# tldr_sam_L1 %>% group_by(Subfamily) %>% summarise(n())
# tldr_sasha %>% group_by(Subfamily) %>% summarise(n())

# unique(tldr_sam_L1$Filter) # i have more here
# unique(tldr_sasha$Filter) #oh you have only subsetted only `passed` filter 

# tldr_sam_L1 %>% group_by(Filter)

#observations; these are unpased, refernce only(NonRef column), 
#now why is not good idea to find methylation of `Chrom:Start-End` ?
#set minimum reads `UsedReads`? or?
# what's `NumSamples`? i used 2 samples but col has 1 in some entries
# head(tldr_sasha)
# head(tldr_sam)[,"SampleReads"]
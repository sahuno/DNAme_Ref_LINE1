#compile somatic inertions for line1
library(tidyverse)
library(data.table)

options(width=250)
sOV051_TN <- "/data1/greenbab/projects/methyl_benchmark_spectrum/ONT_BSseq/ONT_DLP_1stPre/tldr_outputs/results/Spectrum-OV-051_TN/Spectrum-OV-051_T_final_Spectrum-OV-051_N_final.table.txt"
s009T2_009T1_n_2_N <- "/data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/L1_insertions_output/s009TN/009T2_009T1_n_2_modBaseCalls_sorted_dedup_009N_modBaseCalls_sorted_dup.table.txt"
s044T_v14 <- '/data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/L1_insertions_output/044T_v14_modBaseCalls_sorted_dup_044N_v14_modBaseCalls_sorted_dup.table.txt'

paths_list_tldr_table <- c(sOV051_TN,s009T2_009T1_n_2_N,s044T_v14)

list_tldr_table <- lapply(paths_list_tldr_table, function(x){fread(x)})
lapply(list_tldr_table, dim)
# c("sOV051_TN", "s009T2_009T1_n_2_N","s044T_v14")
dt_tldr <- rbindlist(list_tldr_table, idcol = TRUE)
names(dt_tldr)
table(dt_tldr$.id)
dt_tldr <- dt_tldr %>% mutate(`.id` = case_when(`.id` == 1 ~ "sOV051_TN", `.id` == 2 ~ "s009T2_009T1_n_2_N", `.id` == 3 ~ "s044T_v14"))

#i have more line1 entries than you.
L1_pre_stats <- dt_tldr %>% group_by(.id) %>% filter(!is.na(Family) & Family == "L1" & Filter == "PASS") %>% summarise(n()) #filter passed
L1_pre_filtered <- dt_tldr %>% group_by(.id) %>% filter(!is.na(Family) & Family == "L1" & Filter == "PASS") 
L1_Nsamples1_filtered <- L1_pre_filtered %>% filter(NumSamples == 1) 

table(L1_Nsamples1_filtered$.id)

#filters, numb samples = 1, 1 read with non-zero quality should be good

#plot the endTE
plt <- ggplot(L1_Nsamples1_filtered, aes(EndTE)) + geom_histogram()
ggsave(plt, filename="EndTE_spectrum.png")


#start can be anything. less
#endTE should be greater than 5900
L1_Nsamples1_FullLength_filtered <- L1_Nsamples1_filtered %>% filter(EndTE >= 5900) 
dim(L1_Nsamples1_FullLength_filtered) #getting close
table(L1_Nsamples1_FullLength_filtered$Chrom)
#remove unalignme
L1_Nsamples1_FullLength_filtered_standChrom <- L1_Nsamples1_FullLength_filtered %>% filter(!str_detect(Chrom, "_"))
head(L1_Nsamples1_FullLength_filtered_standChrom)[, c("Chrom" , "Start", "End", "StartTE", "EndTE")] %>% mutate(width_ref = End - Start)
dim(L1_Nsamples1_FullLength_filtered_standChrom)

#how many invversions
df_invr <- L1_Nsamples1_FullLength_filtered_standChrom %>% group_by(.id, Inversion) %>% summarise(n()) #filter passed
#plot full length insertions
plt_nFullengths <- ggplot(L1_Nsamples1_FullLength_filtered_standChrom, aes(`.id`)) + geom_bar()
# plt_InsertionLengths <- ggplot(L1_Nsamples1_FullLength_filtered_standChrom, aes(y=LengthIns, x=`.id`)) + geom_col()

plt_InsertionLengths <- ggplot(L1_Nsamples1_FullLength_filtered_standChrom, aes(y=LengthIns, x=`.id`)) + geom_boxplot()

ggsave(plt_nFullengths, filename="nFullLength_InsertionsSpectrum.png")
ggsave(plt_InsertionLengths, filename="nFullLength_InserLengthsSpectrum.png")

pInversions <- ggplot(df_invr, aes(x = .id, y = `n()`, fill = Inversion)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Inversion Counts per ID",
       x = "ID",
       y = "Count",
       fill = "Inversion") +
  theme_minimal()
ggsave(pInversions, filename="nFullLength_InversionsSpectrum.png")

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
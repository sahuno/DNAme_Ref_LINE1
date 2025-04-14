# l1 EM annoations
# system("wget https://www.repeatmasker.org/L1/EMBL/EMBL_L1_annotated.bed.gz")


library(data.table)
library(tidyverse)

D <- fread("mm39.L1EM.bed")

D %>% separate_wider_delim(delim = ".", cols = V4, names = c("name", "score", "chrom", "start", "end", "strand"), too_few = "debug") %>% select(name) %>% group_by(name) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n))

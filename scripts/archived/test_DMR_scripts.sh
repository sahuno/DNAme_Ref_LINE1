#test how dmr is working

modkit dmr multi \
  -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/results/DNAme_overlaps/D-0-1_5000_4000/D-0-1_5000_4000_5mCpG_5hmCpG_sortedBed_minCov20_CpGIslands.bed.gz ctrl1 \
  -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/results/DNAme_overlaps/D-S-1_5000_4000/D-S-1_5000_4000_5mCpG_5hmCpG_sortedBed_minCov20_CpGIslands.bed.gz setdb11 \
  -o dmr_dir_D01DS1 \
  -r /data1/greenbab/database/mm10/mm10_CpGIslands.bed \
  --ref /data1/greenbab/database/mm10/mm10.fa \
  --base C \
  -t 10 \
  -f \
  --log-filepath dmr_multi_D01DS1.log


# ssend modkit dmr /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/test_DMR_scripts.sh


library(data.table)
df <- fread("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/dmr_dir_D01DS1/ctrl1_setdb11.bed")
head(df)
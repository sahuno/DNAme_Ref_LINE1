
# 1. Get the minimal L1 SVs
l1_minimal=/data1/greenbab/projects/methyl_benchmark_spectrum/ONT_BSseq/ONT_DLP_1stPre/L1_DNAme/LINE1_SVs_minimal.bed

#
for f in /data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/hac_5mCG_5hmCG/results/modkit/044T_v14/044T_v14_modpileup_5mC.bed \
/data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/sup_5mCG_5hmCG/results/modkit/044N_v14/044N_v14_modpileup_5mC.bed \
/data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/hac_5mCG_5hmCG/methyl_modkit_mergedData/009T2_009T1_n_2_modBaseCalls_sorted_dedup_5mC.bed \
/data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/009N_rerun/modbasecalls/mergebam_modkit/results/modkit/009N/009N_modpileup_5mC.bed
do
NAME=$(basename -s .bed ${f})
awk 'BEGIN {FS=OFS="\t"} $10 > 15' $f > ${NAME}_minCov15.bed
sort -k1,1 -k2,2n ${NAME}_minCov15.bed > ${NAME}_minCov15.sorted.bed
# bedtools map -a ${l1_minimal} -b ${NAME}_minCov15.sorted.bed -c 11,11,11 -o mean,median,count > ${1}.methylation.bed
done


bedMethylMinCov=/data1/greenbab/projects/methyl_benchmark_spectrum/ONT_BSseq/ONT_DLP_1stPre/L1_DNAme/overlaps_L1_methylation
ROOT_L1_SVs=/data1/greenbab/projects/methyl_benchmark_spectrum/ONT_BSseq/ONT_DLP_1stPre/L1_DNAme
bedtools map -a ${ROOT_L1_SVs}/sample_044T.bed -b ${bedMethylMinCov}/044T_v14_modpileup_5mC_minCov15.sorted.bed -c 11,11,11 -o mean,median,count > 044T.L1_SVs_DNAme.bed
bedtools map -a ${ROOT_L1_SVs}/sample_044N.bed -b ${bedMethylMinCov}/044N_v14_modpileup_5mC_minCov15.sorted.bed -c 11,11,11 -o mean,median,count > 044N.L1_SVs_DNAme.bed
bedtools map -a ${ROOT_L1_SVs}/sample_009T.bed -b ${bedMethylMinCov}/009T2_009T1_n_2_modBaseCalls_sorted_dedup_5mC_minCov15.sorted.bed -c 11,11,11 -o mean,median,count > 009T.L1_SVs_DNAme.bed
bedtools map -a ${ROOT_L1_SVs}/sample_009N.bed -b ${bedMethylMinCov}/009N_modpileup_5mC_minCov15.sorted.bed -c 11,11,11 -o mean,median,count > 009N.L1_SVs_DNAme.bed



# REF=/data1/greenbab/database/hg38/v0/Homo_sapiens_assembly38.fasta
# L1_bed=/data1/greenbab/projects/methyl_benchmark_spectrum/ONT_BSseq/ONT_DLP_1stPre/L1_DNAme/LINE1_SVs.bed
# # bedMethyl009=/data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/bedMethyl_filtered/009T2_009T1_n_2_modBaseCalls_sorted_dedup_5mC_minCov15.bed
# # sort -k1,1 -k2,2n $bedMethyl009 > /data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/bedMethyl_filtered/009T2_009T1_n_2_modBaseCalls_sorted_dedup_5mC_minCov15_sorted.bed
# bedMethyl009=/data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/bedMethyl_filtered/009T2_009T1_n_2_modBaseCalls_sorted_dedup_5mC_minCov15_sorted.bed
# l1_minimal=/data1/greenbab/projects/methyl_benchmark_spectrum/ONT_BSseq/ONT_DLP_1stPre/L1_DNAme/LINE1_SVs_minimal.bed

# # sort -k1,1 -k2,2n $l1_minimal > l1_SVs_minimal.sorted.bed
# # bedtools map -a l1_SVs_minimal.sorted.bed -b ${bedMethyl009} -c 11 -o mean > 009T2_009T1_L1_methylation_new.bed


# # bedtools map -a ${l1_minimal} -b ${bedMethyl009} -c 11 -o mean > 009T2_009T1_L1_methylation_supossedlyMerged.bed

mkdir -p LINE1_promoters


l1Bed=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/SQuireL1MostActive_ValidnNonValidReadCounts_Aggregated_FWD_REV.bed

# l1Bed=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/SQuireRepeatsValidnNonValidReadCounts_Aggregated_FWD_REV.bed #full repeats
mm10Fa=/data1/greenbab/database/mm10/mm10.chrom.sizes

upstreamPromoter=300
nSlidingWindow=100
stepSize=5



#make ${upstreamPromoter} bp upstream of the 5'UTR
bedtools flank -i ${l1Bed} -g $mm10Fa -l 0 -r ${upstreamPromoter} -s > LINE1_promoters/l1_${upstreamPromoter}bp_5UTR_stranded.bed

# internal promoter region
bedtools flank -i "$l1Bed" -g $mm10Fa -l 600 -r 0 -s | \
awk '{start=$2+400; if(start<$2) start=$2; print $1, start, $3, $4, $5, $6}' OFS='\t' > LINE1_promoters/l1_400_600bp5UTR.bed

#make sliding windows of 10bp 
bedtools makewindows -b "$l1Bed" -w ${nSlidingWindow} -s ${stepSize} -i srcwinnum > LINE1_promoters/l1_${nSlidingWindow}bp_windows.bed


#test scripts
# head $l1Bed > head_l1.bed
# bedtools flank -i head_l1.bed -g $mm10Fa -l 0 -r ${upstreamPromoter} -s 
# bedtools flank -i head_l1.bed -g $mm10Fa -l 600 -r 0 -s | \
# awk '{start=$2+400; if(start<$2) start=$2; print $1, start, $3, $4, $5, $6}' OFS='\t'
# bedtools makewindows -b head_l1.bed -w 10 -s 10 -i srcwinnum

# genomeRegions=[300,1000]
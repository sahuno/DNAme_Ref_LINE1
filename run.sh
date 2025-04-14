
## run DNAme statistics for L1 promoters (L1Base)

#save results here;
# mkdir -p DNAme_L1Base_SummaryStats/nonCpGIslands
# mkdir -p DNAme_L1Base_SummaryStats/CpGIslandsOnly

# 400bp - 600bp 
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R \
--bedfiles "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAmeQC_Overlaps_L1Base/results/prepareBedFiles/nonCpGIslands" \
--L1ref "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.400_600bp5UTR_sorted.bed" \
--promoterType "400_600bpPromoter" \
--suffixBed "_minCov10.bed"

# --bedfiles "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/DNAme_bed/version2/results/prepareBedFiles/nonCpGIslands" \

# /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAmeQC_Overlaps_L1Base/results/BedGraphsBigWigs/nonCpGIslands/D-0-1_5000_4000

# whole length
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R \
--bedfiles "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAmeQC_Overlaps_L1Base/results/prepareBedFiles/nonCpGIslands" \
--L1ref "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.head.bed" \
--promoterType "wholeLengthPromoter" \
--suffixBed "_minCov10.bed"


## 100bpPromoter promoter
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R \
--bedfiles "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAmeQC_Overlaps_L1Base/results/prepareBedFiles/nonCpGIslands" \
--L1ref "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.100bp5UTR_sorted.bed" \
--promoterType "100bpPromoter" \
--suffixBed "_minCov10.bed"



#########################################
#### Run for cpg islands
#########################################

# 400bp - 600bp 
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R \
--bedfiles "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAmeQC_Overlaps_L1Base/results/prepareBedFiles/CpGIslandsOnly" \
--L1ref "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.400_600bp5UTR_sorted.bed" \
--promoterType "400_600bpPromoter" \
--suffixBed "_minCov10.bed"

# whole length
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R \
--bedfiles "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAmeQC_Overlaps_L1Base/results/prepareBedFiles/CpGIslandsOnly" \
--L1ref "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.head.bed" \
--promoterType "wholeLengthPromoter" \
--suffixBed "_minCov10.bed"


## 100bpPromoter promoter
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R \
--bedfiles "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAmeQC_Overlaps_L1Base/results/prepareBedFiles/CpGIslandsOnly" \
--L1ref "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.100bp5UTR_sorted.bed" \
--promoterType "100bpPromoter" \
--suffixBed "_minCov10.bed"

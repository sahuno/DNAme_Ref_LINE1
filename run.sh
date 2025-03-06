
## run statistics with promoters

# 400bp - 600bp 
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R \
--bedfiles "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/DNAme_bed/version2/results/prepareBedFiles/nonCpGIslands" \
--L1ref "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.400_600bp5UTR_sorted.bed" \
--promoterType "400_600bpPromoter" \
--suffixBed "_minCov10.bed"



# whole length
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R \
--bedfiles "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/DNAme_bed/version2/results/prepareBedFiles/nonCpGIslands" \
--L1ref "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.head.bed" \
--promoterType "wholeLengthPromoter" \
--suffixBed "_minCov10.bed"


## 100bpPromoter promoter
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R \
--bedfiles "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/DNAme_bed/version2/results/prepareBedFiles/nonCpGIslands" \
--L1ref "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.100bp5UTR_sorted.bed" \
--promoterType "100bpPromoter" \
--suffixBed "_minCov10.bed"
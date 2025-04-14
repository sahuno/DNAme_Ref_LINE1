
#run DNAme Summary Stats: Active Full length LINE1

~/miniforge3/envs/r-env/bin/Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R \
--bedfiles "../prepareDNAme/results/prepareBedFiles/nonCpGIslands" \
--L1ref "../database/L1BaseMmusculus/mmflil1_8438_noHeader.400_600bp5UTR_sorted.bed" \
--promoterType "400_600bpPromoter" \
--suffixBed "_minCov10.bed"
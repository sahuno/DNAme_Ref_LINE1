Author: samuel ahuno
date ; 11/06/2024

purpose: DNAme of reference line 1

Workflows:

```
snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/fullLengthL1_L1Base_DNAmeOverlaps.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/config/slurm --jobs unlimited --cores all --use-conda -np

```

```
sh /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/dDNAme_Ref_LINE1/run.sh
```

# run on slurm cluster

1. filter methyl bed files, sort, index & for CpG islands only

```
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/fullLengthL1_L1Base_DNAmeOverlaps.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs unlimited --cores all --keep-going --forceall -np
```


2. Genomic Overlaps in R;
```
# 400bp - 600bp 
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R \
--bedfiles "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/DNAme_bed/version2/results/prepareBedFiles/nonCpGIslands" \
--L1ref "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.400_600bp5UTR_sorted.bed" \
--promoterType "400_600bpPromoter" \
--suffixBed "_minCov10.bed"


Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R \
--bedfiles "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/DNAme_bed/version2/results/prepareBedFiles/nonCpGIslands" \
--L1ref "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.400_600bp5UTR_sorted.bed" \
--promoterType "400_600bpPromoter" \
--suffixBed "_minCov10.bed"

```

3. locus specific DE RNA-seq & overlaps with DNAme
```
scripts/R/locusSpecificRNA_DNAmeRate_with_fstTables.R
```

TODO: 
1. restrict methylation to 
a. internal promoter (antisense)
b. sense promoter

2. could you correct for read depth instead of remotving cpgs
3. DMR case vrs control
4. compare with RNA


Scripts
1. filter high quality cpgs 
```
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/fullLengthL1_L1Base_DNAmeOverlaps.smk
```

2. plot l1 DNAme distributions
```
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/GenomicOverlapsWithR.R
```

3. locus-specfic line 1 rna vrs dname
```
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/locusSpecificRNA_DNAmeRate.R
```




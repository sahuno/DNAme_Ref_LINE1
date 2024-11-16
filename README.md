author: samuel ahuno
date ; 11/06/2024

purpose: DNAme of reference line 1

# run on slurm cluster

```
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/fullLengthL1_L1Base_DNAmeOverlaps.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs unlimited --cores all --keep-going --forceall -np
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

3. locus -specfic line 1 rna vrs dname
```
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/locusSpecificRNA_DNAmeRate.R
```
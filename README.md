author:samuel ahuno
date ; 11/06/2024

purpose: dnam of reference line 1

#run on slurm cluster
```
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/fullLengthL1_L1Base_DNAmeOverlaps.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs unlimited --cores all --keep-going --forceall -np
```
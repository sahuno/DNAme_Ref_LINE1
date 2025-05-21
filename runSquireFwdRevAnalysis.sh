# run fwd and reverse Squie analysis
#1. generate count matrix for Repeats RNA
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/mergeSquire_RNA.R

#2. Generate DNAme for 
define promoter motifs for LINE1
./data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/Define_RepeatPromoters.sh

## generate DNAme stats
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/locSpecDNAme_overlaps.R

# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/DNAme_of_SquireRepeatsWithValidnNonValidReadCounts.R
# mkdir -p RepeatsDNAmeStats
snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/RepeatsPromotersStats.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/config/slurmMinimal --jobs unlimited --cores all --use-conda -np

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/RepeatsPromotersStats.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/config/slurmMinimal --jobs unlimited --cores all --use-conda -R merge_dna_rna_diagnostics -np


# rm -r .snakemake logs results


# /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/DNAme_of_SquireRepeatsWithValidnNonValidReadCounts.R


#3. merge repeats DNAme table
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/mergeRepeatsDNAmeTable.R


#4. Run locus-specific DE analysis
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/DE_locusSpecificLINE1.R --LocusSpecific TRUE

#5. HeatMaps for DE analysis
Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/DE_RepeatsHeatmaps.R \
--paths_rdata "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/LocusSpecific" --lfcCutoff 0.05


## Upset plot of DE analysis to know sharing
### Run for all repeats
# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/upsetPlots_DE.R \
#   --rna /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecific/resultsAllRepeatsRNA_DE_withMetadata.txt \
#   --dna /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/promoterMethyl/results/merged_repeatsDNAme_RNA_table_WithConds.tsv \
#   --focus CKi_DMSO,SETDB1i-CKi_DMSO,QSTAT-CKi_DMSO,QSTAT_DMSO,SETDB1i_DMSO \
#   --minsets 2 \
#   --prefix minsets2 \
#   --classes SINE,LTR,LINE,tRNA,scRNA,DNA,Other,RC,snRNA,RNA
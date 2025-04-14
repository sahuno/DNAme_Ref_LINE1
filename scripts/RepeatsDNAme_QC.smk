import glob as glob
import os

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/RepeatsDNAme_QC.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/config/slurmMinimal --jobs unlimited --cores all --use-conda -np


DIR="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03262025/prepareDNAme/results/prepareBedFiles/nonCpGIslands/"

files = glob.glob(os.path.join(DIR, "**", "*_sortedBed_minCov5.bed"), recursive=True)
samples = [os.path.basename(f).replace("_CpGs_sortedBed_minCov5.bed", "") for f in files]

rule all:
    input:
        expand("results/{sample}_DNAme_output.txt", sample=samples)

rule repeatsDNAme:
    input:
        lambda wildcards: os.path.join(DIR, f"{wildcards.sample}/{wildcards.sample}_CpGs_sortedBed_minCov5.bed")
    output:
        "results/{sample}_DNAme_output.txt"
    log: "logs/{sample}_log.txt"
    shell:
        """
        ~/miniforge3/envs/r-env/bin/Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/DNAme_of_SquireRepeatsWithValidnNonValidReadCounts.R \
        --RepeatRegions "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/SQuireRepeatsValidnNonValidReadCounts_Aggregated_FWD_REV.fst" \
        --FileDNAme "{input:q}" \
        --sampleName "{wildcards.sample:q}" 2> {log} && touch {output}
        """


# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/DNAme_of_SquireRepeatsWithValidReadCounts.R \
# --RepeatRegions "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/SQuireRepeatsValidReadCounts_Aggregated_FWD_REV.fst" \
# --FileDNAme "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03202025/DNAmeQC_Overlaps_L1Base/results/prepareBedFiles/nonCpGIslands/D-A-3_4000/D-A-3_4000_CpGs_sortedBed_minCov5.bed" \
# --sampleName "D-A-3_4000" && touch results/D-A-3_4000_DNAme_output.txt
#    conda:
#        "r-env"
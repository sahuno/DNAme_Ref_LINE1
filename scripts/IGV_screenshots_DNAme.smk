#bam IGV snapshots
# Usage:
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/IGV_screenshots_DNAme.smk --cores all --jobs unlimited --forcerun --printshellcmds --use-singularity --singularity-args "--bind /data1/greenbab" 

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/IGV_screenshots_DNAme.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/config/cluster_profiles/slurm --jobs 10 --cores all --use-singularity --singularity-args "-B /data1/greenbab" --keep-going --forceall -np

# rm -rf .snakemake
# rm -rf resultsIGV

import os
import glob as glob

fileSuffix="_CpGs_sortedBed_minCov10.bw"
# bam_dir = "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025/ALN/"
bw_dir = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/DNAme_bed/version2/results/BedGraphsBigWigs/nonCpGIslands"
sorted_bws = glob.glob(os.path.join(bw_dir, f"**/*{fileSuffix}"))

# Extract sample names (remove the .bam extension)
samples = [os.path.basename(bw).replace(fileSuffix, "") for bw in sorted_bws]
print(samples)

dirResults = "/results/"
reference_genome='/data1/greenbab/database/mm10/mm10.fa'

# os.path.join(bw_dir,'D-0-1_5000_4000' ,'D-0-1_5000_4000'+fileSuffix)
#test region file
# /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/regions_UID100.txt

rule all:
    input:
        "resultsIGV/done.txt"

rule igvSnapshot:
    input:
        bw = sorted_bws
    output:
        resultsDone = "resultsIGV/done.txt"
    params:
        REF = reference_genome,
        # regions = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/regions_UID100.txt"
        regions = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/regionsIgver/mmflil1_8438_noHeader.sorted.txt"
    threads: 8
    singularity: "/data1/greenbab/users/ahunos/apps/containers/igv_latest.sif"
    shell:
        """
        igver.py --bam {input.bw} -r {params.regions} \
            -o "resultsIGV/" \
                -mph 500 -od expand \
                    --genome 'mm10' && touch {output.resultsDone}
        """

# /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/regionsIgver/mmflil1_8438_noHeader.sorted.txt
        #-mph 500 -od expand \
#bam IGV snapshots
# Usage:
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/igv_batch.py --cores all --jobs unlimited --forcerun --printshellcmds --use-singularity --singularity-args "--bind /data1/greenbab" 

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/igv_batch.py --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/config/cluster_profiles/slurm --jobs 10 --cores all --use-singularity --singularity-args "-B /data1/greenbab" --keep-going --forceall -np

# rm -rf .snakemake
# rm -rf resultsIGV

import os
import glob as glob

bam_dir = "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025/ALN/"
# Collect unsorted BAM files that match the pattern (e.g. *.bam)
sorted_bams = glob.glob(os.path.join(bam_dir, "*_MM38.sorted.bam"))

# Extract sample names (remove the .bam extension)
samples = [os.path.basename(bam).replace("_MM38.sorted.bam", "") for bam in sorted_bams]
print(samples)

dirResults = "/results/"
reference_genome='/data1/greenbab/database/mm10/mm10.fa'

rule all:
    input:
        expand("resultsIGV/{sample}/{sample}.done.txt", sample=samples)

rule igvSnapshot:
    input:
        bam = os.path.join(bam_dir, "{sample}_MM38.sorted.bam")
    output:
        resultsDone = "resultsIGV/{sample}/{sample}.done.txt"
    params:
        REF = reference_genome,
        regions = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/regionsIgver/mmflil1_8438_noHeader.sorted.txt"
    threads: 8
    singularity: "/data1/greenbab/users/ahunos/apps/containers/igv_latest.sif"
    shell:
        """
        igver.py --bam {input.bam} -r {params.regions} \
            -o "resultsIGV/{wildcards.sample}/" \
                -mph 500 -od expand \
                    --genome 'mm10' && touch {output.resultsDone}
        """

# /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/regionsIgver/mmflil1_8438_noHeader.sorted.txt
        #-mph 500 -od expand \
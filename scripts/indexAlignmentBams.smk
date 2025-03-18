# Index BAM files
# Usage:
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/indexAlignmentBams.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/config/cluster_profiles/slurm --jobs 10 --cores all --keep-going --forceall -np

# rm -rf .snakemake

import os
import glob as glob

#bamList=glob.glob("/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025/ALN/*.bam")[1:3]
# baseNamebamList = [file.split("/")[:-1] for file in bamList]

#bam_files = glob.glob("/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025/ALN/*.bam")[1:3]
#sorted_bams = [bam.replace(".bam", ".sorted.bam") for bam in bam_files]

bam_dir = "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025/ALN/"
# Collect unsorted BAM files that match the pattern (e.g. *.bam)
unsorted_bams = glob.glob(os.path.join(bam_dir, "*.bam"))

# Extract sample names (remove the .bam extension)
samples = [os.path.basename(bam).replace(".bam", "") for bam in unsorted_bams]

rule all:
    input:
        # Final target: the index (.bai) files for each sorted BAM file
        expand(os.path.join(bam_dir, "{sample}.sorted.bam"), sample=samples),
        expand(os.path.join(bam_dir, "{sample}.sorted.bam.bai"), sample=samples)

rule sort_bam:
    input:
        bam = os.path.join(bam_dir, "{sample}.bam")
    output:
        sorted_bam = os.path.join(bam_dir, "{sample}.sorted.bam")
    threads: 8
    shell:
        """
        samtools sort -@ {threads} -o {output.sorted_bam} {input.bam}
        """

rule index_bam:
    input:
        sorted_bam = os.path.join(bam_dir, "{sample}.sorted.bam")
    output:
        bai = os.path.join(bam_dir, "{sample}.sorted.bam.bai")
    threads: 8
    shell:
        """
        samtools index -@ {threads} {input.sorted_bam}
        """
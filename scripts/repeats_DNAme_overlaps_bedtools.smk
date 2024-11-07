#input: 
#   1.bam 2.ref genome
# do this for later:
# goal: for any bedfile find the ovelaps
import os
import glob
from datetime import datetime
#purpose: phase snps and 5mc

# Assign the current date to a variable
current_datetime = datetime.now()
formatted_datetime = current_datetime.strftime('%Y%m%d_%H_%M_%S')

parent_dir = "/data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/"


configfile: parent_dir + "config/config.yaml"
# configfile: parent_dir + "config/samples_bams_5000sampleRate.yaml" #mouse samples; test run
configfile: parent_dir + "config/samples_bed.yaml" #mouse samples

#/data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/config/samples_bed_gz.yaml
# teref.mouse.fa
print(config["samples"]) #sanity check

#set species
set_species = "mouse"

#run to downlaod all necesay L1base files
# sh /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/scripts/download_L1Base_referencec.sh
#get basenames of files and append to output
L1Base_hg38=glob.glob('database/L1BaseHsapiens/*.sorted.bed', recursive=True)
L1Base_mm10=glob.glob('database/L1BaseMmusculus/*.sorted.bed', recursive=True)
# L1Base_hg38=['hsflil1_8438_noHeader.sorted.bed',  'hsflnil1_8438_rm_noHeader.sorted.bed',  'hsorf2l1_8438_noHeader.sorted.bed']

rule all:
    input:
        # expand('results/DNAme_overlaps/{samples}/done.{samples}.txt', samples=config["samples"]),
        expand('results/DNAme_overlaps/{samples}/{samples}_DNAme_RepEOverlaps.bed', samples=config["samples"])

rule DNAme_overlaps:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        RepeatsBed=lambda wildcards: L1Base_mm10 if set_species == "mouse" else L1Base_hg38,
        genomeFile=lambda wildcards: config["mm10_chrSizes"] if set_species == "mouse" else config["hg38_chrSizes"],
        threads=12,
        sortBedFiles=False
    # singularity: "/data1/greenbab/projects/yelena_long_read/te.sif"
    output:
        # done='results/DNAme_overlaps/{samples}/done.{samples}.txt',
        overlapsBed='results/DNAme_overlaps/{samples}/{samples}_DNAme_RepEOverlaps.bed'
        sortedBed='results/DNAme_overlaps/{samples}/{samples}_DNAme_RepEOverlaps.bed'
        # overlapsBed='results/DNAme_overlaps/{samples}/{samples}_DNAme_RepEOverlaps.bed'
        # overlapsBed='results/DNAme_overlaps/{samples}/{samples}_DNAme_RepEOverlaps.bed'
    log:
      "logs/DNAme_overlaps/{samples}/{samples}.log"
    run:
        if params.sortBedFiles:
            shell(
                  """
            sort -k1,1 -k2,2n {input} > {input.sortedBed}
            bedtools map -a {params.RepeatsBed[0]} -b {input} -c 11,11,11 -o mean,median,count -g {params.genomeFile} > {output.overlapsBed} 2> {log}
            """)
        else:
            shell("""
                bedtools map -a {params.RepeatsBed[0]} -b {input} -c 11,11,11 -o mean,median,count -g {params.genomeFile} > {output.overlapsBed} 2> {log}
            """)



    # bedtools map -a {params.RepeatsBed[1]} -b {input} -c 11,11,11 -o mean,median,count -g {params.genomeFile} > {output.overlapsBed} 2> {log}
    #     bedtools map -a {params.RepeatsBed[2]} -b {input} -c 11,11,11 -o mean,median,count -g {params.genomeFile} > {output.overlapsBed} 2> {log}


# rule DNAme_overlaps:
#     input:
#         lambda wildcards: config["samples"][wildcards.samples]
#     params:
#         reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
#         # RepeatsBed=lambda wildcards: config["RepE_mm10"] if set_species == "mouse" else config["consRef_human"],
#         genomeFile=lambda wildcards: config["mm10_chrSizes"] if set_species == "mouse" else config["hg38_chrSizes"],
#         threads=12
#     # singularity: "/data1/greenbab/projects/yelena_long_read/te.sif"
#     output:
#         done='results/DNAme_overlaps/{samples}/done.{samples}.txt',
#         overlapsBed='results/DNAme_overlaps/{samples}/{samples}_DNAme_RepEOverlaps.bed'
#     log:
#       "logs/DNAme_overlaps/{samples}/{samples}.log"
#     shell:
#         """
#         bedtools map -a {params.RepeatsBed} -b {input} -c 11,11,11 -o mean,median,count -g {params.genomeFile} > {output.overlapsBed} 2> {log}
#         """

#run interactively
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/scripts/repeats_DNAme_overlaps_bedtools.smk --cores 12 --forcerun --use-singularity --singularity-args "\"--bind /data1/greenbab\"" -np

#run on slurm
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/scripts/repeats_DNAme_overlaps_bedtools.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs 10 --cores all --use-singularity --singularity-args "'-B "/data1/greenbab"'" --keep-going --forceall -np

#methylation of full length line1- L1Base
# l1DB=/data1/greenbab/database/L1Base/L1base2_NCBIm38/mmflil1_8438_noHeader.sorted.bed
# D014000=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/modkit/D-0-1_4000/D-0-1_4000_modpileup_combined.bed.gz
# gfile=/data1/greenbab/database/mm10/chrom.sizes
# bedtools map -a ${l1DB} -b $D014000 -c 11,11,11 -o mean,median,count -g $gfile > D014000.mmflil1_DNAme.bed


# bedtools map -a /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed -b /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/modkit/D-0-1_5000_4000/D-0-1_5000_4000_modpileup_combined.bed.gz -c 11,11,11 -o mean,median,count -g /data1/greenbab/database/mm10/mm10.sorted.chrom.sizes > D-0-1_5000_4000_mmflil1_8438.bed


# sort -k1,1 -k2,2n /data1/greenbab/database/mm10/mm10.chrom.sizes >  /data1/greenbab/database/mm10/mm10.sorted.chrom.sizes

##get stand chrp=o
# awk '$0 !~ "_" {print $0}' /data1/greenbab/database/mm10/mm10.sorted.chrom.sizes | sort -k1,1 -k2,2n > /data1/greenbab/database/mm10/mm10.sorted.standard.chrom.sizes

#try with all chrom
# bedtools map -a /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed -b /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/modkit/D-0-1_5000_4000/D-0-1_5000_4000_modpileup_combined.bed.gz -c 11,11,11 -o mean,median,count -g  /data1/greenbab/database/mm10/mm10.sorted.chrom.sizes > D-0-1_5000_4000_mmflil1_8438.bed


##try with standard chrom files
# bedtools map -a /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed -b /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/modkit/D-0-1_5000_4000/D-0-1_5000_4000_modpileup_combined.bed.gz -c 11,11,11 -o mean,median,count -g  /data1/greenbab/database/mm10/mm10.sorted.standard.chrom.sizes > D-0-1_5000_4000_mmflil1_8438.bed



# bedtools map -a /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed -b /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/modkit/D-0-1_5000_4000/D-0-1_5000_4000_modpileup_combined.bed -c 11,11,11 -o mean,median,count -g /data1/greenbab/database/mm10/mm10.sorted.chrom.sizes > D-0-1_5000_4000_mmflil1_8438.bed
#cat /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed | cut -f 1 | sort | uniq | less
#cat /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed | cut -f 1 | sort | uniq | less
# cat /data1/greenbab/database/mm10/mm10.sorted.chrom.sizes | cut -f 1 | sort | uniq | less



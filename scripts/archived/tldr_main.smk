

#input: 
#   1.bam 2.ref genome

import os
import glob
from datetime import datetime
#purpose: phase snps and 5mc

# Assign the current date to a variable
current_datetime = datetime.now()
formatted_datetime = current_datetime.strftime('%Y%m%d_%H_%M_%S')

parent_dir = "/data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/"
#set species
set_species = "mouse"
configfile: parent_dir + "config/config.yaml"
# configfile: parent_dir + "config/samples_bams_5000sampleRate.yaml" #mouse samples; test run
configfile: parent_dir + "config/samples_bams_all.yaml" #mouse samples

# teref.mouse.fa
print(config["samples"]) #sanity check

#/data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/config/samples_bams_all.yaml

rule all:
    input:
        expand('results/tldr/{samples}/done.{samples}.txt', samples=config["samples"])

rule tldr:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        te_re=lambda wildcards: config["consRef_mm10"] if set_species == "mouse" else config["consRef_human"],
        software_dir=config["software_dir"],
        out_dir='results/call_snps_indels/{samples}/',
        threads=12,
    singularity: "/data1/greenbab/projects/yelena_long_read/te.sif"
    output:
        done='results/tldr/{samples}/done.{samples}.txt'
    log:
      "logs/tldr/{samples}/{samples}.log"
    shell:
        """
        tldr --bams {input} --outbase {wildcards.samples} --procs 10 --elts /tldr/ref/teref.human.fa --ref {params.reference_genome} --min_te_len 100 --max_cluster_size 100 --trdcol && touch {output.done} 2> {log}
        """


#run interactively
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/scripts/tldr_main.smk --cores 12 --forcerun --use-singularity --singularity-args "\"--bind /data1/greenbab\"" -np

#run on slurm
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/scripts/tldr_main.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs 10 --cores all --use-singularity --singularity-args "'-B "/data1/greenbab"'" --keep-going --forceall -np



# Snakefile
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/dmr_multi_groups.smk \
# --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3" \
# --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/config/slurmMinimal --jobs unlimited --cores all

import os
import os, glob

# 0) path to your ModKit Singularity image
MODKIT_SIF = "/data1/greenbab/users/ahunos/apps/containers/modkit_latest.sif"

# 1) load treatments from sample.yaml
configfile: "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/configs/treatments_dmr_groups.yaml"
SAMPLES = config["samples"]


# 0) find all of your region-BEDs and make a name→path map
REF_BEDS = glob.glob(
    "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/"
    "DNAme_Ref_LINE1/outputs/LINE1_promoters/first_*.bed"
)
REF_DICT = {
    os.path.splitext(os.path.basename(p))[0]: p
    for p in REF_BEDS
}
REF_NAMES = list(REF_DICT.keys())


# 2) fixed inputs & params
GENOME_CHROM_SIZES = "/data1/greenbab/database/mm10/mm10.sorted.chrom.sizes"
GENOME_FASTA       = "/data1/greenbab/database/mm10/mm10.fa"
LINE1_PROMOTERS    = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LINE1_promoters/first_100bp.bed"

# control_samples = {
#     "D-0-2_5000_4000": "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/results/Unphased_modkit_5mC_pileup/D-0-2_5000_4000/D-0-2_5000_4000_unphased_filtered_sorted.bed.gz",
#     "D-0-1_5000_4000": "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/results/Unphased_modkit_5mC_pileup/D-0-1_5000_4000/D-0-1_5000_4000_unphased_filtered_sorted.bed.gz",
#     "D-0-3_4000":      "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/results/Unphased_modkit_5mC_pileup/D-0-3_4000/D-0-3_4000_unphased_filtered_sorted.bed.gz"
# }

grouped_ctrl_samples =  {
    "DMSO": "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/results/Unphased_modkit_5mC_pileup/D-0-2_5000_4000/D-0-2_5000_4000_unphased_filtered_sorted.bed.gz",
    "DMSO": "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/results/Unphased_modkit_5mC_pileup/D-0-1_5000_4000/D-0-1_5000_4000_unphased_filtered_sorted.bed.gz",
    "DMSO":      "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/results/Unphased_modkit_5mC_pileup/D-0-3_4000/D-0-3_4000_unphased_filtered_sorted.bed.gz"
}


# build the control flags string
CONTROL_PARAMS = " ".join(f"-s {path} {name}"
                          for name, path in control_samples.items())

merged_CONTROL_PARAMS = " ".join(f"-s {path} {name}"
                          for name, path in grouped_ctrl_samples.items())
# list of treatment groups
TREATMENT_GROUPS = list(SAMPLES.keys())
rule all:
    input:
        expand("results/dmr/{group}/{ref}.done",
               group=TREATMENT_GROUPS,
               ref=REF_NAMES)

rule run_dmr:
    """
    Run modkit dmr multi for each treatment group AND each region‐BED.
    """
    output:
        outdir = directory("results/dmr/{group}/{ref}"),
        done   = touch("results/dmr/{group}/{ref}.done")
    params:
        control   = CONTROL_PARAMS,
        treat     = lambda wc: " ".join(
                        f"-s {path} {name}"
                        for name, path in SAMPLES[wc.group].items()
                     ),
        regionbed = lambda wc: REF_DICT[wc.ref]
    log:
        "logs/dmr_{group}_{ref}.log"
    threads: 8
    singularity:
        MODKIT_SIF
    shell:
        r"""
        modkit dmr multi \
          {params.control} \
          {params.treat} \
          --out-dir {output.outdir} \
          --ref {GENOME_FASTA} \
          --regions-bed {params.regionbed} \
          --base C \
          -t {threads} \
          -f \
          --header \
          --log-filepath {log}

        touch {output.done}
        """

#DO NOT RUN
# sort tabix any outputs in the outdir
        #   sort -k1,1 -k2,2n {output.segmentation_out} | bgzip -c > {output.segmentation_outGz}
        #   tabix -p bed {output.segmentation_outGz}


# rule run_merge_dmr:
#     """
#     Run modkit dmr multi for each treatment group AND each region‐BED.
#     """
#     output:
#         outdir = directory("results/dmr/{group}/{ref}"),
#         done   = touch("results/dmr/{group}/{ref}.done")
#     params:
#         control   = merged_CONTROL_PARAMS,
#         treat     = lambda wc: " ".join(
#                         f"-s {path} {name}"
#                         for name, path in SAMPLES[wc.group].items()
#                      ),
#         regionbed = lambda wc: REF_DICT[wc.ref]
#     log:
#         "logs/dmr_{group}_{ref}.log"
#     threads: 8
#     singularity:
#         MODKIT_SIF
#     shell:
#         r"""
#         modkit dmr multi \
#           {params.control} \
#           {params.treat} \
#           --out-dir {output.outdir} \
#           --ref {GENOME_FASTA} \
#           --regions-bed {params.regionbed} \
#           --base C \
#           -t {threads} \
#           -f \
#           --header \
#           --log-filepath {log}

#         touch {output.done}
#         """
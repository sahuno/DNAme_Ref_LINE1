# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/merged_dmr_multi.smk --cores 12 --forcerun -np


# Snakefile
import os
import glob

# 0) your ModKit Singularity image
MODKIT_SIF = "/data1/greenbab/users/ahunos/apps/containers/modkit_latest.sif"

# 1) load treatments from YAML (dict of dicts)
configfile: "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/configs/treatments_dmr_groups.yaml"
SAMPLES = config["samples"]

# 2) discover all of your region‐BEDs
REF_BEDS = glob.glob(
    "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/"
    "DNAme_Ref_LINE1/outputs/LINE1_promoters/first_*.bed"
)
REF_DICT = { os.path.splitext(os.path.basename(p))[0]: p for p in REF_BEDS }
REF_NAMES = list(REF_DICT.keys())

# 3) fixed inputs
GENOME_FASTA = "/data1/greenbab/database/mm10/mm10.fa"


# 4) controls, as a list of (label, path) so you can give them all "DMSO"
grouped_ctrl_samples = [
    ("DMSO", "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/results/Unphased_modkit_5mC_pileup/D-0-2_5000_4000/D-0-2_5000_4000_unphased_filtered_sorted.bed.gz"),
    ("DMSO", "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/results/Unphased_modkit_5mC_pileup/D-0-1_5000_4000/D-0-1_5000_4000_unphased_filtered_sorted.bed.gz"),
    ("DMSO", "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/results/Unphased_modkit_5mC_pileup/D-0-3_4000/D-0-3_4000_unphased_filtered_sorted.bed.gz")
]
MERGED_CONTROL_PARAMS = " ".join(
    f"-s {path} {label}" for label, path in grouped_ctrl_samples
)

# 5) our top‐level wildcards
TREATMENT_GROUPS = list(SAMPLES.keys())

rule all:
    input:
        expand("results/dmr/{group}/{ref}.done",
               group=TREATMENT_GROUPS,
               ref=REF_NAMES)

rule run_dmr:
    """
    Run modkit dmr multi for each (group, region) pair,
    merging all controls as DMSO and merging all treatment beds by group name.
    """
    output:
        outdir = directory("results/dmr/{group}/{ref}"),
        done   = touch("results/dmr/{group}/{ref}.done")
    params:
        control   = MERGED_CONTROL_PARAMS,
        # re-label every sample in this group as "{group}"
        treat     = lambda wc: " ".join(
                         f"-s {path} {wc.group}"
                         for path in SAMPLES[wc.group].values()
                     ),
        regionbed = lambda wc: REF_DICT[wc.ref]
    log:
        "logs/dmr_{group}_{ref}.log"
    threads: 8
    singularity:
        MODKIT_SIF
    shell:
        r"""
        echo "Running modkit dmr multi for {wildcards.group} and {wildcards.ref}"
        echo "Control samples: {params.control}"
        echo "Treatment samples: {params.treat}"
        echo "##########################################"

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

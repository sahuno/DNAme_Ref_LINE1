# Snakefile
# --use-singularity --singularity-args "\"--bind /data1/greenbab\"" 
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/get_promoter_DNAme.smk --cores 12 --forcerun -np
import os, glob
from snakemake.io import glob_wildcards, expand

# ─── CONFIG ────────────────────────────────────────────────────────────────────
GENOME_SIZES = "/data1/greenbab/database/mm10/mm10.sorted.chrom.sizes"

UNPHASED_DIR = (
    "/data1/greenbab/projects/triplicates_epigenetics_diyva/"
    "DNA/preprocessed/methylTables/"
    "results/Unphased_modkit_5mC_pileup"
)

# 1) Either hard-code your list of REF_BED paths:
#  REF_BEDS = [
#      "/data1/greenbab/users/ahunos/apps/"
#      "workflows/methylation_workflows/DNAme_Ref_LINE1/"
#      "outputs/repeatsBed/Repeats_DE_UpDown_LINE.bed",
#      "/path/to/another/refA.bed",
#      "/path/to/another/refB.bed",
#  ]

#  OR auto-discover them:
import glob as glob

REF_BEDS = glob.glob(
    "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LINE1_promoters/first_*.bed"
)

# Build a map from a short name → full path
REF_DICT = {
    os.path.splitext(os.path.basename(p))[0]: p
    for p in REF_BEDS
}
REF_NAMES = list(REF_DICT.keys())


# ─── SAMPLE DISCOVERY ──────────────────────────────────────────────────────────
SAMPLES = glob_wildcards(
    os.path.join(
        UNPHASED_DIR,
        "{sample}",
        "{sample}_unphased_filtered_sorted.bed"
    )
).sample


# ─── FINAL TARGETS ─────────────────────────────────────────────────────────────
rule all:
    input:
        expand(
            "results/{ref}/{sample}.bed",
            ref=REF_NAMES,
            sample=SAMPLES
        ),
        # plus the single merged file:
        "results/merged_all_samples.bed",
        # per-ref merged files
        expand("results/{ref}/merged_{ref}_samples.bed", ref=REF_NAMES),


# ─── MAPPING RULE ───────────────────────────────────────────────────────────────
rule map_repeats:
    """
    For each (ref,name) and sample run:
      bedtools map -a <refBed> -b <sample>.bed \
                   -c 11,11,11 -o mean,median,count \
                   -g <GENOME_SIZES> > results/{ref}/{sample}.bed
    """
    input:
        refBed      = lambda wc: REF_DICT[wc.ref],
        bedSample   = lambda wc: os.path.join(
            UNPHASED_DIR,
            wc.sample,
            f"{wc.sample}_unphased_filtered_sorted.bed"
        ),
        genomeSizes = GENOME_SIZES
    output:
        "results/{ref}/{sample}.bed"
    shell:
        """
        mkdir -p $(dirname {output})
        /data1/greenbab/users/ahunos/apps/bedtools2/bin/bedtools map \
            -a {input.refBed} \
            -b {input.bedSample} \
            -c 11,11,11 \
            -o mean,median,count \
            -g {input.genomeSizes} \
        > {output}
        """

# ─── PER-REF MERGE ──────────────────────────────────────────────────────────────
rule merge_per_ref:
    input:
        lambda wc: expand(
            "results/{ref}/{sample}.bed",
            ref=wc.ref,
            sample=SAMPLES
        )
    output:
        "results/{ref}/merged_{ref}_samples.bed"
    shell:
        """
        mkdir -p $(dirname {output})
        > {output}
        for f in {input}; do
            sample=$(basename $f .bed)
            # prefix each line with sample name
            awk -v s="$sample" '{{ print s"\t"$0 }}' "$f"
        done >> {output}
        """

# ─── GLOBAL MERGE ───────────────────────────────────────────────────────────────
rule merge_all:
    input:
        expand("results/{ref}/{sample}.bed", ref=REF_NAMES, sample=SAMPLES)
    output:
        "results/merged_all_samples.bed"
    shell:
        """
        > {output}
        for f in {input}; do
            sample=$(basename $f .bed)
            awk -v s="$sample" '{{ print s"\t"$0 }}' "$f"
        done >> {output}
        """
import glob as glob
import os
#how to run
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/RepeatsPromotersStats.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/config/slurmMinimal --jobs unlimited --cores all --use-conda -np


DIR="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/clean_analysis03262025/prepareDNAme/results/prepareBedFiles/nonCpGIslands/"

files = glob.glob(os.path.join(DIR, "**", "*_sortedBed_minCov5.bed"), recursive=True)
samples = [os.path.basename(f).replace("_CpGs_sortedBed_minCov5.bed", "") for f in files]
RULENAME = "repeatsDNAme"

rule all:
    input:
        # Repeat DNAme summary table
        expand("results/{rulename}/{sample}/{sample}_repeatsDNAme.tsv", sample=samples, rulename=RULENAME),
        # Plots
        expand("results/{rulename}/{sample}/{sample}_repeats_methylation.pdf", sample=samples, rulename=RULENAME),
        expand("results/{rulename}/{sample}/{sample}_L1_methylation.pdf", sample=samples, rulename=RULENAME),
        expand("results/{rulename}/{sample}/{sample}_L1MD_methylation_pval.pdf", sample=samples, rulename=RULENAME),
        expand("results/{rulename}/{sample}/{sample}_DNA_MethylStats_ggpairs.pdf", sample=samples, rulename=RULENAME),
        expand("results/{rulename}/{sample}/{sample}_DNA_MethylStats_ggpairs.png", sample=samples, rulename=RULENAME),
        expand("results/{rulename}/{sample}/{sample}_DNA_MethylStats_ggpairsCombineValidNonValidReadCounts.png", sample=samples, rulename=RULENAME),
        expand("results/{rulename}/{sample}/{sample}_DNA_MethylStats_ggpairsConsensus.png", sample=samples, rulename=RULENAME),
        # The done file to mark completion
        expand("results/{rulename}/{sample}/{sample}_DNAme_output.txt", sample=samples, rulename=RULENAME),
        "results/mergeRNADNAme/merged_repeatsDNAme_RNA_table.tsv",
        "results/merge_dna_rna_diagnostics/merge_dna_rna_diagnostics_done.txt"
#        "results/merge_dna_rna_diagnostics/boxplot_fracMethyl_geom_mean_by_sample.png",
#        "results/merge_dna_rna_diagnostics/boxplot_fracMethyl_geom_mean_by_consensus.png",
#        "results/merge_dna_rna_diagnostics/boxplot_fracMethyl_geom_mean_by_consensusCond.png",
#        "results/merge_dna_rna_diagnostics/scatter_fpkm_vs_fracMethyl_geom_mean.png"
#        "results/merge_dna_rna_diagnostics/merged_repeatsDNAme_RNA_table_WithConds.tsv",

rule repeatsDNAme:
    input:
        lambda wildcards: os.path.join(DIR, f"{wildcards.sample}/{wildcards.sample}_CpGs_sortedBed_minCov5.bed")
    output:
        doneTxt="results/{rulename}/{sample}/{sample}_DNAme_output.txt",
        repeat_summary="results/{rulename}/{sample}/{sample}_repeatsDNAme.tsv",
        repeat_plot="results/{rulename}/{sample}/{sample}_repeats_methylation.pdf",
        L1_plot="results/{rulename}/{sample}/{sample}_L1_methylation.pdf",
        L1MD_plot="results/{rulename}/{sample}/{sample}_L1MD_methylation_pval.pdf",
        ggpairs_pdf="results/{rulename}/{sample}/{sample}_DNA_MethylStats_ggpairs.pdf",
        ggpairs_png="results/{rulename}/{sample}/{sample}_DNA_MethylStats_ggpairs.png",
        ggpairs_combine="results/{rulename}/{sample}/{sample}_DNA_MethylStats_ggpairsCombineValidNonValidReadCounts.png",
        ggpairs_consensus="results/{rulename}/{sample}/{sample}_DNA_MethylStats_ggpairsConsensus.png"
    log: "logs/{rulename}/{sample}/{sample}_log.txt"
    params:
        bedRepeat="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LINE1_promoters/l1_300bp_5UTR_stranded.bed",
        scriptR="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/R/locSpecDNAme_overlaps.R",
        outdir="results/{rulename}/{sample}",
        rnaSeqCountsTable="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/squireRepeats_AggregatedFWDnREV_masterTable.fst"
    shell:
        """
        mkdir -p {params.outdir}
        ~/miniforge3/envs/r-env/bin/Rscript {params.scriptR} \
        --RepeatRegions {params.bedRepeat} \
        --FileDNAme "{input:q}" \
        --sampleName "{wildcards.sample:q}" \
        --outDir {params.outdir} \
        --RNASeqCountsPath {params.rnaSeqCountsTable} \
        --runPlots TRUE 2> {log} && touch {output.doneTxt}
        """

rule mergeRNADNAme:
    input:
        expand("results/{rulename}/{sample}/{sample}_repeatsDNAme.tsv", sample=samples, rulename=RULENAME)
    output:
        "results/mergeRNADNAme/merged_repeatsDNAme_RNA_table.tsv"
    log: "logs/mergeRNADNAme/merge_repeatsDNAme_RNA_table.log"
    shell:
        """
        mkdir -p results/mergeRNADNAme
        awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input} > {output} 2> {log}
        """

rule merge_dna_rna_diagnostics:
    input:
        merged_path="results/mergeRNADNAme/merged_repeatsDNAme_RNA_table.tsv",
        metadata="/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_recoded.csv"
    output:
        dir=directory("results/merge_dna_rna_diagnostics"),
        doneTxt="results/merge_dna_rna_diagnostics/merge_dna_rna_diagnostics_done.txt"
    params:
        script="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/mergeDNA_RNA_diagnostics.R"
    log:
        "logs/merge_dna_rna_diagnostics/merge_dna_rna_diagnostics.log"
    shell:
        """
        ~/miniforge3/envs/r-env/bin/Rscript {params.script} \
            --mergedPath {input.merged_path} \
            --metadata {input.metadata} \
            --outDir {output.dir} && touch {output.doneTxt}
        """
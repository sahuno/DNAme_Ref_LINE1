parent_dir = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/"
configfile: parent_dir + "config/config.yaml"
configfile: parent_dir + "config/samples_bed.yaml" #mouse samples

minCoverages = [5, 10, 15, 20]


#to Run: 
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/fullLengthL1_L1Base_DNAmeOverlaps.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs unlimited --cores all --use-conda -np

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/fullLengthL1_L1Base_DNAmeOverlaps.smk --jobs 10 --cores all --keep-going --forceall -np

#decision
#1: what about other stats
# bedtools map -o collapse #to get all values which could then be used for other stats

#to rerun without interfeerance
# rm -rf .snakemake logs results 


###TODO: 
#how about promoter line1?

import glob
def getSpecificSamplesPaths(filterKeyword):
    return {v for k,v in config["samples"].items() if filterKeyword in k}


#glob recursiveltly and filter files that have that pattern
# def get_Samples_CpGIslandsBed(minCov):
#     pattern = f'results/DNAme_overlaps/*/overlaps/*_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}_CpGIslands.bed'
#     return glob.glob(pattern)

# def get_Samples(filterKeyword):
#     pattern = f'results/DNAme_overlaps/*/overlaps/*minCov{filterKeyword}*.bed'
#     listFiles=glob.glob(pattern, recursive=True)
#     return listFiles

def get_Overlaps_minCov(filterKeyword, filterKeyword2):
    pattern = f'results/DNAme_overlaps/*/overlaps/*minCov{filterKeyword}*{filterKeyword2}.bed'
    listFiles=glob.glob(pattern, recursive=True)
    return listFiles

# get_Overlaps_minCov(filterKeyword=15, filterKeyword2="CpGIslands")


rule all:
    input:
        expand('results/DNAme_overlaps/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed', samples=config["samples"],minCov = minCoverages),
        expand('results/DNAme_overlaps/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed', samples=config["samples"], minCov = minCoverages),
        expand('results/DNAme_overlaps/{samples}/overlaps/{samples}_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}.bed', samples=config["samples"],minCov=minCoverages),
        expand('results/DNAme_overlaps/{samples}/overlaps/{samples}_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}_CpGIslands.bed', samples=config["samples"],minCov=minCoverages),
        #use for dmr
        expand('results/DNAme_overlaps/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed.gz', samples=config["samples"],minCov = minCoverages),
        expand('results/DNAme_overlaps/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed.gz.tbi', samples=config["samples"],minCov = minCoverages),
        expand('results/DNAme_overlaps/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed.gz', samples=config["samples"],minCov = minCoverages),
        expand('results/DNAme_overlaps/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed.gz.tbi', samples=config["samples"], minCov = minCoverages),
        expand("results/figures/CGI/done.{minCov}.txt", minCov = minCoverages),
        expand("results/figures/nonCGI/done.{minCov}.txt", minCov = minCoverages)

rule DNAme_overlaps:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    params:
        RepeatsBed="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed",
        genomeFile="/data1/greenbab/database/mm10/mm10.sorted.chrom.sizes",
        cpgIsland="/data1/greenbab/database/mm10/mm10_CpGIslands.bed",
        threads=12,
        minCov=minCoverages,
        sortBedFiles=True
    output:
        sortedBed='results/DNAme_overlaps/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed',
        sortedBedCpGIslandsOnly='results/DNAme_overlaps/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed',
        sortedBedGz='results/DNAme_overlaps/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed.gz',
        sortedBedGzTbi='results/DNAme_overlaps/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed.gz.tbi',
        sortedBedCpGIslandsOnlyGz='results/DNAme_overlaps/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed.gz',
        sortedBedCpGIslandsOnlyGzTbi='results/DNAme_overlaps/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed.gz.tbi',

        L1overlapsCpGIslandsOnly='results/DNAme_overlaps/{samples}/overlaps/{samples}_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}_CpGIslands.bed',
        L1overlapsDNAme='results/DNAme_overlaps/{samples}/overlaps/{samples}_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}.bed'
    log:
      "logs/DNAme_overlaps/{samples}/{samples}_minCov{minCov}.log"
    run:
        if params.sortBedFiles:
            shell("""
                awk -v min_cov="{params.minCov}" 'BEGIN { OFS = "\t" } ($10 > min_cov) {{$11=$11/100; print}}' "{input}" | sort -k1,1 -k2,2n > {output.sortedBed}
                bgzip -k {output.sortedBed} && tabix -p bed {output.sortedBedGz}
                
                bedtools intersect -a {output.sortedBed} -b {params.cpgIsland} > {output.sortedBedCpGIslandsOnly}
                bgzip -k {output.sortedBedCpGIslandsOnly} && tabix -p bed {output.sortedBedCpGIslandsOnlyGz}

                bedtools map -a {params.RepeatsBed} -b {output.sortedBedCpGIslandsOnly} -c 11,11,11,11,11 -o min,max,mean,median,count -g {params.genomeFile} > {output.L1overlapsCpGIslandsOnly} 2> {log}
                bedtools map -a {params.RepeatsBed} -b {output.sortedBed} -c 11,11,11,11,11 -o min,max,mean,median,count -g {params.genomeFile} > {output.L1overlapsDNAme} 2> {log}
            """)
        else:
            shell("""
                bedtools map -a {params.RepeatsBed} -b {input} -c 11,11,11,11,11 -o min,max,mean,median,count -g {params.genomeFile} > {output.overlapsBed} 2> {log}
            """)

#maybe glob.glob the directly
rule plotRegions:
    input:
        inL1overlapsCpGIslandsOnly=lambda wildcards: get_Overlaps_minCov(filterKeyword=wildcards.minCov, filterKeyword2="CpGIslands"),
        #'results/DNAme_overlaps/{samples}/overlaps/{samples}_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}_CpGIslands.bed',
        inL1overlapsDNAme=lambda wildcards: get_Overlaps_minCov(filterKeyword=wildcards.minCov, filterKeyword2="")
    output:
        CgiPlots="results/figures/CGI/done.{minCov}.txt",
        nonCgiPlots="results/figures/nonCGI/done.{minCov}.txt"
    conda: "r-env"
    log:
      logCGI="logs/figures/nonCGI/plotLog_CGI_minCov.{minCov}.txt",
      lognoCGI="logs/figures/nonCGI/plotLog_noCGI_minCov.{minCov}.txt"

    shell:
        """ 
        Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/plot_within_workflow.R --files "{input.inL1overlapsCpGIslandsOnly:q}" 2> {log.logCGI} && touch {output.CgiPlots}
        Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/plot_within_workflow.R --files "{input.inL1overlapsDNAme:q}" 2> {log.lognoCGI} && touch {output.nonCgiPlots}
        """

# Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/plot_within_workflow.R --files {input:q} --metadata 


# /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/results/DNAme_overlaps/D-0-1_5000_4000/D-0-1_5000_4000_5mCpG_5hmCpG_sortedBed_minCov20_CpGIslands.bed.gz ctrl1
# /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/results/DNAme_overlaps/D-S-1_5000_4000/D-S-1_5000_4000_5mCpG_5hmCpG_sortedBed_minCov20_CpGIslands.bed.gz setdb11


# bedMinCov=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/results/DNAme_overlaps/D-0-1_5000_4000/D-0-1_5000_4000_5mCpG_5hmCpG_sortedBed_minCov10.bed

# bedMinCovCGIs=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/results/DNAme_overlaps/D-0-1_5000_4000/D-0-1_5000_4000_5mCpG_5hmCpG_sortedBed_minCov10_CpGIslands.bed

#run independelbt 
# bgzip -k ${bedMinCov} -o D-0-1_5000_4000_5mCpG_5hmCpG_sortedBed_minCov10.bed.gz
# tabix -p bed D-0-1_5000_4000_5mCpG_5hmCpG_sortedBed_minCov10.bed.gz

# # pipe or combien commands
# bgzip -k ${bedMinCov} && \
# tabix -p bed ${bedMinCov}.gz




# rule DMR:
#     input:
#         controls=lambda wildcards: getSpecificSamplesPaths('D-0-'),
#         cases=lambda wildcards: getSpecificSamplesPaths('D-S-')
#         # controls=lambda wildcards: getSpecificSamples('D-0-')[wildcards.samples],
#         # cases=lambda wildcards: getSpecificSamples('D-S-')[wildcards.samples]
#     output:
#         doneFile='results/casesControlsJobs/doneFile.txt'
#     shell:
#         """ 
#         echo Controls: {input.controls} >> {output.doneFile}
#         echo "\n cases paths\n" >> {output.doneFile}
#         echo Cases: {input.cases} >> {output.doneFile}
#         """

# inputBedfile=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/modkit/D-0-1_5000_4000/D-0-1_5000_4000_modpileup_combined.bed
# cpgIsland=/data1/greenbab/database/mm10/mm10_CpGIslands.bed
# bedtools intersect -a ${inputBedfile} -b ${cpgIsland} > D-0-1_5000_4000_cggisland.bed

# less /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/D-0-1_5000_4000_cggisland.bed

# # map cpg methylation values
# bedtools map -a /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed \
# -b /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/D-0-1_5000_4000_cggisland.bed \
# -c 11,11,11,11,11 -o min,max,mean,median,count \
# -g /data1/greenbab/database/mm10/mm10.sorted.chrom.sizes > D-0-1_5000_4000_RepeatsDNAmeIncggisland.bed

                # sort -k1,1 -k2,2n {input} > {output.sortedBed}
## Add rule to plot all overlaps
## add command to do cpg only
# find -type f -name "*.polars"


#plot with python
# import polars as pl
# paths_ls = ["/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/results/DNAme_overlaps/D-0-1_5000_4000/overlaps/D-0-1_5000_4000_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov15.bed",
# "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/sandbox/results/DNAme_overlaps/D-S-1_5000_4000/overlaps/D-S-1_5000_4000_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov15.bed"]

# pls_df = pl.scan_csv(paths_ls)
# plsHead = pls_df.head(10)#.collect()
# plsHead.collect()
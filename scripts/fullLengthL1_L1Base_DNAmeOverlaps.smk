parent_dir = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/"
configfile: parent_dir + "config/config.yaml"
configfile: parent_dir + "config/samples_bed.yaml" #mouse samples

minCoverages = [5, 10, 15, 20]


#to Run: 
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/fullLengthL1_L1Base_DNAmeOverlaps.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/config/slurm --jobs unlimited --cores all --use-conda -np

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/fullLengthL1_L1Base_DNAmeOverlaps.smk --jobs 10 --cores all --keep-going --forceall -np

#decision
#1: what about other stats
# bedtools map -o collapse #to get all values which could then be used for other stats

#to rerun without interfeerance
# rm -rf .snakemake logs results *.png

#l .snakemake/slurm_logs/rule_DNAme_overlaps/D-Q-2_4000_5/

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

# def get_Overlaps_minCov(filterKeyword, filterKeyword2):
#     pattern = f'results/DNAme_overlaps/*/overlaps/*minCov{filterKeyword}*{filterKeyword2}.bed'
#     listFiles=glob.glob(pattern, recursive=True)
#     return listFiles

##takes paths and mean coverage values 
# def get_Overlaps_minCov(paths, minCoverage):
#     pattern = f"{paths}/*_minCov{minCoverage}*.bed"
#     listFiles=glob.glob(pattern, recursive=True)
#     if not listFiles:
#         print(f"Warning: No files found for pattern {pattern}")
#     return listFiles

##takes paths and mean coverage values 
def get_Overlaps_minCov(base_path_pattern, minCoverage):
    # Construct the pattern with the wildcard
    pattern = f"{base_path_pattern}/*_minCov{minCoverage}*.bed"
    listFiles = glob.glob(pattern, recursive=True)
    if not listFiles:
        print(f"Warning: No files found for pattern {pattern}")
    return listFiles

# get_Overlaps_minCov(filterKeyword=15, filterKeyword2="CpGIslands")


rule all:
    input:
        expand('results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed', samples=config["samples"],minCov = minCoverages),
        expand('results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed', samples=config["samples"], minCov = minCoverages),
        #use for dmr
        expand('results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed.gz', samples=config["samples"],minCov = minCoverages),
        expand('results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed.gz.tbi', samples=config["samples"],minCov = minCoverages),
        expand('results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed.gz', samples=config["samples"],minCov = minCoverages),
        expand('results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed.gz.tbi', samples=config["samples"], minCov = minCoverages),
        #CpGs across the full lengths
        expand('results/OverlapsL1PromoterDNAme/{samples}/fullLength/CpGIs/{samples}_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}_CpGIslands.bed', samples=config["samples"],minCov=minCoverages),
        expand('results/OverlapsL1PromoterDNAme/{samples}/fullLength/nonCGI/{samples}_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}.bed', samples=config["samples"],minCov=minCoverages),
        #L1 promoter upstream 
        expand('results/OverlapsL1PromoterDNAme/{samples}/100bp5UTR/CpGIs/{samples}_100bp5UTR_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}_CpGIslands.bed', samples=config["samples"], minCov = minCoverages),
        expand('results/OverlapsL1PromoterDNAme/{samples}/100bp5UTR/nonCGI/{samples}_100bp5UTR_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}.bed', samples=config["samples"], minCov = minCoverages),
        expand('results/OverlapsL1PromoterDNAme/{samples}/400600bp5UTR/CpGIs/{samples}_400600bp5UTR_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}_CpGIslands.bed', samples=config["samples"], minCov = minCoverages),
        expand('results/OverlapsL1PromoterDNAme/{samples}/400600bp5UTR/nonCGI/{samples}_400600bp5UTR_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}.bed', samples=config["samples"], minCov = minCoverages),
        ##plot results 
        expand("results/figures/fullLength/CpGIs/done.{minCov}.txt", minCov=minCoverages),
        expand("results/figures/fullLength/nonCGI/done.{minCov}.txt", minCov=minCoverages),
        expand("results/figures/100bp5UTR/CpGIs/done.{minCov}.txt", minCov=minCoverages),
        expand("results/figures/100bp5UTR/nonCGI/done.{minCov}.txt", minCov=minCoverages),
        expand("results/figures/400600bp5UTR/CpGIs/done.{minCov}.txt", minCov=minCoverages),
        expand("results/figures/400600bp5UTR/nonCGI/done.{minCov}.txt", minCov=minCoverages)

        # expand("results/figures/fullLength/CpGIs/done.{minCov}.txt", minCov = minCoverages),
        # expand("results/figures/fullLength/nonCGI/done.{minCov}.txt", minCov = minCoverages)

rule prepareBedFiles:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    params:
        RepeatsFulllengthBed="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed",
        promoter100bp="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.100bp5UTR_sorted.bed",
        promoter400_600bp="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.400_600bp5UTR_sorted.bed",
        genomeFile="/data1/greenbab/database/mm10/mm10.sorted.chrom.sizes",
        cpgIsland="/data1/greenbab/database/mm10/mm10_CpGIslands.bed",
        threads=12,
        minCov=minCoverages,
        sortBedFiles=True
    output:
        sortedBed='results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed',
        sortedBedCpGIslandsOnly='results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed',
        sortedBedGz='results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed.gz',
        sortedBedGzTbi='results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed.gz.tbi',
        sortedBedCpGIslandsOnlyGz='results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed.gz',
        sortedBedCpGIslandsOnlyGzTbi='results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed.gz.tbi',

    log:
      "logs/prepareBedFiles/{samples}/{samples}_minCov{minCov}.log"
    run:
        if params.sortBedFiles:
            shell("""
                awk -v min_cov='{wildcards.minCov}' 'BEGIN {{ OFS = "\\t" }} ($10 > min_cov) {{$11=$11/100; print}}' "{input}" | sort -k1,1 -k2,2n > {output.sortedBed}
                bgzip -k {output.sortedBed} && tabix -p bed {output.sortedBedGz}
                
                bedtools intersect -a {output.sortedBed} -b {params.cpgIsland} > {output.sortedBedCpGIslandsOnly}
                bgzip -k {output.sortedBedCpGIslandsOnly} && tabix -p bed {output.sortedBedCpGIslandsOnlyGz}
            """)
        else:
            shell("""
                bedtools map -a {params.RepeatsBed} -b {input} -c 11,11,11,11,11 -o min,max,mean,median,count -g {params.genomeFile} > {output.overlapsBed} 2> {log}
            """)

                # bedtools map -a {params.promoter100bp} -b {output.sortedBedCpGIslandsOnly} -c 11,11,11,11,11 -o min,max,mean,median,count -g {params.genomeFile} > {output.L1_100bp5UTR_overlapsCpGIslandsOnly} 2> {log}
                # bedtools map -a {params.promoter400_600bp} -b {output.sortedBed} -c 11,11,11,11,11 -o min,max,mean,median,count -g {params.genomeFile} > {output.L1_400_600bp5UTR_overlapsDNAme} 2> {log}

        # L1_100bp5UTR_overlapsCpGIslandsOnly='results/DNAme_overlaps/{samples}/overlaps/100bp5UTR/{samples}_100bp5UTR_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}_CpGIslands.bed',
        # L1_400_600bp5UTR_overlapsDNAme='results/DNAme_overlaps/{samples}/overlaps/400600bp5UTR/{samples}_400600bp5UTR_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}.bed'


rule OverlapsL1PromoterDNAme:
    input:
        sortedBed='results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed',
        sortedBedCpGIslandsOnly='results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed',
        sortedBedGz='results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed.gz',
        sortedBedGzTbi='results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}.bed.gz.tbi',
        sortedBedCpGIslandsOnlyGz='results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed.gz',
        sortedBedCpGIslandsOnlyGzTbi='results/prepareBedFiles/{samples}/{samples}_5mCpG_5hmCpG_sortedBed_minCov{minCov}_CpGIslands.bed.gz.tbi',

    output:
        L1overlapsCpGIslandsOnly='results/OverlapsL1PromoterDNAme/{samples}/fullLength/CpGIs/{samples}_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}_CpGIslands.bed',
        L1overlapsDNAme='results/OverlapsL1PromoterDNAme/{samples}/fullLength/nonCGI/{samples}_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}.bed',
        ##overlaps with L1 promoters; 100bp upstream & 400-600bp upstream 
        L1_100bp5UTR_overlapsCpGIslandsOnly='results/OverlapsL1PromoterDNAme/{samples}/100bp5UTR/CpGIs/{samples}_100bp5UTR_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}_CpGIslands.bed',
        L1_100bp5UTR_overlapsnCpGIslandsOnly='results/OverlapsL1PromoterDNAme/{samples}/100bp5UTR/nonCGI/{samples}_100bp5UTR_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}.bed',
        L1_400_600bp5UTR_overlapsCGIs_DNAme='results/OverlapsL1PromoterDNAme/{samples}/400600bp5UTR/CpGIs/{samples}_400600bp5UTR_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}_CpGIslands.bed',
        L1_400_600bp5UTR_overlapsnCGIs_DNAme='results/OverlapsL1PromoterDNAme/{samples}/400600bp5UTR/nonCGI/{samples}_400600bp5UTR_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}.bed'
    params:
        RepeatsFulllengthBed="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed",
        promoter100bp="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.100bp5UTR_sorted.bed",
        promoter400_600bp="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.400_600bp5UTR_sorted.bed",
        genomeFile="/data1/greenbab/database/mm10/mm10.sorted.chrom.sizes",
        cpgIsland="/data1/greenbab/database/mm10/mm10_CpGIslands.bed",
        threads=12,
        minCov=minCoverages,
        sortBedFiles=True
    shell:
        """
        bedtools map -a {params.RepeatsFulllengthBed} -b {input.sortedBedCpGIslandsOnly} -c 11,11,11,11,11 -o min,max,mean,median,count -g {params.genomeFile} > {output.L1overlapsCpGIslandsOnly}
        bedtools map -a {params.RepeatsFulllengthBed} -b {input.sortedBed} -c 11,11,11,11,11 -o min,max,mean,median,count -g {params.genomeFile} > {output.L1overlapsDNAme}
        bedtools map -a {params.promoter100bp} -b {input.sortedBedCpGIslandsOnly} -c 11,11,11,11,11 -o min,max,mean,median,count -g {params.genomeFile} > {output.L1_100bp5UTR_overlapsCpGIslandsOnly}
        bedtools map -a {params.promoter100bp} -b {input.sortedBed} -c 11,11,11,11,11 -o min,max,mean,median,count -g {params.genomeFile} > {output.L1_100bp5UTR_overlapsnCpGIslandsOnly}
        bedtools map -a {params.promoter400_600bp} -b {input.sortedBedCpGIslandsOnly} -c 11,11,11,11,11 -o min,max,mean,median,count -g {params.genomeFile} > {output.L1_400_600bp5UTR_overlapsCGIs_DNAme}
        bedtools map -a {params.promoter400_600bp} -b {input.sortedBed} -c 11,11,11,11,11 -o min,max,mean,median,count -g {params.genomeFile} > {output.L1_400_600bp5UTR_overlapsnCGIs_DNAme}
        """






#maybe glob.glob the directly
rule plotRegions:
    input:
        L1fullLengthsOverlapsCGIs=lambda wildcards: get_Overlaps_minCov(base_path_pattern="results/OverlapsL1PromoterDNAme/*/fullLength/CpGIs", minCoverage=wildcards.minCov),
        L1fullLengthsOverlapsnCGIs=lambda wildcards: get_Overlaps_minCov(base_path_pattern="results/OverlapsL1PromoterDNAme/*/fullLength/nonCGI", minCoverage=wildcards.minCov),
        #'results/DNAme_overlaps/{samples}/overlaps/{samples}_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}_CpGIslands.bed',
        # inL1overlapsDNAme=lambda wildcards: get_Overlaps_minCov(filterKeyword=wildcards.minCov, filterKeyword2="")
        l1_100bp5UTRCGIs=lambda wildcards: get_Overlaps_minCov(base_path_pattern="results/OverlapsL1PromoterDNAme/*/100bp5UTR/CpGIs", minCoverage=wildcards.minCov),
        l1_100bp5UTRnCGIs=lambda wildcards: get_Overlaps_minCov(base_path_pattern="results/OverlapsL1PromoterDNAme/*/100bp5UTR/nonCGI", minCoverage=wildcards.minCov),
        L1_400_600bp5UTRCGIs=lambda wildcards: get_Overlaps_minCov(base_path_pattern="results/OverlapsL1PromoterDNAme/*/400600bp5UTR/CpGIs", minCoverage=wildcards.minCov),
        L1_400_600bp5UTRnCGIs=lambda wildcards: get_Overlaps_minCov(base_path_pattern="results/OverlapsL1PromoterDNAme/*/400600bp5UTR/nonCGI", minCoverage=wildcards.minCov)

    output:
        outFLCgiPlots="results/figures/fullLength/CpGIs/done.{minCov}.txt",
        outFLnonCgiPlots="results/figures/fullLength/nonCGI/done.{minCov}.txt",
        out100CGI="results/figures/100bp5UTR/CpGIs/done.{minCov}.txt",
        out100noCGI="results/figures/100bp5UTR/nonCGI/done.{minCov}.txt",
        out4_600CGI="results/figures/400600bp5UTR/CpGIs/done.{minCov}.txt",
        out4_600nCGI="results/figures/400600bp5UTR/nonCGI/done.{minCov}.txt",
        # "results/figures/*/CpGIs/done.{minCov}.txt",
        # "results/figures/*/nonCGI/done.{minCov}.txt",
    conda: "r-env"
    log:
      logFLCGI="logs/figures/fullLength/CpGIs/plotLog_CGI_minCov.{minCov}.txt",
      logFLnoCGI="logs/figures/fullLength/nonCGI/plotLog_noCGI_minCov.{minCov}.txt",
#add  
      log100CGI="logs/figures/100bp5UTR/CpGIs/plotLog_CGI_minCov.{minCov}.txt",
      log100noCGI="logs/figures/100bp5UTR/nonCGI/plotLog_noCGI_minCov.{minCov}.txt",
      log4_600CGI="logs/figures/400600bp5UTR/CpGIs/plotLog_CGI_minCov.{minCov}.txt",
      log4_600noCGI="logs/figures/400600bp5UTR/nonCGI/plotLog_noCGI_minCov.{minCov}.txt"

    shell:
        """ 
        Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/plot_within_workflow.R --files "{input.L1fullLengthsOverlapsCGIs:q}" --outdir "results/figures/fullLength/CpGIs" 2> {log.logFLCGI} && touch {output.outFLCgiPlots}
        Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/plot_within_workflow.R --files "{input.L1fullLengthsOverlapsnCGIs:q}" --outdir "results/figures/fullLength/nonCGI" 2> {log.logFLnoCGI} && touch {output.outFLnonCgiPlots}

        Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/plot_within_workflow.R --files "{input.l1_100bp5UTRCGIs:q}" --outdir "results/figures/100bp5UTR/CpGIs" 2> {log.log100CGI} && touch {output.out100CGI}
        Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/plot_within_workflow.R --files "{input.l1_100bp5UTRnCGIs:q}" --outdir "results/figures/100bp5UTR/nonCGI" 2> {log.log100noCGI} && touch {output.out100noCGI}
        Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/plot_within_workflow.R --files "{input.L1_400_600bp5UTRCGIs:q}" --outdir "results/figures/400600bp5UTR/CpGIs" 2> {log.log4_600CGI} && touch {output.out4_600CGI}
        Rscript /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/plot_within_workflow.R --files "{input.L1_400_600bp5UTRnCGIs:q}" --outdir "results/figures/400600bp5UTR/nonCGI" 2> {log.log4_600noCGI} && touch {output.out4_600nCGI}

        """
# 'results/OverlapsL1PromoterDNAme/{samples}/fullLength/CpGIs/{samples}_5mCpG_5hmCpG_DNAme_mmflil1_8438_Overlaps_minCov{minCov}_CpGIslands.bed'


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

#convert 5mc coverage files to bismark coverage files

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/convert2BismarkCoverage.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/config/slurm --jobs unlimited --cores all --use-conda -np

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/convert2BismarkCoverage.smk --jobs 10 --cores all --keep-going --forceall -np


# rm -rf .snakemake logs results

parent_dir = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/"
configfile: parent_dir + "config/samples_5mC_bed_finalSet.yaml" #mouse samples

rule all:
    input: 
        expand("results/convert2BismarkCoverage/{samples}.cov.txt.gz", samples=config['samples'])

rule convert2BismarkCoverage:
    input: 
        lambda wildcards: config["samples"][wildcards.samples]
    output:
        "results/convert2BismarkCoverage/{samples}.cov.txt.gz"
    shell:
        """ 
        awk 'BEGIN {{ OFS = "\\t" }} {{print $1, $2, $3, $11, $12, $13}}' {input} | sort -k1,1 -k2,2n >  {output}
        """



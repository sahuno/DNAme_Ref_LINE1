#set parent direct
DIR=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1
mkdir -p ${DIR}/outputs/clean_analysis03262025/prepareDNAme




mamba activate snakemake
cd ${DIR}/outputs/clean_analysis03262025/prepareDNAme
snakemake -s ${DIR}/scripts/fullLengthL1_L1Base_DNAmeOverlaps.smk --workflow-profile ${DIR}/config/slurm --jobs unlimited --cores all --use-conda -np #remove -np to run for real
cd ..


mkdir -p ${DIR}/outputs/clean_analysis03262025/DNAme_L1Base_SummaryStats
cd ${DIR}/outputs/clean_analysis03262025/DNAme_L1Base_SummaryStats

sh /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/runL1Base_DNAmeSummaryStats.sh


#how to run tldr
# ssend numpyro tldr /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/scripts/run_tldr.sh

s044T_v14=/data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/hac_5mCG_5hmCG/results/mark_duplicates/044T_v14/044T_v14_modBaseCalls_sorted_dup.bam
s044N_v14=/data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/sup_5mCG_5hmCG/results/mark_duplicates/044N_v14/044N_v14_modBaseCalls_sorted_dup.bam
s009T2_009T1=/data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/hac_5mCG_5hmCG/mergedBams/009T2_009T1_n_2_modBaseCalls_sorted_dedup.bam
s009N=/data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/009N_rerun/modbasecalls/mergebam_modkit/results/mark_duplicates/009N/009N_modBaseCalls_sorted_dup.bam

# insert your parameters
# TUMOR=/data1/greenbab/projects/methyl_benchmark_spectrum/ONT_BSseq/ONT_DLP_1stPre/ont_bam/Spectrum-OV-051_T_final.bam
# NORMAL=/data1/greenbab/projects/methyl_benchmark_spectrum/ONT_BSseq/ONT_DLP_1stPre/ont_bam/Spectrum-OV-051_N_final.bam

REF=/data1/greenbab/database/hg38/v0/Homo_sapiens_assembly38.fasta
# singularity exec --no-home -B /data1/greenbab /data1/greenbab/projects/yelena_long_read/te.sif tldr --bams ${TUMOR},${NORMAL} --procs 10 --elts /tldr/ref/teref.human.fa --ref ${REF} --min_te_len 100 --max_cluster_size 100 --trdcol --detail_output --methylartist --keep_pickles

#run with sammples
# singularity exec --no-home -B /data1/greenbab /data1/greenbab/projects/yelena_long_read/te.sif tldr --bams ${s044T_v14},${s044N_v14} --procs 10 --elts /tldr/ref/teref.human.fa --ref ${REF} --min_te_len 100 --max_cluster_size 100 --trdcol --detail_output --methylartist --keep_pickles
singularity exec --no-home -B /data1/greenbab /data1/greenbab/projects/yelena_long_read/te.sif tldr --bams ${s009T2_009T1},${s009N} --procs 10 --elts /tldr/ref/teref.human.fa --ref ${REF} --min_te_len 100 --max_cluster_size 100 --trdcol --detail_output --methylartist --keep_pickles


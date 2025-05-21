singularity exec -B /data1/greenbab/ /data1/greenbab/users/ahunos/apps/containers/igv_latest.sif \
igver.py --bam /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/results/QSTAT_group_5mc/tmp/*.bw -r /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/igvr_resultsDTSampleUpDown.txt \
            -o "igv_of_DE/" \
                -mph 500 -od expand \
                    --genome 'mm10'
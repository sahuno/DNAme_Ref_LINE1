singularity run -B /data1/greenbab docker://shahcompbio/igv igver.py \
--bam /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/mark_duplicates/D-0-1_5000_4000/D-0-1_5000_4000_modBaseCalls_sorted_dup.bam \
/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/mark_duplicates/D-0-2_5000_4000/D-0-2_5000_4000_modBaseCalls_sorted_dup.bam \
/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/mark_duplicates/D-0-3_4000/D-0-3_4000_modBaseCalls_sorted_dup.bam \
/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/mark_duplicates/D-A-1_4000/D-A-1_4000_modBaseCalls_sorted_dup.bam \
/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/mark_duplicates/D-A-2_4000/D-A-2_4000_modBaseCalls_sorted_dup.bam \
/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/mark_duplicates/D-A-3_4000/D-A-3_4000_modBaseCalls_sorted_dup.bam \
-r /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/regions_UID100.txt \
-o /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/igv_screenShots \
-mph 500 -od squish \
--genome 'mm10'



### run on bigwigs
find /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/BedGraphsBigWigs/ -name '*_minCov10.bw'

singularity run -B /data1/greenbab docker://shahcompbio/igv igver.py \
--bam /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/BedGraphsBigWigs/D-0-3_4000/D-0-3_4000_5mCpG_5hmCpG_sortedBed_minCov10.bw \
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/BedGraphsBigWigs/D-A-1_4000/D-A-1_4000_5mCpG_5hmCpG_sortedBed_minCov10.bw \
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/BedGraphsBigWigs/D-A-3_4000/D-A-3_4000_5mCpG_5hmCpG_sortedBed_minCov10.bw \
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/BedGraphsBigWigs/D-0-2_5000_4000/D-0-2_5000_4000_5mCpG_5hmCpG_sortedBed_minCov10.bw \
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/BedGraphsBigWigs/D-A-2_4000/D-A-2_4000_5mCpG_5hmCpG_sortedBed_minCov10.bw \
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/BedGraphsBigWigs/D-0-1_5000_4000/D-0-1_5000_4000_5mCpG_5hmCpG_sortedBed_minCov10.bw \
-r /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/regions_UID100.txt \
-o /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/igv_screenShots/bw \
-mph 500 -od squish \
--genome 'mm10'



## run for all line1
#make 
awk '{print $1":"$2"-"$3"\t"$4}' /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/head_mmflil1_8438_noHeader.sorted.bed \
> /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/head_regions_mmflil1_8438_noHeader.sorted.bed

awk '{print $1":"$2"-"$3"\t"$4}' /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed \
> /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/regions_mmflil1_8438_noHeader.sorted.bed


singularity run -B /data1/greenbab docker://shahcompbio/igv igver.py \
--bam /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/BedGraphsBigWigs/D-0-3_4000/D-0-3_4000_5mCpG_5hmCpG_sortedBed_minCov10.bw \
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/BedGraphsBigWigs/D-A-1_4000/D-A-1_4000_5mCpG_5hmCpG_sortedBed_minCov10.bw \
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/BedGraphsBigWigs/D-A-3_4000/D-A-3_4000_5mCpG_5hmCpG_sortedBed_minCov10.bw \
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/BedGraphsBigWigs/D-0-2_5000_4000/D-0-2_5000_4000_5mCpG_5hmCpG_sortedBed_minCov10.bw \
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/BedGraphsBigWigs/D-A-2_4000/D-A-2_4000_5mCpG_5hmCpG_sortedBed_minCov10.bw \
/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/BedGraphsBigWigs/D-0-1_5000_4000/D-0-1_5000_4000_5mCpG_5hmCpG_sortedBed_minCov10.bw \
-r /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.100bp5UTR_sorted.bed \
-o /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/igv_screenShots/mmflil1_100bp \
-mph 500 -od squish \
--genome 'mm10'

# singularity run -B /data1/greenbab docker://shahcompbio/igv igver.py --help
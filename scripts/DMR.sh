
#purpose: dmr on all pairs of samples 5-Aza & DMSO

#look for all files
#find . -name "*_minCov10.bed.gz" | xargs realpath
# find . -name *_h_CG0_*
D01=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/prepareBedFiles/D-0-1_5000_4000/D-0-1_5000_4000_5mCpG_5hmCpG_sortedBed_minCov10.bed.gz
DA1=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/prepareBedFiles/D-A-1_4000/D-A-1_4000_5mCpG_5hmCpG_sortedBed_minCov10.bed.gz
D02=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/prepareBedFiles/D-0-2_5000_4000/D-0-2_5000_4000_5mCpG_5hmCpG_sortedBed_minCov10.bed.gz
DA2=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/prepareBedFiles/D-A-2_4000/D-A-2_4000_5mCpG_5hmCpG_sortedBed_minCov10.bed.gz
D03=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/prepareBedFiles/D-0-3_4000/D-0-3_4000_5mCpG_5hmCpG_sortedBed_minCov10.bed.gz
DA3=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/prepareBedFiles/D-A-3_4000/D-A-3_4000_5mCpG_5hmCpG_sortedBed_minCov10.bed.gz


mmflil1=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.head.bed
genomeChromSizes=/data1/greenbab/database/mm10/mm10.sorted.chrom.sizes
genomeFasta=/data1/greenbab/database/mm10/mm10.fa

modkit dmr multi \
  -s ${D01} D01 \
  -s ${DA1} DA1 \
  -s ${D02} D02 \
  -s ${DA2} DA2 \
  -s ${D03} D03 \
  -s ${DA3} DA3 \
  --out-dir dmr_D013_DA13_dir \
  -r ${mmflil1} \
  --ref ${genomeFasta} \
  --base C \
  -t 10 \
  -f \
  --log-filepath dmr_multiD013_DA13.log


#bedD01=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/prepareBedFiles/D-0-1_5000_4000/D-0-1_5000_4000_5mCpG_5hmCpG_sortedBed_minCov5.bed.gz
#bedDA1=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/prepareBedFiles/D-A-1_4000/D-A-1_4000_5mCpG_5hmCpG_sortedBed_minCov5.bed.gz
#modkit localise ${bedD01} --regions ${mmflil1} --genome-sizes ${genomeChromSizes} --chart "bedD01_mmflil1.html"

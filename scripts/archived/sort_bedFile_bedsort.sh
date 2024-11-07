#bedtools sort 

bedNotZipped=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/modkit/D-0-1_5000_4000/D-0-1_5000_4000_modpileup_combined.bed

sort -k1,1 -k2,2n $bedNotZipped > bedfileD015000_4000_modpileup_combined.sorted.bed

# sortBed -i $bedNotZipped > bedfileD015000_4000_modpileup_combined.sorted.bed

bedtools map -a /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed -b /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/bedfileD015000_4000_modpileup_combined.sorted.bed -c 11,11,11 -o mean,median,count -g /data1/greenbab/database/mm10/mm10.sorted.chrom.sizes > D-0-1_5000_4000_mmflil1_8438_sortBedBedtoolsV2.bed

# ssend numpyro lgphase /data1/greenbab/users/ahunos/apps/workflows/LINE1_insertions/scripts/sort_bedFile_bedsort.sh
sq=/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/CT/squire_te_fwd.tsv
#less $sq

chrs="chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chrX|chrY"
cut --complement -f 1 "$sq" | grep "L1" | awk -F'|' '{OFS="\t"; print $1, $2, $3, $4, $5, $6}' | grep -E $chrs > squireL1.bed

#create a file with only cordinates to test repeatmasker igv
cat squireL1.bed | cut -f 1-6 > squireL1_cordinatesOnly.bed
#head squireL1.bed

bedtools intersect -b squireL1.bed -a /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed -wa -wb > overlapsSquireL1.bed

# #get the complement of sq
# cut --complement -f 1 $sq | grep "L1" | head
# chrs="chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chrX, chrY"
# cut --complement -f 1 $sq | head grep chrs

# #can either 1. merge l 
# awk '{}' 


# cut --complement -f 1 $sq | split

# cut --complement -f 1 "$sq" | grep -E "chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chrX|chrY" | grep

# #all line 1 
# # cut --complement -f 1 "$sq" | grep "L1" | awk -F'|' '{OFS="\t"; print $1, $2, $3, $4, $5, $6}'  | grep -E "chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chrX|chrY" | head


# i want rna-cordinates `SquireL1.bed` that aligns ~80% with any full length l1 from l1 base
bedtools intersect -a SquireL1.bed -b l1Base_entireLength.bed -f 0.8 -wb > overlapsSquireRNA_l1Base.bed
less overlapsSquireRNA_l1Base.bed

# grep "1.21e+08" l1Base_entireLength.bed
# grep "1.21e+08" SquireL1.bed

#download L1Base databases
# Homo sapiens (Human): NCBI38

#create dir to download files into
mkdir -p database/L1BaseHsapiens
cd database/L1BaseHsapiens

#  Human Full-Length, Intact LINE-1 Elements [FLI-L1] (Ens84.38)
wget https://l1base.charite.de/BED/hsflil1_8438.bed 

#  Human ORF2 Intact LINE-1 Elements [ORF2-L1] (Ens84.38)
wget https://l1base.charite.de/BED/hsorf2l1_8438.bed

# Human Full-Length >4500nt LINE-1 Elements [FLnI-L1] (Ens84.38)
wget https://l1base.charite.de/BED/hsflnil1_8438_rm.bed


#sort bedfiles
for i in $(ls)
do 
Sname=$(basename -s .bed $i)
#remove header
sed '1d' ${i} | sort -k1,1 -k2,2n >  ${Sname}_noHeader.sorted.bed
# echo "sed '1d' ${i} | sort -k1,1 -k2,2n >  ${Sname}_noHeader.sorted.bed"
done

#return to database/ directory
cd ..


###mouse
mkdir -p L1BaseMmusculus
cd L1BaseMmusculus

wget https://l1base.charite.de/BED/mmflil1_8438.bed
wget https://l1base.charite.de/BED/mmorf2l1_8438.bed
wget https://l1base.charite.de/BED/mmflnil1_8438_rm.bed


for i in $(ls)
do 
Sname=$(basename -s .bed $i)
sed '1d' ${i} | sort -k1,1 -k2,2n >  ${Sname}_noHeader.sorted.bed
# echo "sed '1d' ${i} | sort -k1,1 -k2,2n >  ${Sname}_noHeader.sorted.bed"
done
# sort -k1,1 -k2,2n mmflil1_8438_noHeader.bed >  mmflil1_8438_noHeader.sorted.bed

#retun to parent dir
cd ../..

# l1File=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed
head $l1File > mmflil1_8438_noHeader.sorted.head.bed
less mmflil1_8438_noHeader.sorted.head.bed

awk 'BEGIN { OFS = "\t" } {
if ($6 == "+"){
    $2 = $2 + 400;
    $3 = $2 + 200
} else if ($6 == "-") {
    $3 = $3 - 400; 
    $2 = $3 - 200
} 
print $0
}' mmflil1_8438_noHeader.sorted.head.bed


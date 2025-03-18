#here is script to create the regions
#input is the bed file with the regions
#output is the `.txt` file for igv

input=$1
output=$2
awk '{print $1":"$2"-"$3"\t"$4}' ${input} > ${output}



#examples
#sh /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/igvScreenshots/createRegions.sh ../mmflil1_8438_noHeader.sorted.bed mmflil1_8438_noHeader.sorted.txt

#sh /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/igvScreenshots/createRegions.sh ../mmflil1_8438_noHeader.100bp5UTR_sorted.bed mmflil1_8438_noHeader.100bp5UTR_sorted.txt
regionsBed=$1
genes=$2

# regionsBed=/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/mergedSquireRepRNA/SQuireL1MostActive_ValidnNonValidReadCounts_Aggregated_FWD_REV.bed

# genes=/data1/greenbab/database/Gencode_Genes/gencode.v47.genes.annotation.sorted.gff3

bedtools intersect -u -a ${regionBed} -b ${genes} > intragenic.bed

  git clone https://github.com/gpertea/gffread
  cd gffread
  make release

#how to run tldr
#how to run
# sgpu numpyro l1_basecalling /data1/greenbab/projects/methyl_benchmark_spectrum/ONT_BSseq/ONT_DLP_1stPre/L1_DNAme/L1_basecalling_mapping.sh 4

#basecalling rese
#### dorado appraoch
OV_044T=/data1/greenbab/projects/methyl_benchmark_spectrum/data/raw/pod5/split_reads/results/pod5subset/split_by_sample_rate/044T_v14/sample_rate-4000.pod5
OV_044N=/data1/greenbab/projects/methyl_benchmark_spectrum/data/raw/pod5/split_reads/results/pod5subset/split_by_sample_rate/044N_v14/sample_rate-4000.pod5
# REF=/data1/greenbab/database/hg38/v0/Homo_sapiens_assembly38.fasta
L1_REF=/data1/greenbab/projects/DNAme_SPECTRUM/data/preprocessed/teref.ont.human.fa

#create chromSizes
cut -f1,2 ${L1_REF}.fai > ${L1_REF}.sizes.genome



#try basecalling with non-reference genomes
# dorado basecaller hac,5mCG_5hmCG@latest ${OV_044T} --device "cuda:all" --emit-sam --reference ${L1_REF} --verbose | samtools sort --threads 12 -o OV_044T_modBase.bam --write-index
# | samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} --write-index  2> {log} || dorado basecaller hac,5mCG_5hmCG@latest {input} --device {params.device} --emit-sam --reference {params.reference_genome} --verbose | samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} --write-index  2> {log}

dorado basecaller hac,5mCG_5hmCG@latest ${OV_044N} --device "cuda:all" --emit-sam --reference ${L1_REF} --verbose | samtools sort --threads 12 -o OV_044N_modBase.bam --write-index



#this gives reads that maps to this TE regions. why don't we do this for te insertions?

#results
OV_044T_BAM=/data1/greenbab/projects/DNAme_SPECTRUM/data/preprocessed/OV_044T/OV_044T_modBase.bam
OV_044N_BAM=/data1/greenbab/projects/DNAme_SPECTRUM/data/preprocessed/OV_044N/OV_044N_modBase.bam
# samtools view -H /data1/greenbab/projects/DNAme_SPECTRUM/data/preprocessed/OV_044T_modBase.bam
OV_044T_pileup=modkit_out/OV_044T/OV_044T_pileup.bed.gz
OV_044N_pileup=modkit_out/OV_044N/OV_044N_pileup.bed.gz




# modkit summary --threads 12 ${OV_044T_BAM} --log-filepath OV_044T_BAM.summary_log > OV_044T_BAM.summary_txt

# modkit pileup --threads 12 --bedgraph ${OV_044T_BAM} modkit_out/OV_044T/ --prefix OV_044T --cpg --combine-mods --ref ${L1_REF}
# modkit pileup --threads 12 --bedgraph ${OV_044T_BAM} modkit_out/OV_044T_BAM/ --prefix OV_044T_BAM --cpg --combine-strands --combine-mods --ref ${L1_REF}
# modkit pileup --threads 12 --bedgraph ${OV_044N_BAM} modkit_out/OV_044N/ --prefix OV_044N --cpg --combine-strands --combine-mods --ref ${L1_REF}


modkit pileup ${OV_044T_BAM} - \
  --cpg \
  --ref ${L1_REF} \
  --threads 12 \
  --log-filepath log.txt | bgzip -c > ${OV_044T_pileup}
tabix -p bed ${tumor_pileup}


##sort and convert to
for i in modkit_out/OV_044T_BAM/*.bedgraph
do
echo $i
SName=$(basename -s .bedgraph $i)
# echo $SName
    awk -v min_cov=10 '$5 > min_cov {print $1, $2, $3, $4}' "${i}" | sort -k1,1 -k2,2n >  ${i}.sorted.bedgraph
    bedGraphToBigWig ${i}.sorted.bedgraph ${L1_REF}.sizes.genome ${i}.sorted.bw
done


# l modkit_out/OV_044T_BAM/
for i in modkit_out/OV_044T_BAM/*_combined.bedgraph
do
echo $i
SName=$(basename -s .bedgraph $i)
# echo $SName
    awk -v min_cov=10 '$5 > min_cov {print $1, $2, $3, $4}' "${i}" | sort -k1,1 -k2,2n >  ${i}.sorted.bedgraph
    bedGraphToBigWig ${i}.sorted.bedgraph ${L1_REF}.sizes.genome ${i}.sorted.bw
done


for i in modkit_out/OV_044N/*.bedgraph
do
echo $i
# SName=$(basename -s .bedgraph $i)
# echo $SName
    awk -v min_cov=10 '$5 > min_cov {print $1, $2, $3, $4}' "${i}" | sort -k1,1 -k2,2n >  ${i}.sorted.bedgraph
    bedGraphToBigWig ${i}.sorted.bedgraph ${L1_REF}.sizes.genome ${i}.sorted.bw
done


# awk -v min_cov="10" '$5 > min_cov{{print $1, $2, $3, $4}}' "${i}" | sort -k1,1 -k2,2n - | bedGraphToBigWig - ${L1_REF}.sizes.genome ${i}.sorted.bw


awk -v min_cov="10" '$5 > min_cov{{print $1, $2, $3, $4}}' "modkit_out/OV_044T_BAM/*.bedgraph" | sort -k1,1 -k2,2n >  
# l - | tee *.bedgraph awk -v min_cov="10" '$5 > min_cov{{print $1, $2, $3, $4}}' "modkit_out/OV_044T_BAM/*.bedgraph" | sort -k1,1 -k2,2n >  


# # insert your parameters; fasta and fastq
# OV_3T2=/data1/greenbab/projects/DNAme_SPECTRUM/data/raw/SPECTRUM_III_050322/3T2/20220503_1329_3A_PAK09677_b17c678d/fastq_pass
# OV_3N=/data1/greenbab/projects/DNAme_SPECTRUM/data/raw/SPECTRUM_III_050322/3N/20220503_1704_2H_PAM56774_6e6d7872/fastq_pass
# OUTDIR=/data1/greenbab/projects/DNAme_SPECTRUM/data/raw/gather_fastq
# cat ${OV_3T2}/*.fastq > ${OUTDIR}/${OV_3T2}.fastq
# cat ${OV_3N}/*.fastq > ${OUTDIR}/${OV_3N}.fastq

# REF=/data1/greenbab/database/hg38/v0/Homo_sapiens_assembly38.fasta

# #this wouldn't work; cos i don't have acess to fast5
# nanopolish index -d fast5_files/ output.fastq

# # BAMLIST_TUMOR=/data1/greenbab/projects/methyl_benchmark_spectrum/ONT_BSseq/ONT_DLP_1stPre/tldr_outputs/scripts/sample_bams_tumorOnly.txt
# singularity exec --no-home -B /data1/greenbab /data1/greenbab/projects/yelena_long_read/te.sif tldr --bams ${TUMOR},${NORMAL} --procs 10 --elts /tldr/ref/teref.human.fa --ref ${REF} --min_te_len 100 --max_cluster_size 100 --trdcol



# # #try basecalling with non-reference genomes
# # dorado basecaller sup,5mCG_5hmCG@latest {input} --device {params.device} --emit-sam --reference {params.reference_genome} --verbose | samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} --write-index  2> {log} || dorado basecaller hac,5mCG_5hmCG@latest {input} --device {params.device} --emit-sam --reference {params.reference_genome} --verbose | samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} --write-index  2> {log}


#
for i in /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/modkit/*/*_combined.bed
do
echo $i
# SName=$(basename -s .bedgraph $i)
echo "processing file $SName"
bgzip -k ${i}
tabix -p bed ${i}.gz
done


modkit localise --regions /data1/greenbab/database/L1base2_NCBIm38/mmflil1_8438.bed --genome-sizes /data1/greenbab/database/mm10/chrom.sizes /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/modkit/D-0-1_5000_4000/D-0-1_5000_4000_modpileup_combined.bed.gz --min-coverage 10 --out-file D-0-1_5000_4000.localise.txt --name D-0-1_5000_4000_localise --chart D-0-1_5000_4000_localise_chart
modkit localise --regions /data1/greenbab/database/L1base2_NCBIm38/mmflil1_8438_UID1.bed --genome-sizes /data1/greenbab/database/mm10/chrom.sizes /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/modkit/D-0-1_5000_4000/D-0-1_5000_4000_modpileup_combined.bed.gz --min-coverage 10 --name D-0-1_5000_4000_mmflil1_8438_UID1.txt --chart D-0-1_5000_4000_mmflil1_8438_UID1.html



#### get methylation at repeats

rmskmm10=/data1/greenbab/database/RepeatMaskerDB/ucsc/mm10/rmsk_mm10.bed
awk '$4 ~ /^L1/' $rmskmm10 > /data1/greenbab/database/RepeatMaskerDB/ucsc/mm10/rmsk_L1_mm10.bed

awk '{print >> $4".txt"}' /data1/greenbab/database/RepeatMaskerDB/ucsc/mm10/rmsk_L1_mm10.bed

awk '$4 ~ /^L1/ {print > $4".txt"}' /data1/greenbab/database/RepeatMaskerDB/ucsc/mm10/rmsk_L1_mm10.bed

awk '$4 ~ /^L1MdMus/' /data1/greenbab/database/RepeatMaskerDB/ucsc/mm10/rmsk_L1_mm10.bed





D01_4000=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/modkit/D-0-1_4000/D-0-1_4000_modpileup_combined.bed.gz
modkit stats --regions ${rmskmm10} --out-table o44N_v14.tsv ${D01_4000}

modkit stats --regions /data1/greenbab/database/RepeatMaskerDB/ucsc/mm10/L1/L1Md_F.txt --out-table D01_4000_L1Md_F.tsv ${D01_4000}



### do for l1base
l1baseFL=/data1/greenbab/database/L1Base/L1base2_NCBIm38/mmflil1_8438.bed
sed '1d' $l1baseFL > mmflil1_8438_noHeader.bed 
modkit stats --regions /data1/greenbab/database/L1Base/L1base2_NCBIm38/mmflil1_8438_noHeader.bed --out-table D01_4000_mmflil1_8438.tsv ${D01_4000}

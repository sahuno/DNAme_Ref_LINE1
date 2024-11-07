find . -type f -name "*_modpileup_combined.bed.gz*" | xargs realpath > tmp.yaml && sed -i "1i samples:" -

find . -type f -name "*_modpileup_combined.bed.gz" | xargs realpath | tee samples_TriEpi_modkit_modpileup_combined.bed.gz.yaml && sed -i 's/^/    sampleName: /' samples_TriEpi_modkit_modpileup_combined.bed.gz.yaml && sed -i '1i samples:' samples_TriEpi_modkit_modpileup_combined.bed.gz.yaml


find . -type f -name "*_modpileup_combined.bed.gz*" | xargs realpath | tee samples_TriEpi_modkit_modpileup_combined.bed.gz.yaml && sed -i '1i samples:' samples_TriEpi_modkit_modpileup_combined.bed.gz.yaml



find . -type f -name "*_modpileup_combined.bed.gz" | xargs realpath | tee samples_TriEpi_modkit_modpileup_combined.bed.gz.yaml && sed -i 's/^/sampleName: /' samples_TriEpi_modkit_modpileup_combined.bed.gz.yaml 
sed -E 's|sampleName: (.*)|sampleName: $(basename \1)|e' samples_TriEpi_modkit_modpileup_combined.bed.gz.yaml 


# && sed -i '1i samples:' samples_TriEpi_modkit_modpileup_combined.bed.gz.yaml

| sed -i '1i samples:' samples_TriEpi_modkit_modpileup_combined.bed.gz.yaml > samples_TriEpi_modkit_modpileup_combined.bed.gz.final.yaml
# rm samples_TriEpi_modkit_modpileup_combined.bed.gz.yaml samples_TriEpi_modkit_modpileup_combined.bed.gz.final.yaml


sort -k1,1 -k2,2n /data1/greenbab/database/mm10/mm10.chrom.sizes >  /data1/greenbab/database/mm10/mm10.sorted.chrom.sizes

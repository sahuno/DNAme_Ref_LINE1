#cpu2
#npo
import glob
import polars as pl
import pandas as pd
import pyranges
import pyfaidx
import pyranges as pr


##load bed files
paths_bed = glob.glob("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/rerun_test/results/prepareBedFiles/*/*minCov10.bed", recursive=True)



beds = read_modkit_pileups(paths_bed).head(500000).collect()
bedp = beds.to_pandas()
p=pr.PyRanges(bedp)

mm10fullLengths = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/database/L1BaseMmusculus/mmflil1_8438_noHeader.sorted.bed"
read_csv()
##import full length 

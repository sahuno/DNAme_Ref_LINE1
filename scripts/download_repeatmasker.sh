#!/bin/bash

# Check input argument
if [[ "$1" != "human" && "$1" != "mouse" ]]; then
  echo "Usage: $0 [human|mouse]"
  exit 1
fi

# Set variables
RESOURCES_DIR="resources"

if [ "$1" == "human" ]; then
  SPECIES="human"
  ASSEMBLY="hg38"
elif [ "$1" == "mouse" ]; then
  SPECIES="mouse"
  ASSEMBLY="mm10"
fi

# Set output directory
OUTDIR="${RESOURCES_DIR}/${SPECIES}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Download annotation and alignment files
echo "Downloading RepeatMasker files for $SPECIES ($ASSEMBLY)..."
wget -N "https://www.repeatmasker.org/genomicDatasets/${ASSEMBLY}.fa.out.gz"
wget -N "https://www.repeatmasker.org/genomicDatasets/${ASSEMBLY}.fa.align.gz"

# Unzip
gunzip -f "${ASSEMBLY}.fa.out.gz"
gunzip -f "${ASSEMBLY}.fa.align.gz"

echo "Downloaded and extracted to: $OUTDIR"



#usage
# chmod +x download_repeatmasker.sh

# # For human
# ./download_repeatmasker.sh human

# # For mouse
# ./download_repeatmasker.sh mouse
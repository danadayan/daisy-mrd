#!/bin/bash

set -euo pipefail

# Check number of arguments
if [ "$#" -ne 3 ]; then
    echo "Error: Expected exactly 3 arguments, got $#." >&2
    echo "Usage: $0 [input_vcf] [output_vcf] [vep cache directory]" >&2
    exit 1
fi

INPUT_VCF="$1"
OUTPUT_VCF="$2"
VEP_CACHE_HOST="$3"

# Resolve absolute paths
INPUT_VCF="$(realpath "$INPUT_VCF")"
OUTPUT_VCF="$(realpath -m "$OUTPUT_VCF")"

# Check input exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF does not exist: $INPUT_VCF" >&2
    exit 1
fi

# Split paths into directories + basenames
INPUT_DIR="$(dirname "$INPUT_VCF")"
INPUT_BASE="$(basename "$INPUT_VCF")"

OUTPUT_DIR="$(dirname "$OUTPUT_VCF")"
OUTPUT_BASE="$(basename "$OUTPUT_VCF")"

# Create output directory if needed
mkdir -p "$OUTPUT_DIR"

# Helper: insert suffix before .vcf or .vcf.gz
add_suffix() {
    local path="$1"
    local suffix="$2"

    if [[ "$path" == *.vcf.gz ]]; then
        echo "${path%.vcf.gz}${suffix}.vcf.gz"
    elif [[ "$path" == *.vcf ]]; then
        echo "${path%.vcf}${suffix}.vcf"
    else
        echo "Error: Expected a .vcf or .vcf.gz file name: $path" >&2
        exit 1
    fi
}

# Intermediate filenames
FILTERED_VCF="$(add_suffix "$OUTPUT_VCF" "_filtered")"
#FILTERED_VEP_ANNOTATED_VCF="$(add_suffix "$FILTERED_VCF" "_vep")"

#echo "FILTERED_VCF=$FILTERED_VCF"
#echo "FILTERED_VEP_ANNOTATED_VCF=$FILTERED_VEP_ANNOTATED_VCF"
#echo "FILTERED_VEP_GNOMAD_ANNOTATED_VCF=$FILTERED_VEP_GNOMAD_ANNOTATED_VCF"

# Step 1: keep only PASS variants
bcftools view -f PASS "$INPUT_VCF" -o "$FILTERED_VCF"

# Docker paths
#VEP_CACHE_HOST="/home/avraham/MaruvkaLab/pedaml/vep_data"
VEP_CACHE_CONTAINER="/input_vep_cache"

# Mount the output directory so Docker can read/write all intermediate/output files
docker run \
  --user "$(id -u):$(id -g)" \
  --rm \
  -v "$VEP_CACHE_HOST:$VEP_CACHE_CONTAINER" \
  -v "$OUTPUT_DIR:/work" \
  ensemblorg/ensembl-vep \
  vep \
  --input_file "/work/$(basename "$FILTERED_VCF")" \
  --output_file "/work/$(basename "$OUTPUT_BASE")" \
  --vcf \
  --cache \
  --offline \
  --dir_cache "$VEP_CACHE_CONTAINER" \
  --assembly GRCh38 \
  --hgvs \
  --everything \
  --no_stats \
  --symbol \
  --canonical \
  --pick \
#  --custom "$VEP_CACHE_CONTAINER/af-only-gnomad.hg38.vcf.gz,gnomAD,vcf,exact,0,AF,AC"

rm "$FILTERED_VCF"
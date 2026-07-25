#!/bin/bash
# ==============================================================================
# SCRIPT NAME: plink_to_vcf_prep.sh
# Pre-imputation data preparation: PLINK to VCF conversion (Autosomes only)
# ==============================================================================

# Define your input and output paths generically
INPUT_PREFIX="/path/to/input/my_dataset"
OUTPUT_PREFIX="/path/to/output/my_dataset_processed"

# Create output parent directory automatically if it does not exist
mkdir -p $(dirname "$OUTPUT_PREFIX")

# =============
# 01. plink2vcf
# =============

# 1. Convert PLINK binary to raw VCF text (Autosomes 1-22 only)
plink --bfile "$INPUT_PREFIX" --chr 1-22 --recode vcf --snps-only just-acgt --out "$OUTPUT_PREFIX"

# Converts flat text layout into an optimization-friendly block-gzip format (-Oz)
bcftools view "$OUTPUT_PREFIX.vcf" -Oz -o "$OUTPUT_PREFIX.vcf.gz"

# 3. Sorts spatial order by genomic position across all target autosomes
bcftools sort "$OUTPUT_PREFIX.vcf.gz" -Oz -o "$OUTPUT_PREFIX.sort.vcf.gz"

# 4. Builds a CSI/BAI index tracking genomic chunks for downstream analysis tools
bcftools index "$OUTPUT_PREFIX.sort.vcf.gz"

# 5. Extract text, prepend 'chr' string to chromosomes, and compress via bgzip
zcat "$OUTPUT_PREFIX.sort.vcf.gz" | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | bgzip > "$OUTPUT_PREFIX.sort.chr.vcf.gz"

# 6. Finalize files and generate spatial coordinate index table via tabix
tabix "$OUTPUT_PREFIX.sort.chr.vcf.gz"


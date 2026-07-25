#!/bin/bash
# ==============================================================================
# Script Name: run_smr_pipeline.sh
# Description: Genome-wide bash automation script for SMR 
#             (Summary-data-based Mendelian Randomization) across chromosomes 1 to 22.
# ==============================================================================

# ==============================================================================
# 1. USER CONFIGURATION PANEL
# ==============================================================================

# --- Executable Configuration ---
# Path to your local SMR tool binary file
SMR_EXEC="path/to/smr_binary_directory/smr"

# --- Reference & Input Data Paths ---
BFILE_PATH="path/to/reference/genotype_prefix"   # e.g., 1000 Genomes reference (prefix only)
GWAS_FILE="path/to/gwas_summary_stats.ma"        # GWAS summary dataset file (.ma)
EQTL_DIR="path/to/beqtl_summary_directory"       # Root directory holding binary eQTL data
OUTPUT_DIR="outputs/SMR_results"                 # Target path to save generated records

# --- File Naming Prefixes ---
EQTL_PREFIX="YOUR_QTL_FILE_PREFIX_chr"           # String pattern before the chromosome number
OUT_PREFIX="YOUR_OUTPUT_FILE_PREFIX_chr"         # String pattern for your output files

# ==============================================================================
# 2. PIPELINE INITIALIZATION & VALIDATION
# ==============================================================================

# Ensure target storage output directory exists
mkdir -p "$OUTPUT_DIR"

# Verify that the SMR software binary is executable before launching the loop
if [ ! -x "$SMR_EXEC" ]; then
    echo "Error: SMR executable not found or not executable at: $SMR_EXEC"
    exit 1
fi

# ==============================================================================
# 3. GENOME-WIDE CHROMOSOME EXECUTION LOOP (1 to 22)
# ==============================================================================
for chr in {1..22}; do
    echo "========================================="
    echo "Starting SMR analysis for Chromosome: $chr"
    echo "========================================="

    # Dynamically build file strings using parameters
    CURRENT_BEQTL="${EQTL_DIR}/${EQTL_PREFIX}${chr}"
    CURRENT_OUT="${OUTPUT_DIR}/${OUT_PREFIX}${chr}"

    # Execute SMR tool binary program
    "$SMR_EXEC" \
        --bfile "$BFILE_PATH" \
        --gwas-summary "$GWAS_FILE" \
        --beqtl-summary "$CURRENT_BEQTL" \
        --peqtl-smr "$P_EQTL_SMR" \
        --out "$CURRENT_OUT"

    echo "Finished SMR analysis for Chromosome: $chr"
done

echo "========================================="
echo "[SUCCESS] SMR genome-wide analysis complete for all chromosomes."
echo "========================================="


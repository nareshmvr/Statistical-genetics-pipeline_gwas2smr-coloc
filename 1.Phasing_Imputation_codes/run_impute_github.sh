#!/bin/bash
# ==============================================================================
# SCRIPT NAME: run_impute.sh
# High-Throughput Imputation Workspace using Impute5 (Autosomes 1-22)

# ==============================================================================
# 1. TOOL PATHS
# Adjust these paths to point to your local Impute5 installation binaries.
# ==============================================================================

CHUNKER="./imp5Chunker_v1.2.0_static"
IMPUTE5="./impute5_v1.2.0_static"

# ==============================================================================
# 2. USER-DEFINED INPUT/OUTPUT DIRECTORIES
# Change these relative or absolute paths to match your system architecture.
# ==============================================================================
# Input paths
REF_DIR="/path/to/reference/1KG3_AF/"                      # Path to reference panel (e.g., 1000 Genomes BCF files)
MAP_DIR="/path/to/reference/shapeit_REF/b38"               # Path to genetic maps (e.g., SHAPEIT4/Impute5 b38 maps)
IP_PHASED_DATA="/path/to/project/data/02.shapeit4_phased_data/"        # Path to your target cohort's phased BCF data

# Output paths
COORD_DIR="/path/to/project/output/coordinates_chr"       # Directory to save generated chunk coordinates 
OUT_IMPUTE="/path/to/project/output/imputed_data"         # Directory to save final imputed BCF files

# Create output folders if they do not exist
mkdir -p "$COORD_DIR" "$OUT_IMPUTE"

# Define the template name of your cohort/dataset
COHORT="your_cohort_name"

# ==============================================================================
# 3. PIPELINE EXECUTION (Chromosomes 1-22)
# ==============================================================================
for CHR in {1..22}; do
    echo "=================================================="
    echo "Running Imputation Workspace for Chromosome: chr${CHR}"
    echo "=================================================="

    # Set up dynamic paths for the current chromosome
    REFERENCE="${REF_DIR}/1000GP_chr${CHR}_sorted.bcf"  # e.g., 1kg3_chr1.bcf
    MAP="${MAP_DIR}/chr${CHR}.b38.gmap.gz"              # e.g., chr1.b38.gmap.gz
    IP_PHASED="${IP_PHASED_DATA}/${COHORT}_chr${CHR}_ref_phased.bcf" # e.g., sample1_chr1_phased.bcf
    OUT_COORDINATES="${COORD_DIR}/coordinates_chr${CHR}.txt"#e.g.,

    # Step A: Run imp5Chunker to split the chromosome into manageable chunks
    $CHUNKER --h "$REFERENCE" --g "$IP_PHASED" --r "chr${CHR}" --o "$OUT_COORDINATES"
    echo "Generated chunk coordinate file: $OUT_COORDINATES"

    # Verify coordinate file exists and has data before proceeding
    if [ ! -s "$OUT_COORDINATES" ]; then
        echo "WARNING: Coordinate file empty or missing for chr${CHR}. Skipping."
        continue
    fi

    # Step B: Read chunks dynamically into memory instead of opening the file thousands of times
    # This single while loop replaces the slow nested awk/for loops
    while read -r col1 col2 buffer_reg impute_reg rest; do
        
        # Skip potential empty lines or comments
        if [ -z "$impute_reg" ] || [ -z "$buffer_reg" ]; then
            continue
        fi

        echo " -> Processing Chunk: Impute=$impute_reg | Buffer=$buffer_reg"

        # Step C: Execute Impute5 per chunk
        $IMPUTE5 --h "$REFERENCE" --m "$MAP" --g "$IP_PHASED" --r "$impute_reg" --buffer-region "$buffer_reg" --o "${OUT_IMPUTE}/${impute_reg}_imputed.bcf"

    done < "$OUT_COORDINATES"

done


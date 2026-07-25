#!/usr/bin/env 

# ==============================================================================
# SCRIPT NAME: run_phasing.sh
# DESCRIPTION: Automates genome phasing across chromosomes 1-22 using SHAPEIT4.
# ==============================================================================

THREADS=32
SAMPLE="sample_id" # Example: NIA1234, project_alpha, etc.

# Directories (Replace with your actual paths)
MAP_DIR="/path/to/genetic_maps"       # e.g., /home/user/ref/shapeit_b38
REF_DIR="/path/to/reference_panels"   # e.g., /home/user/ref/1000Genomes
DATA_DIR="/path/to/input_data"       # e.g., /home/user/project/raw_vcfs
OUT_DIR="/path/to/output"             # e.g., /home/user/project/phased_results

# --- LOOP THROUGH CHROMOSOMES ---
for CHR in {1..22}; do
    echo "Processing Chromosome ${CHR}..."

    # Define File Paths (Examples based on standard naming conventions)
    MAP="${MAP_DIR}/chr${CHR}.b38.gmap.gz"      # e.g., chr1.b38.gmap.gz
    REF="${REF_DIR}/panel_chr${CHR}.vcf.gz"     # e.g., 1kg3_chr1.vcf.gz
    INP="${DATA_DIR}/${SAMPLE}.vcf.gz"          # e.g., sample1.vcf.gz
    
    OUT="${OUT_DIR}/${SAMPLE}_chr${CHR}.bcf"    # e.g., sample1_chr1.bcf
    LOG="${OUT_DIR}/${SAMPLE}_chr${CHR}.log"    # e.g., sample1_chr1.log

    # Run Phasing- Computes haplotype phasing utilizing target panels and matching genetic maps
    shapeit4 --input "$INP" --reference "$REF" --map "$MAP" --region "chr${CHR}" --output "$OUT" --log "$LOG" --thread $THREADS

    # Run Indexing-Builds CSI index for rapid downstream programmatic querying (.bcf.csi)
    bcftools index "$OUT" --thread $THREADS
done

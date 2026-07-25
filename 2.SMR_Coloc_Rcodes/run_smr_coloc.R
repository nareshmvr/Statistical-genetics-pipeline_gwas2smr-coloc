#!/usr/bin/env Rscript
# ==============================================================================
# Script Name: run_smr_coloc.R
# Description: Automates chromosome-by-chromosome SMR-filtering and Coloc 
#              colocalization analysis across a genome-wide scale.
# ==============================================================================

# Load required libraries
library(dplyr)
library(coloc)
library(data.table)
library(openxlsx)

# --- 1. CONFIGURATION & USER INPUTS ---
# Update these paths or parameter values to fit your directory structure
DIR_SMR   <- "path/to/SMR_results/"       # Folder containing .smr files
DIR_EQTL  <- "path/to/cis_eQTL_data/"     # Folder containing eQTL .txt files
FILE_GWAS <- "path/to/gwas_sumstats.tsv"  # Path to your clean GWAS summary stats
DIR_OUT   <- "outputs/"                   # Folder where results will be saved

# Parameter settings
P_SMR_CUTOFF  <- xx 			  # eg.. 5e-8  SMR significance cutoff
P_HEIDI_CUTOFF <- xx		  	  # eg ..5e-2  HEIDI p-value threshold (p > threshold passes)
WINDOW_SIZE   <- xxxxx                    #  eg ..500000 500kb window around the gene probe
MIN_SNPS      <- xx                       # eg..50     Minimum overlapping SNPs required

# Study-specific sample sizes (Update these for your specific study!)
GWAS_N    <- xxx               # eg .. 25000 Case-control sample size
GWAS_S    <-xxx                # eg .. 0.60 Case proportion (cases / total N)
EQTL_N    <-xxx                # eg .. 1500 eQTL sample size
GWAS_TYPE <- "xx"   	       # e.g., "cc" for case-control, "quant" for continuous traits
EQTL_TYPE   <- "xx"            # e.g., "quant" for continuous molecular QTL traits

# Create output directory if it doesn't exist
if(!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

# --- 2. MAIN LOOP (CHROMOSOMES 1 TO 22) ---
all_chr_results <- list()

# Load GWAS once outside the loop to save massive processing time
message("Loading GWAS data...")
gwas <- fread(FILE_GWAS)

for(chr in 1:22) {
  message(paste("=== Processing chromosome:", chr, "==="))
   
  # Construct paths dynamically using configured directories
  smr_path  <- file.path(DIR_SMR, paste0("SMR_PREFIX_chr", chr, ".smr"))
  eqtl_path <- file.path(DIR_EQTL, paste0("EQTL_PREFIX_chr", chr, ".txt"))

  # Validation check
  if(!file.exists(smr_path) | !file.exists(eqtl_path)) {
    warning(paste("Files missing for chr", chr, "- skipping chromosome."))
    next
  }

  # Load and filter SMR results
  smr_res <- fread(smr_path) %>% filter(p_SMR < P_SMR_CUTOFF & p_HEIDI > P_HEIDI_CUTOFF)
  
  if(nrow(smr_res) == 0) {
    message(paste("No genes passed SMR filters for chromosome:", chr))
    next
  }
  
  # Load eQTL dataset for current chromosome
  eqtl <- fread(eqtl_path)
  message(paste("Loaded SMR genes:", nrow(smr_res), "| eQTL rows:", nrow(eqtl)))

  # --- 3. COLOC ANALYSIS PER PRIORITIZED GENE ---
  for(i in 1:nrow(smr_res)) {
    gene_name <- smr_res$probeID[i]
    chr_hit   <- smr_res$ProbeChr[i]
    pos_hit   <- smr_res$Probe_bp[i]

    # Subset datasets to local genomic window
    gwas_sub <- gwas %>% filter(CHR == chr_hit, BP >= pos_hit - WINDOW_SIZE, BP <= pos_hit + WINDOW_SIZE)
    eqtl_sub <- eqtl %>% filter(Probe == gene_name)
     
    # Align and deduplicate SNPs between datasets
    gwas_clean <- gwas_sub %>% filter(SNP %in% intersect(gwas_sub$SNP, eqtl_sub$SNP)) %>%  group_by(SNP) %>% slice_min(P, with_ties = FALSE) %>% ungroup() %>% arrange(SNP)
    
    eqtl_clean <- eqtl_sub %>% filter(SNP %in% gwas_clean$SNP) %>% group_by(SNP) %>% slice_min(p, with_ties = FALSE) %>% ungroup() %>% arrange(SNP)

    # Skip if too few SNPs overlap
    if(nrow(gwas_clean) < MIN_SNPS) next

    # Format Coloc datasets
    dataset1 <- list(snp = gwas_clean$SNP, pvalues = gwas_clean$P, MAF = gwas_clean$FRQ, s = GWAS_S, N = GWAS_N, type = GWAS_TYPE)
    dataset2 <- list(snp = eqtl_clean$SNP, pvalues = eqtl_clean$p, type = EQTL_TYPE, N = EQTL_N, MAF = eqtl_clean$Freq)  
    
    # Run Colocalization
    res <- coloc.abf(dataset1, dataset2)

    # Extract lead variant data
    lead_snp_idx <- which.min(gwas_clean$P)
    lead_rsid    <- gwas_clean$SNP[lead_snp_idx]
    lead_pos     <- gwas_clean$BP[lead_snp_idx]

    # Combine results
    res_row           <- as.data.frame(t(res$summary))
    res_row$Gene      <- gene_name
    res_row$Chr       <- chr
    res_row$Lead_SNP  <- lead_rsid
    res_row$Lead_BP   <- lead_pos
    res_row$NSNPs     <- nrow(gwas_clean)
  
    all_chr_results[[paste0(gene_name, "_", chr)]] <- res_row
  }
}

# --- 4. EXPORT COMPILED RESULTS ---
if(length(all_chr_results) > 0) {
  final_table <- do.call(rbind, all_chr_results) %>% arrange(desc(PP.H4.abf))
  
  # File exports
  tsv_output  <- file.path(DIR_OUT, "SMR_Coloc_Results_Summary.tsv")
  xlsx_output <- file.path(DIR_OUT, "SMR_Coloc_Results_Summary.xlsx")
  
  fwrite(final_table, tsv_output, sep = "\t")
  write.xlsx(final_table, xlsx_output)
  
  message(paste("Analysis complete! Results successfully saved to:", DIR_OUT))
} else {
  warning("No significant colocalizations found across any chromosomes.")
}


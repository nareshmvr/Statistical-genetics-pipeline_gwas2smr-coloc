# ==============================================================================
# SCRIPT: Multi-QTL Integrated Manhattan Plot
# DESCRIPTION: Merges GWAS summary statistics with colocalization data to 
#              generate a publication-ready Manhattan plot using 'topr'.
# ==============================================================================

# 1. Load Required Libraries ---------------------------------------------------
library(data.table)
library(topr)
library(ggplot2)
library(ggrepel)

# 2. Define File Paths (Modify these for your project) -------------------------
coloc_file <- "path/to/coloc_rsid_hits_final.txt"
gwas_file  <- "path/to/gwas_summary_stats.tsv"
output_plot <- "path/to/output_manhattan_plot.png"

# 3. Load Datasets -------------------------------------------------------------
# Colocalization data must contain: rsid, Gene, Multi_QTL
coloc_data <- fread(coloc_file)

# GWAS data must contain mapping information and p-values
gwas_data  <- fread(gwas_file, select = c("SNP", "CHR", "BP", "P"))

# 4. Format and Merge Data -----------------------------------------------------
# Rename GWAS columns to match 'topr' package requirements
setnames(gwas_data, old = c("SNP", "BP", "CHR"), new = c("ID", "POS", "CHROM"))

# Merge datasets keeping all GWAS variants
merged_data <- merge(gwas_data, coloc_data, by.x = "ID", by.y = "rsid", all.x = TRUE)

# Sort by p-value ascending (strongest associations first)
merged_data <- merged_data[order(merged_data$P)]

# 5. Subset Data by QTL Classification ----------------------------------------
# Single-modality QTL subsets
eQTL     <- merged_data[Multi_QTL == "eQTL", ]
caQTL    <- merged_data[Multi_QTL == "CaQTL", ]
lncQTL   <- merged_data[Multi_QTL == "lncQTL", ]
mQTL     <- merged_data[Multi_QTL == "mQTL", ]
pQTL     <- merged_data[Multi_QTL == "pQTL", ]
sQTL     <- merged_data[Multi_QTL == "sQTL", ]
metabQTL <- merged_data[Multi_QTL == "metabQTL", ]

# Variants showing zero colocalization
Non_significant_loci <- merged_data[is.na(Multi_QTL), ]

# Multi-QTL subsets (variants overlapping multiple biological modalities)
e_ca_m_QTL  <- merged_data[Multi_QTL == "e_ca_m_QTL", ] # eg., coloc hits found in eQTL,mQTL and caQTL
e_m_QTL     <- merged_data[Multi_QTL == "e_m_QTL", ]    #eg., coloc hits found in eQTL,mQTL and caQTL  
e_p_me_QTL  <- merged_data[Multi_QTL == "e_p_me_QTL", ] #e.g., coloc hits found in eQTL,pQTL and meQTL
e_s_QTL     <- merged_data[Multi_QTL == "e_s_QTL", ]    #e.g., coloc hits found in eQTL and sQTL
m_s_QTL     <- merged_data[Multi_QTL == "m_s_QTL", ]     #e.g., coloc hits found in mQTL and seQTL
p_metab_QTL <- merged_data[Multi_QTL == "p_metab_QTL", ]  #e.g., coloc hits found in pQTL and meQTL


# Combine all intersecting multi-QTL tracks into a single visual tier
multi_QTL_combined <- rbind(e_ca_m_QTL, e_m_QTL, e_p_me_QTL, e_s_QTL, m_s_QTL, p_metab_QTL)

p <- manhattan(
  list(eQTL, caQTL, lncQTL, mQTL, pQTL, sQTL, metabQTL, multi_QTL_combined, Non_significant_loci),
  color = colors_list,
  legend_labels = legend_names,
  ntop=14,                       # Top variants to explicitly display
  use_shades = TRUE,             # Alternate background shades per chromosome
  shades_alpha = 0.4,            
  sign_thresh = c(5e-08, 5e-05), # Genome-wide and suggestive thresholds
  sign_thresh_color = c("red", "blue"), 
  ymax = 18,                     # Cap visual height for scale readability
  annotate = 5e-04,              # P-value cutoff for text labels
  label_fontface = "bold",
  alpha = 0.7
)

# 8. Refine Plot Layout and Typography ----------------------------------------
p_final <- p + labs(title    = "Study Manhattan Plot",subtitle = "Colocalized single and multi-QTL hits across molecular layers (PP.H4 > 0.70)") +theme(plot.title    = element_text(size = 18, face = "bold"),legend.text   = element_text(size = 12),legend.title  = element_text(size = 14, face = "bold"),legend.position = "right") + guides(color = guide_legend(title = "QTL Category"))

# 9. Save Visualization -------------------------------------------------------

ggsave(filename = output_plot,plot = p_final,width = 14, height = 8,dpi = 300)

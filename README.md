#Scripts for pre-GWAS data processing (VCF conversion, SHAPEIT4 phasing, IMPUTE5 imputation), followed by post-GWAS SMR, ColocR
integration, and Manhattan data visualization

## Repository Structure
### 1. Phasing & Imputation Codes
Located in `1.Phasing_Imputation_codes/`
* `run_plink_to_vcf_conversion_github.sh`: Converts PLINK binary formats to VCF.
* `run_phase_shapeit4_github.sh`: Performs haplotype phasing using SHAPEIT4.
* `run_impute_github.sh`: Runs imputation steps.
* ### 2. SMR & Colocalization
Located in `2.colocR_codes/`
* `run_smr_pipeline.sh`: Runs the SMR analysis pipeline.
* `run_smr_coloc.R`: Performs Summary-data-based Mendelian Randomization (SMR) and genetic colocalization using the `coloc` R package.
## 3. Plots
Located in `3.Plots/`
* `topr_manhattan_coloc_map_github.r`: Generates publication-ready Manhattan plots and regional colocalization maps using `topr`.

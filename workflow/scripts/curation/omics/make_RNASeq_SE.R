# workflow/scripts/curation/omics/make_RNASeq_SE.R
snakemake@source("../helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(fs)
  library(SummarizedExperiment)
  library(dplyr)
})

message("🚀 Loading RNASeq matrix: ", snk$input$raw)

# 1. Build the base SE from the TPM matrix
# Note: feature_col must match the row name column in your rna_tpm_hugo.csv
rna <- wide_assay_tsv_to_se(
  path = snk$input$raw,
  datatype = "RNASeq",
  feature_col = "hugo.id", 
  row_data_name = "hugo.id"
)
se <- rna$se

# 2. Load and Attach Local Metadata
message("📂 Attaching local feature metadata: ", snk$input$feature_meta)
feature_df <- read.delim(
  snk$input$feature_meta, 
  sep = "\t", 
  header = TRUE, 
  check.names = FALSE, 
  stringsAsFactors = FALSE
)

# 1. Perform the match
matched_indices <- match(rownames(se), feature_df$hugo.id)
# 2. Reorder with drop = FALSE to force R to keep the data.frame structure
feature_df_ordered <- feature_df[matched_indices, , drop = FALSE]
# 3. Explicitly set rownames
# This is where your error was occurring; drop = FALSE ensures dimensions exist
rownames(feature_df_ordered) <- rownames(se)

# Update rowData
rowData(se) <- S4Vectors::DataFrame(feature_df_ordered)

# 3. Add Experiment Metadata
metadata(se) <- list(
  name = "UHN Breast PDX",
  lab = "Dave Cescon Lab",
  genome = "GRCh38",
  source_metadata = snk$input$feature_meta
)

# 4. Save Outputs
dir_create(path_dir(snk$output$se))
message("💾 Saving SE object to: ", snk$output$se)
saveRDS(se, snk$output$se)

# Save the clean matrix for compatibility
write_matrix_tsv(
  matrix = rna$matrix,
  row_ids = rownames(se),
  row_id_name = "hgnc_symbol",
  output_path = snk$output$matrix
)

message("✅ RNASeq SE Curation Complete (Using Local Metadata).")
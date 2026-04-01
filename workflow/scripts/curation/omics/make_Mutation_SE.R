snakemake@source("../helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(fs)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(dplyr)
})

message("🚀 Loading detailed mutation calls: ", snk$input$raw)

# 1. Load the Wide TSV
# check.names=FALSE is critical for IDs like '01_1104_945W_RES'
mut_wide <- read.delim(
  snk$input$raw, 
  sep = "\t", 
  header = TRUE, 
  check.names = FALSE, 
  stringsAsFactors = FALSE
)

# 2. Extract the Matrix
# Use hugo.id as rownames and remove the first column from the matrix
gene_ids <- mut_wide$hugo.id
mutation_mat <- as.matrix(mut_wide[, -1, drop = FALSE])
rownames(mutation_mat) <- gene_ids
sample_ids <- colnames(mutation_mat)

# 3. Clean strings (Optional but recommended)
# Replaces NA with empty strings to keep the matrix clean
mutation_mat[is.na(mutation_mat)] <- ""

# 4. Prepare Metadata
message("📂 Aligning feature metadata: ", snk$input$feature_meta)
feature_df <- read.delim(snk$input$feature_meta, sep = "\t", check.names = FALSE)

# Match rows and maintain dimensions
matched_idx <- match(rownames(mutation_mat), feature_df$hugo.id)
feature_df_ordered <- feature_df[matched_idx, , drop = FALSE]
rownames(feature_df_ordered) <- rownames(mutation_mat)

# 5. Build SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(variant_type = mutation_mat), # Detailed types preserved here
  rowData = S4Vectors::DataFrame(feature_df_ordered),
  colData = S4Vectors::DataFrame(sample_annotation_from_ids(sample_ids))
)

S4Vectors::metadata(se) <- list(
  datatype = "mutation(tumor-only WES)",
  genome = "GRCh38",
  source = snk$input$raw
)

# 6. Save Outputs
dir_create(path_dir(snk$output$se))
message("💾 Saving SE to: ", snk$output$se)
saveRDS(se, snk$output$se)

# Save matrix back to TSV for compatibility
write.table(
  mut_wide, 
  snk$output$matrix, 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)

message("✅ Mutation SE Curation Complete (Detailed Profile).")
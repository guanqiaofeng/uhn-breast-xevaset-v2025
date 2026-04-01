# workflow/scripts/curation/omics/make_CNV_SE.R
snakemake@source("../helpers.R")
snk <- parse_snakemake()

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(S4Vectors)
  library(fs)
})

message("🚀 Loading CNV matrix: ", snk$input$raw)

# 1. Use the helper function (now that we know it returns a list)
cnv_data <- wide_assay_tsv_to_se(
  path = snk$input$raw,
  datatype = "CNV(tumor-only WES)",
  feature_col = "hugo.id",
  assay_name = "copy_number"
)
se <- cnv_data$se

# 2. Attach Feature Metadata (using the TSV logic we perfected)
message("📂 Attaching CNV feature metadata...")
feature_df <- read.delim(snk$input$feature_meta, sep = "\t", check.names = FALSE)

# Match rows and use drop = FALSE to prevent dimension loss
matched_idx <- match(rownames(se), feature_df$hugo.id)
feature_df_ordered <- feature_df[matched_idx, , drop = FALSE]
rownames(feature_df_ordered) <- rownames(se)

rowData(se) <- S4Vectors::DataFrame(feature_df_ordered)

# 3. Add CNV-specific Metadata
metadata(se) <- list(
  name = "UHN Breast PDX CNV",
  genome = "GRCh38",
  type = "Relative Copy Number"
)

# 4. Save
dir_create(path_dir(snk$output$se))
saveRDS(se, snk$output$se)

# Save matrix for downstream GISTIC or visualization
write_matrix_tsv(
  matrix = cnv_data$matrix,
  row_ids = rownames(se),
  row_id_name = "hugo.id",
  output_path = snk$output$matrix
)

message("✅ CNV SE Curation Complete.")
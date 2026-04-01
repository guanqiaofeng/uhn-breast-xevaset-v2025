suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(dplyr)
})

snakemake@source("../helpers.R") 
snk <- parse_snakemake()

message("📊 Starting MultiAssayExperiment construction (Auto-generating Metadata)...")

# 1. Load Assays Dynamically
assay_paths <- unname(unlist(snk$input$assays))
assays <- lapply(assay_paths, readRDS)

assay_names <- vapply(
  seq_along(assays),
  function(idx) {
    # 1. Get the name from metadata or filename
    name <- S4Vectors::metadata(assays[[idx]])$datatype
    if (is.null(name) || !nzchar(name)) {
      name <- path_ext_remove(path_file(assay_paths[[idx]]))
    }
    
    # 2. CLEANUP LOGIC:
    # Remove "_SE" from the end
    name <- gsub("_SE$", "", name)
    
    # Remove everything in parentheses (like "(tumor-only WES)")
    name <- gsub("\\(.*\\)", "", name)
    
    # Trim any leftover whitespace
    name <- trimws(name)
    
    return(name)
  },
  character(1)
)
names(assays) <- make.unique(assay_names)
message("✅ Loaded assays: ", paste(names(assays), collapse = ", "))

# 2. Generate the Sample Map
# This links each SE column name to a 'primary' ID (the model name)
sample_maps <- lapply(names(assays), function(assay_name) {
  se <- assays[[assay_name]]
  
  # We assume the colnames of your curated SEs are already the clean MODEL_IDs
  # from your Python scripts (e.g., REF036)
  data.frame(
    assay = assay_name,
    primary = colnames(se),
    colname = colnames(se),
    stringsAsFactors = FALSE
  )
})
sample_map <- do.call(rbind, sample_maps)

# 3. Create Dynamic Master Metadata (colData)
# Gather all unique model IDs across all experiments
unique_models <- unique(sample_map$primary)

sample_df <- data.frame(
  sample_id = unique_models,
  patient_id = unique_models,
  tissue = "BRCA",
  tissue.name = "Breast Cancer",
  stringsAsFactors = FALSE
)
rownames(sample_df) <- sample_df$sample_id

# Add Data Availability Flags (has_rnaseq, has_cnv, has_mutation)
for (assay_type in names(assays)) {
  # Create a clean column name (e.g., "RNASeq_SE" -> "has_rnaseq")
  clean_name <- tolower(gsub("_SE|\\(.*\\)", "", assay_type))
  col_name <- paste0("has_", clean_name)
  
  # Identify models present in this specific assay
  present_models <- sample_map$primary[sample_map$assay == assay_type]
  
  # Assign TRUE/FALSE based on presence
  sample_df[[col_name]] <- sample_df$sample_id %in% present_models
}

message("📝 Generated metadata with availability flags for ", nrow(sample_df), " models.")

# Output Metadata Table for Records
metadata_output_path <- file.path(dirname(snk$output$mae), "model_metadata_final.tsv")
data.table::fwrite(sample_df, metadata_output_path, sep = "\t")
message("💾 Metadata table exported to: ", metadata_output_path)

# 4. Build and Save the MAE
mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = assays,
  colData = S4Vectors::DataFrame(sample_df),
  sampleMap = sample_map,
  metadata = list(
    dataset = "UHN_Breast_XevaSet",
    version = "2025_v1",
    genome = "GRCh38",
    curated_date = Sys.Date()
  )
)

message("💾 Saving MAE: ", ncol(mae), " models across ", length(assays), " assays.")

dir_create(path_dir(snk$output$mae))
saveRDS(mae, snk$output$mae)

message("🏁 MultiAssayExperiment construction finished successfully.")
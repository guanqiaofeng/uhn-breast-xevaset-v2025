setwd("/Users/guanqiaofeng/Documents/BHK/XevaDB/curation_pipeline/uhn_breast_xevaset_v2025")

# --- 1. SETUP & LIBRARIES ---
suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(MultiAssayExperiment)
  library(Xeva)
})

# --- 2. SWITCH: SNAKEMAKE VS. MANUAL ---
if (exists("snakemake")) {
  # Mode A: Running via Snakemake
  mae_path           <- snakemake@input[["multiAssayExperiment"]]
  model_info_path    <- snakemake@input[["modelInfo"]]
  drug_info_path     <- snakemake@input[["treatmentMetadata"]]
  exp_path           <- snakemake@input[["experiment"]]
  design_path        <- snakemake@input[["expDesign"]]
  batch_metrics_path <- snakemake@input[["batchMetrics"]]
  model_metrics_path <- snakemake@input[["modelMetrics"]]
  mod_map_path       <- snakemake@input[["modToBiobaseMap"]]
  
  output_xeva        <- snakemake@output[["xeva"]]
} else {
  # Mode B: Manual Testing (Update these paths to your local files)
  message("🛠️ Manual test mode detected. Setting local paths...")
  base_dir           <- "data/procdata/drug_screen/build"
  mae_path           <- "data/procdata/omics/UHN_Breast_MAE_v2025.rds"
  model_info_path    <- file.path(base_dir, "model.csv")
  drug_info_path     <- file.path(base_dir, "drug_annotated.csv")
  exp_path           <- file.path(base_dir, "experiment.csv")
  design_path        <- file.path(base_dir, "expDesign.csv")
  batch_metrics_path <- file.path(base_dir, "batch_metrics.csv")
  model_metrics_path <- file.path(base_dir, "model_metrics.csv")
  mod_map_path       <- file.path(base_dir, "modToBiobaseMap.csv")
  
  output_xeva        <- "data/results/UHN_Breast_XevaSet_v2025.rds"
}

# --- 3. DATA LOADING ---
message("📊 Loading curated data...")
model_df        <- data.table::fread(model_info_path) |> as.data.frame()
drug_df         <- data.table::fread(drug_info_path)  |> as.data.frame()
experiment_df   <- data.table::fread(exp_path)        |> as.data.frame()
exp_design      <- data.table::fread(design_path)     |> as.data.frame()
modToBiobaseMap <- data.table::fread(mod_map_path)    |> as.data.frame()
mae             <- readRDS(mae_path)
mol_list <- lapply(MultiAssayExperiment::experiments(mae), function(x) {
  if (inherits(x, "SummarizedExperiment")) {
    # Helper to convert SE to ExpressionSet if needed
    return(as(x, "ExpressionSet")) 
  }
  return(x)
})

## replace 'drug.id' to 'drug' in column name
if("drug.id" %in% colnames(experiment_df)) {
  setnames(experiment_df, "drug.id", "drug")
}

## "Intersection" check
#retained_drugs <- intersect(unique(drug_df$drug.id), unique(experiment_df$drug.id))
#experiment_df <- experiment_df[experiment_df$drug.id %in% retained_drugs, ]

## Re-format expDesign to match your working version
message("🏗️ Formatting expDesign as list-of-lists...")
exp_design_list <- lapply(1:nrow(exp_design), function(i) {
  treat_ids <- unlist(strsplit(as.character(exp_design$treatment[i]), ","))
  ctrl_ids  <- unlist(strsplit(as.character(exp_design$control[i]), ","))
  # Trim whitespace just in case there are spaces after the commas
  treat_ids <- trimws(treat_ids)
  ctrl_ids  <- trimws(ctrl_ids)
  
  list(
    batch.name = exp_design$batch.name[i],
    treatment  = treat_ids,
    control    = ctrl_ids
  )
})

names(exp_design_list) <- exp_design$batch.name

# Load Python "Gold Standard" metrics
py_batch_metrics <- fread(batch_metrics_path) |> as.data.frame()
py_model_metrics <- fread(model_metrics_path) |> as.data.frame()

# --- 4. CREATE XEVSET ---
message("🏗️ Initializing XevaSet...")

pdxe <- createXevaSet(
  name              = "UHN_Breast_v2025",
  model             = model_df,
  drug              = drug_df,
  experiment        = experiment_df,
  expDesign         = exp_design_list,
  modToBiobaseMap   = modToBiobaseMap,
  molecularProfiles = mol_list
)

# --- 5. INJECT PYTHON METRICS ---
message("💉 Injecting Python-calculated metrics...")

# A. Model-level injection (mRECIST, etc.)
# We must ensure the rownames match Xeva's expectation
# rownames(py_model_metrics) <- py_model_metrics$model.id
# slot(pdxe, "sensitivity")$model <- py_model_metrics

pdxe <- setResponse(pdxe, res.measure = "slope")
pdxe <- setResponse(pdxe, res.measure = "AUC")
pdxe <- setResponse(pdxe, res.measure = "angle")
pdxe <- setResponse(pdxe, res.measure = "abc")
pdxe <- setResponse(pdxe, res.measure = "mRECIST")
pdxe <- setResponse(pdxe, res.measure = "TGI")  
pdxe <- setResponse(pdxe, res.measure = "lmm") 

# B. Batch-level injection (ABC, Angle, etc.)
# rownames(py_batch_metrics) <- py_batch_metrics$batch.name
# slot(pdxe, "sensitivity")$batch <- py_batch_metrics

# --- 6. VALIDATION & SAVE ---
final_counts <- table(slot(pdxe, "sensitivity")$batch$mRECIST)
message("\n✅ Final XevaSet mRECIST Summary:")
print(final_counts)

# Create directory if it doesn't exist and save
results_dir <- path_dir(output_xeva)
dir_create(results_dir)

# Export the sensitivity tables actually stored in the XevaSet as CSV
message("💾 Exporting model- and batch-level sensitivity tables...")
model_sensitivity <- sensitivity(pdxe, type = "model")
batch_sensitivity <- sensitivity(pdxe, type = "batch")
fwrite(model_sensitivity, file.path(results_dir, "UHN_Breast_XevaSet_v2025_model_sensitivity.csv"))
fwrite(batch_sensitivity, file.path(results_dir, "UHN_Breast_XevaSet_v2025_batch_sensitivity.csv"))

saveRDS(pdxe, output_xeva)
message("\n🚀 XevaSet successfully saved to: ", output_xeva)

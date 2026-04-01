import os
import logging
import pandas as pd
from pathlib import Path

# 1. ENVIRONMENT & PATHS
PROJECT_ROOT = Path(os.getenv("PIXI_PROJECT_ROOT", Path(__file__).parents[4]))
RAW_DIR = PROJECT_ROOT / "data" / "rawdata"
PROC_DIR = PROJECT_ROOT / "data" / "procdata"
LOG_DIR = PROJECT_ROOT / "logs" / "omics"

LOG_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "format_cnv.log"),
        logging.StreamHandler()
    ]
)

def run_cnv_curation():
    logging.info("🚀 Starting CNV Curation (Zero Filtering / Raw Preservation)")

    # Path Definitions
    cnv_path = RAW_DIR / "omics" / "20251121_cnv_matrix.tsv"
    model_meta_path = RAW_DIR / "omics" / "omics_raw_metadata_v2025.xlsx"
    info_path = RAW_DIR / "additional" / "geneInfo.tab"
    
    output_matrix = PROC_DIR / "omics" / "cnv_hugo.tsv"
    output_features = PROC_DIR / "omics" / "feature_metadata_cnv.tsv"

    # 2. LOAD DATA
    df = pd.read_csv(cnv_path, sep='\t')
    
    # Identify sample columns
    sample_cols = [c for c in df.columns if c != 'gene']
    
    # Ensure numeric types while preserving NaNs
    df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors='coerce')

    # Load Annotation for metadata consistency
    info = pd.read_csv(
        info_path, 
        sep='\t', 
        skiprows=1, 
        names=['ensembl.id', 'hugo.id', 'description']
    )

    model_meta = pd.read_excel(
        model_meta_path,
        sheet_name="WES",
        usecols=["in_XevaSet_v2025", "MODEL_ID", "SAMPLE_ID_ORIGINAL"]
    )

    # 3. CLEAN MODEL NAMES & FILTER
    logging.info("🧹 Filtering and renaming columns based on metadata")

    # Filter for samples intended for the 2025 XevaSet
    valid_meta = model_meta[model_meta["in_XevaSet_v2025"].str.upper() == "Y"].copy()
    logging.info(f"📋 Metadata check: Found {len(valid_meta)} samples marked as 'Y' in Excel.")

    # Create a mapping dictionary: {Original_ID: New_Model_ID}
    rename_map = dict(zip(valid_meta["SAMPLE_ID_ORIGINAL"], valid_meta["MODEL_ID"]))

    # Identify which columns to keep and which are missing
    original_cols = set(df.columns)
    requested_cols = set(rename_map.keys())
    
    found_cols = original_cols.intersection(requested_cols)
    missing_cols = requested_cols - original_cols

    if missing_cols:
        logging.warning(f"⚠️ {len(missing_cols)} samples in metadata were NOT found in the CSV headers: {list(missing_cols)[:5]}...")
    
    logging.info(f"🧬 Matching: {len(found_cols)} samples successfully matched between CSV and Metadata.")

    # Subset the expression matrix
    keep_cols = ['gene'] + list(found_cols)
    df_filtered = df[keep_cols].copy()

    # Rename the columns to the clean MODEL_ID
    df_filtered = df_filtered.rename(columns=rename_map)
    logging.info(f"✨ Column renaming complete. New matrix width: {df_filtered.shape[1]} (including gene column).")

    # 4. IDENTIFIER ALIGNMENT (No Filtering)
    # Rename 'gene' to 'hugo.id' to match your RNAseq structure
    df_filtered = df_filtered.rename(columns={'gene': 'hugo.id'})
    
    # If there are duplicate Hugo IDs, we keep them all or drop exact duplicates
    # Genomic matrices occasionally have duplicate rows; let's keep it safe.
    df_filtered = df_filtered.drop_duplicates().copy()

    # 5. FINAL FORMATTING
    # Redefine sample_cols based on the NEW column names (excluding hugo.id)
    new_sample_cols = [c for c in df_filtered.columns if c != 'hugo.id']
    final_matrix = df_filtered[['hugo.id'] + new_sample_cols]

    # Prepare Feature Metadata
    # Merge on 'hugo.id' to get Ensembl IDs and Descriptions for all genes present
    feature_metadata = final_matrix[['hugo.id']].merge(info, on='hugo.id', how='left')

    # 6. EXPORT
    PROC_DIR.joinpath("omics").mkdir(parents=True, exist_ok=True)
    
    # Save with 'NA' representation to match your preferred raw format
    final_matrix.to_csv(output_matrix, index=False, na_rep='NA', sep='\t')
    feature_metadata.to_csv(output_features, index=False, na_rep='NA', sep='\t')

    logging.info(f"✅ Success. CNV Matrix preserved with {len(final_matrix)} genes.")
    logging.info(f"Matrix saved to: {output_matrix}")

if __name__ == "__main__":
    try:
        run_cnv_curation()
    except Exception as e:
        logging.error(f"CNV Workflow failed: {e}", exc_info=True)
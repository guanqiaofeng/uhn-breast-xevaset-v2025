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
        logging.FileHandler(LOG_DIR / "format_rnaseq.log"),
        logging.StreamHandler()
    ]
)

def run_curation():
    logging.info("🚀 Starting RNAseq Curation Workflow")

    # Path Definitions
    rna_path = RAW_DIR / "omics" / "20251007_gene_tpm_normalized_matrix.csv"
    model_meta_path = RAW_DIR / "omics" / "omics_raw_metadata_v2025.xlsx"
    info_path = RAW_DIR / "additional" / "geneInfo.tab"
    
    # Output Definitions
    output_matrix = PROC_DIR / "omics" / "rna_tpm_hugo.tsv"
    output_features = PROC_DIR / "omics" / "feature_metadata_rna.tsv"

    # 2. LOAD DATA
    df = pd.read_csv(rna_path) 
    if 'gene' not in df.columns:
        # Assumes first column is the Ensembl ID if not named 'gene'
        df = df.rename(columns={df.columns[0]: 'gene'})

    info = pd.read_csv(
        info_path, 
        sep='\t', 
        skiprows=1, 
        names=['ensembl.id', 'hugo.id', 'description']
    )

    model_meta = pd.read_excel(
        model_meta_path,
        sheet_name="RNAseq",
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

    # 4. MAPPING & CLEANING
    mapping_dic = info.set_index('ensembl.id')['hugo.id'].to_dict()
    df_filtered['hugo.id'] = df_filtered['gene'].map(mapping_dic)

    unmapped = df_filtered[df_filtered['hugo.id'].isna()]
    if not unmapped.empty:
        logging.warning(f"Found {len(unmapped)} Ensembl IDs with no Hugo mapping.")
        unmapped[['gene']].to_csv(LOG_DIR / "unmapped_genes.txt", index=False)

    # 5.1. DEDUPLICATION
    numeric_cols = df_filtered.select_dtypes(include=['number']).columns
    df_filtered['expression_sum'] = df_filtered[numeric_cols].sum(axis=1)

    df_clean = df_filtered.dropna(subset=['hugo.id'])
    df_deduplicated = (
        df_clean.sort_values("expression_sum", ascending=False)
        .drop_duplicates(subset="hugo.id", keep="first")
    ).copy()

    # 5.2 FILTER ZERO-EXPRESSION GENES
    sample_cols = [c for c in df_deduplicated.columns if c not in ['gene', 'hugo.id', 'expression_sum']]
    
    initial_count = len(df_deduplicated)
    df_deduplicated = df_deduplicated[df_deduplicated[sample_cols].sum(axis=1) > 0].copy()
    
    removed_count = initial_count - len(df_deduplicated)
    logging.info(f"Filtered out {removed_count} genes with 0 expression across all samples.")

    # 6. FINAL FORMATTING
    # Prepare the Expression Matrix (Hugo ID as first column)
    matrix_cols = ['hugo.id'] + sample_cols
    final_matrix = df_deduplicated[matrix_cols]

    # Prepare Feature Metadata
    feature_metadata = df_deduplicated[['gene', 'hugo.id']].rename(columns={'gene': 'ensembl.id'})
    feature_metadata = feature_metadata.merge(info, on=['ensembl.id', 'hugo.id'], how='left')

    # 7. EXPORT
    output_matrix.parent.mkdir(parents=True, exist_ok=True)
    final_matrix.to_csv(output_matrix, index=False, sep='\t')
    feature_metadata.to_csv(output_features, index=False, sep='\t')

    logging.info(f"✅ Success. Matrix shape: {final_matrix.shape}")
    logging.info(f"Matrix saved to: {output_matrix}")

if __name__ == "__main__":
    try:
        run_curation()
    except Exception as e:
        logging.error(f"Workflow failed: {e}", exc_info=True)
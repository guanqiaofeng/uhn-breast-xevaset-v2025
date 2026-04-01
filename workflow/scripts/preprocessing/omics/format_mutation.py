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
        logging.FileHandler(LOG_DIR / "format_mutation.log"),
        logging.StreamHandler()
    ]
)

def run_mutation_curation():
    logging.info("🚀 Starting Mutation Curation (Output: TSV for comma-safety)")

    # Path Definitions
    maf_path = RAW_DIR / "omics" / "cohort_20260217_gnomad0002_vaf04_alt5_updaterescue_allgermline.maf.oncokb.txt"
    model_meta_path = RAW_DIR / "omics" / "omics_raw_metadata_v2025.xlsx"
    info_path = RAW_DIR / "additional" / "geneInfo.tab"
    
    # Using .tsv for the matrix to avoid CSV delimiter conflicts
    output_matrix = PROC_DIR / "omics" / "mutation_hugo.tsv"
    output_features = PROC_DIR / "omics" / "feature_metadata_mutation.tsv"

    # 2. LOAD MAF DATA
    logging.info(f"Reading MAF: {maf_path.name}")
    df = pd.read_csv(
        maf_path, 
        sep='\t', 
        low_memory=False, 
        usecols=['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Variant_Classification']
    )

    model_meta = pd.read_excel(
        model_meta_path,
        sheet_name="WES",
        usecols=["in_XevaSet_v2025", "MODEL_ID", "SAMPLE_ID_ORIGINAL"]
    )

    # 3. PIVOT TO WIDE FORMAT
    # Comma-joining multiple mutations per gene/sample
    logging.info("Pivoting MAF to Wide-Format...")
    mutation_matrix = df.pivot_table(
        index='Hugo_Symbol',
        columns='Tumor_Sample_Barcode',
        values='Variant_Classification',
        aggfunc=lambda x: ",".join(sorted(set(x.astype(str))))
    )
    mutation_matrix.index.name = 'hugo.id'

    # 4. CLEAN MODEL NAMES & FILTER
    logging.info("🧹 Filtering and renaming columns based on metadata")

    valid_meta = model_meta[model_meta["in_XevaSet_v2025"].str.upper() == "Y"].copy()
    rename_map = dict(zip(valid_meta["SAMPLE_ID_ORIGINAL"], valid_meta["MODEL_ID"]))

    # Subset to samples found in metadata
    found_cols = [c for c in mutation_matrix.columns if c in rename_map]
    mutation_matrix_filtered = mutation_matrix[found_cols].copy()

    # Rename columns
    mutation_matrix_filtered = mutation_matrix_filtered.rename(columns=rename_map)

    # CHECK FOR DUPLICATES: If multiple barcodes mapped to the same MODEL_ID
    if mutation_matrix_filtered.columns.duplicated().any():
        logging.info("🧬 Multiple barcodes found for the same MODEL_ID. Merging mutations...")
        
        # Group by column name and join unique mutation strings with commas
        mutation_matrix_filtered = mutation_matrix_filtered.T.groupby(level=0).agg(
            lambda x: ",".join(sorted(set(",".join(x.astype(str)).split(","))))
        ).T
        
        # Clean up strings (remove 'nan' or empty values if they were joined)
        mutation_matrix_filtered = mutation_matrix_filtered.replace(r'^nan,|,nan$|^nan$', '', regex=True)

    # Now sorting will work because labels are unique
    mutation_matrix_filtered = mutation_matrix_filtered.reindex(
        sorted(mutation_matrix_filtered.columns), axis=1
    )

    logging.info(f"✨ Final mutation matrix width: {mutation_matrix_filtered.shape[1]} unique models.")
    mutation_matrix = mutation_matrix_filtered

    # 5. FEATURE METADATA ALIGNMENT
    logging.info("Mapping Hugo IDs to Ensembl for metadata...")
    info = pd.read_csv(
        info_path, 
        sep='\t', 
        skiprows=1, 
        names=['ensembl.id', 'hugo.id', 'description']
    )
    feature_list = mutation_matrix.index.to_frame(index=False)
    feature_metadata = feature_list.merge(info, on='hugo.id', how='left')

    # 6. EXPORT
    PROC_DIR.joinpath("omics").mkdir(parents=True, exist_ok=True)
    
    # Export Matrix as TSV (sep='\t')
    mutation_matrix.to_csv(output_matrix, sep='\t', index=True)
    # Export Metadata as CSV
    feature_metadata.to_csv(output_features, index=False, sep='\t')

    logging.info(f"✅ Success. Mutation Matrix shape: {mutation_matrix.shape}")
    logging.info(f"TSV Matrix saved to: {output_matrix}")

if __name__ == "__main__":
    try:
        run_mutation_curation()
    except Exception as e:
        logging.error(f"Mutation Workflow failed: {e}", exc_info=True)
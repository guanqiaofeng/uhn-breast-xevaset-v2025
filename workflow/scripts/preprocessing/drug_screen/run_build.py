import sys
import logging
import pandas as pd
from pathlib import Path

# SETUP PROJECT PATHS
PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.append(str(PROJECT_ROOT))

from workflow.utilities.build_function import (
    build_model, build_drug, build_experiment, 
    build_expDesign, build_modToBiobaseMap_pdata
)

def main():
    # 1. DIRECTORY CONFIG
    DATA_DIR = PROJECT_ROOT / "data"
    PROC_DIR = DATA_DIR / "procdata" / "drug_screen"
    OMICS_DIR = DATA_DIR / "procdata" / "omics"
    OUT_DIR = PROC_DIR / "build"
    LOG_DIR = PROJECT_ROOT / "logs" / "drug_screen"
    
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    # 2. LOGGING CONFIGURATION
    log_file = LOG_DIR / "build.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ],
        force=True
    )

    logging.info("🏗️  Starting XevaSet Build Phase...")
    logging.info(f"Log file: {log_file}")

    # 3. LOAD DATA FROM UPDATE STEP
    input_file = PROC_DIR / "update" / "update_model_drug_batch_data.tsv"
    if not input_file.exists():
        logging.error(f"❌ Critical Error: Input file not found at {input_file}")
        return
    
    df = pd.read_csv(input_file, sep="\t", low_memory=False)

    # 4. GENERATE CORE COMPONENTS
    logging.info("Step 3.1-3.4: Generating Model, Drug, Experiment, and Design tables...")
    model_df = build_model(df)
    drug_df = build_drug(df)
    exp_df = build_experiment(df)
    design_df = build_expDesign(df)

    # 5. GENERATE OMICS MAPPING
    logging.info("Step 3.5: Mapping Model IDs to Omics Sample IDs...")
    map_df, pdata_dict = build_modToBiobaseMap_pdata(df, OMICS_DIR)

    # 6. SAVE OUTPUTS (Standard CSV for Xeva Compatibility)
    logging.info(f"Saving finalized CSV files to {OUT_DIR}...")
    model_df.to_csv(OUT_DIR / "model.csv", index=False)
    drug_df.to_csv(OUT_DIR / "drug.csv", index=False)
    exp_df.to_csv(OUT_DIR / "experiment.csv", index=False)
    design_df.to_csv(OUT_DIR / "expDesign.csv", index=False)
    
    if not map_df.empty:
        map_df.to_csv(OUT_DIR / "modToBiobaseMap.csv", index=False)
        for dtype, pdata in pdata_dict.items():
            pdata.to_csv(OUT_DIR / f"{dtype}_pdata.csv", index=False)

    # 7. FINAL INTEGRITY SUMMARY
    logging.info("=" * 45)
    logging.info("       BUILD STEP: FINAL SUMMARY")
    logging.info("=" * 45)
    logging.info(f"Models Ready:        {len(model_df)}")
    logging.info(f"Drugs Registered:    {len(drug_df)}")
    logging.info(f"Total Observations:  {len(exp_df)}")
    logging.info(f"Treatment Batches:   {len(design_df)}")
    
    if not map_df.empty:
        logging.info("-" * 45)
        logging.info("MULTI-OMIC COVERAGE SUMMARY:")
        
        # Pivot to see which models have which data types
        map_df = map_df[['biobase.id', 'mDataType']].drop_duplicates()
        coverage = map_df.groupby('biobase.id')['mDataType'].apply(list).reset_index()
        
        total_models = len(model_df)
        has_rna = coverage['mDataType'].apply(lambda x: 'RNASeq' in x).sum()
        has_mut = coverage['mDataType'].apply(lambda x: 'mutation' in x).sum()
        has_cnv = coverage['mDataType'].apply(lambda x: 'CNV' in x).sum()
        
        all_three = coverage['mDataType'].apply(lambda x: all(t in x for t in ['RNASeq', 'mutation', 'CNV'])).sum()

        logging.info(f"  - Models with RNASeq:   {has_rna} / {total_models}")
        logging.info(f"  - Models with Mutation: {has_mut} / {total_models}")
        logging.info(f"  - Models with CNV:      {has_cnv} / {total_models}")
        logging.info(f"  - Models with ALL 3:    {all_three} (Gold Standard)")

if __name__ == "__main__":
    main()
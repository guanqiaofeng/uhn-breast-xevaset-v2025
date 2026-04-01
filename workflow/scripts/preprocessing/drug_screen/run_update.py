import sys
import logging
import argparse
import pandas as pd
from pathlib import Path

# SETUP PROJECT PATHS
PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.append(str(PROJECT_ROOT))

from workflow.utilities.update_function import (
    update_model_data, map_model, update_drug_data, update_batch_data
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--modeltype', default='both', choices=['resistant', 'parental', 'both'])
    args = parser.parse_args()

    # PATHS
    PROC_DIR = PROJECT_ROOT / "data" / "procdata" / "drug_screen"
    ADD_DIR = PROJECT_ROOT / "data" / "rawdata" / "additional"
    OUT_DIR = PROC_DIR / "update"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # LOGGING
    LOG_DIR = PROJECT_ROOT / "logs" / "drug_screen"
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_file = LOG_DIR / "update.log"

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),  # Saves to logs/drug_screen/update.log
            logging.StreamHandler(sys.stdout) # Prints to your MacBook terminal
        ],
        force=True # Ensures configuration is applied even if logging was started elsewhere
    )
    
    logging.info(f"🚀 Updating Drug Screen Data ({args.modeltype})...")
    logging.info(f"Log file: {log_file}")

    # 1. LOAD PROCESSED DATA
    input_tsv = PROC_DIR / "process" / "processed_data.tsv"
    df = pd.read_csv(input_tsv, sep="\t")

    # 2. RUN UPDATE STEPS (DataFrame Pipeline)
    df_updated = update_model_data(df)
    
    df_mapped = map_model(
        df_updated, 
        str(ADD_DIR / "all_modelid_sampleid_modeltype_mapping.xlsx"), 
        args.modeltype
    )
    
    df_drug, df_filtered = update_drug_data(
        df_mapped, 
        str(ADD_DIR / "all_drug.xlsx")
    )
    
    final_df = update_batch_data(df_drug)

    # 3. SAVE
    final_output = OUT_DIR / "update_model_drug_batch_data.tsv"
    filtered_output = OUT_DIR / "drug_filtered_out.tsv"
    
    final_df.to_csv(final_output, sep="\t", index=False)
    df_filtered.to_csv(filtered_output, sep="\t", index=False)
    
    # --- DATA INTEGRITY SUMMARY ---
    logging.info("=" * 45)
    logging.info("       UPDATE STEP: FINAL SUMMARY")
    logging.info("=" * 45)
    logging.info(f"Input Rows:          {len(df)}")
    logging.info(f"Final Rows:          {len(final_df)}")
    logging.info(f"Filtered Out:        {len(df_filtered)} rows (Check drug_filtered_out.tsv)")
    logging.info("-" * 45)
    logging.info(f"Unique Models:       {final_df['modelID'].nunique()}")
    logging.info(f"Unique Samples:      {final_df['sampleID'].nunique()}")
    logging.info(f"Unique Batches:      {final_df['batch'].nunique()}")
    logging.info(f"Unique Mice:         {final_df['model.id'].nunique()}")
    logging.info("-" * 45)
    
    # Check for any "UNKNOWN" tags from the regex step
    unknowns = final_df['modelID'].str.contains("UNKNOWN").sum()
    if unknowns > 0:
        logging.warning(f"⚠️ ALERT: {unknowns} rows contain 'UNKNOWN' IDs. Regex may have failed.")
    else:
        logging.info("✅ Regex Extraction: All IDs parsed successfully.")
        
    logging.info(f"Results saved to: {OUT_DIR}")
    logging.info("=" * 45)

if __name__ == "__main__":
    main()
import os
import pandas as pd
import sys
import logging
import yaml
from pathlib import Path

# 1. SETUP PROJECT PATHS
# Script is at: workflow/scripts/preprocessing/drug_screen/run_process.py
# Go up 4 levels to reach the root: uhn_breast_xevaset_v2025/
PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.append(str(PROJECT_ROOT))

# Updated imports from your new workflow.utilities structure
from workflow.utilities.process_function import (
    process_all_files, 
    save_processed_data, 
    save_qc_report
)

def setup_logging(log_dir: Path):
    """Initializes logging to both a permanent log file and the terminal."""
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "process.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ],
        force=True
    )
    return log_file

def main():
    # 2. LOAD PIPELINE CONFIGURATION
    config_path = PROJECT_ROOT / "config" / "pipeline.yaml"
    
    if not config_path.exists():
        print(f"❌ Error: Configuration file not found at {config_path}")
        sys.exit(1)
        
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    
    # Extract parameters from YAML
    # Defaulting to 14 if the key is missing
    min_days = config.get("drug_screen", {}).get("processing", {}).get("min_treatment_days", 14)

    # 3. DIRECTORY CONFIGURATION
    RAW_DIR = PROJECT_ROOT / "data" / "rawdata"
    PROC_DIR = PROJECT_ROOT / "data" / "procdata"
    LOG_DIR = PROJECT_ROOT / "logs" / "drug_screen"
    
    # Input: Result from the Scan step
    metadata_file = PROC_DIR / "drug_screen" / "scan" / "all_file_scan.json"
    drug_screen_dir = RAW_DIR / "drug_screen"
    
    # Output: Process directory
    output_dir = PROC_DIR / "drug_screen" / "process"
    output_dir.mkdir(parents=True, exist_ok=True)

    # 4. INITIALIZE LOGGING
    log_path = setup_logging(LOG_DIR)
    logging.info("🚀 Starting Drug Screen Data Processing...")
    logging.info(f"Threshold: Including replicates with ≥ {min_days} days of data.")

    # 5. EXECUTION
    # Convert Paths to strings for compatibility with the utility functions
    all_data, qc_issues = process_all_files(
        str(metadata_file), 
        str(drug_screen_dir), 
        min_days=min_days
    )

    # Define output file paths
    processed_tsv = output_dir / "processed_data.tsv"
    qc_report_json = output_dir / "qc_report.json"

    # 6. SAVE OUTPUTS & REPORTING
    save_processed_data(all_data, str(processed_tsv))
    save_qc_report(qc_issues, str(qc_report_json))
    
    logging.info("-" * 30)
    logging.info("       PROCESS COMPLETE")
    logging.info("-" * 30)
    logging.info(f"Master TSV: {processed_tsv}")
    logging.info(f"QC Report:  {qc_report_json}")
    logging.info(f"All logs:   {log_path}")
    logging.info("-" * 30)

    # --- 7. DATA INTEGRITY SUMMARY ---
    if all_data:
        # We re-load or use the concatenated data for a final check
        master_df = pd.concat(all_data, ignore_index=True)
        
        logging.info("=" * 35)
        logging.info("       DATA INTEGRITY REPORT")
        logging.info("=" * 35)
        logging.info(f"Total Observations: {len(master_df)}")
        logging.info(f"Unique Samples:     {master_df['sampleID'].nunique()}")
        logging.info(f"Drugs Found:        {', '.join(master_df['drug'].unique())}")
        logging.info(f"Time Range:         {master_df['time'].min()} to {master_df['time'].max()} days")
        
        missing_v = master_df['volume'].isna().sum()
        if missing_v == 0:
            logging.info("✅ Volume Imputation: SUCCESS (No missing values)")
        else:
            logging.warning(f"⚠️ Volume Imputation: INCOMPLETE ({missing_v} missing values)")
        logging.info("=" * 35)
    
    logging.info(f"All process logs saved to: {log_path}")

if __name__ == "__main__":
    main()
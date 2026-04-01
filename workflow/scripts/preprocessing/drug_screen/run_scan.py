import os
import sys
import logging
import pandas as pd
from pathlib import Path

# 1. SETUP PROJECT PATHS
# Adjust parents[4] to reach the root 'uhn_breast_xevaset_v2025'
PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.append(str(PROJECT_ROOT))

from workflow.utilities.scan_function import collect_metadata, save_metadata

def setup_logging(log_dir: Path):
    """Sets up logging to both the console and a file in logs/."""
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "scan.log"
    
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
    # 2. CONFIGURATION
    RAW_DIR = PROJECT_ROOT / "data" / "rawdata"
    PROC_DIR = PROJECT_ROOT / "data" / "procdata"
    LOG_DIR = PROJECT_ROOT / "logs" / "drug_screen"
    
    DRUG_SCREEN_DIR = RAW_DIR / "drug_screen"
    ADDITIONAL_DIR = RAW_DIR / "additional"
    OUTPUT_DIR = PROC_DIR / "drug_screen" / "scan"
    
    # Input files
    remove_list_path = ADDITIONAL_DIR / "models_to_remove.xlsx"
    
    # 3. INITIALIZE LOGGING
    log_path = setup_logging(LOG_DIR)
    logging.info(f"🚀 Starting Drug Screen Scan. Project Root: {PROJECT_ROOT}")
    
    # 4. EXECUTION
    logging.info(f"Scanning directory: {DRUG_SCREEN_DIR}")
    metadata = collect_metadata(str(DRUG_SCREEN_DIR), str(remove_list_path))

    # Ensure the output directory exists
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Define output paths
    output_json = OUTPUT_DIR / "all_file_scan.json"
    output_excel = OUTPUT_DIR / "all_file_scan.xlsx"
    output_tsv = OUTPUT_DIR / "all_file_scan.tsv"
    
    # Save outputs
    save_metadata(metadata, str(output_json), str(output_excel), str(output_tsv))

    # 5. GENERATE LOGGED SUMMARY (For PPT and Paper prep)
    df = pd.DataFrame(metadata)
    
    logging.info("-" * 30)
    logging.info("       SCAN SUMMARY")
    logging.info("-" * 30)
    logging.info(f"Total Files Found:   {len(df)}")
    
    # Model Type Breakdown
    if not df.empty:
        type_counts = df['model_type'].value_counts()
        for mtype, count in type_counts.items():
            logging.info(f"  - {mtype.capitalize()}: {count}")

        # Tag Analysis (Duplicates/Exclusions)
        duplicates = df['tag'].str.contains('duplicated_file', na=False).sum()
        excluded = df['tag'].str.contains('in_remove_list', na=False).sum()
        mouse = df['tag'].str.contains('mouse', na=False).sum()
        
        logging.info(f"Duplicates Flagged:  {duplicates}")
        logging.info(f"Excluded Models:     {excluded}")
        logging.info(f"Mouse Files Found:   {mouse}")
        
        # This is the "Gold Number" for your 81 models
        clean_count = len(df) - duplicates - excluded
        logging.info(f"✅ Clean Models Ready: {clean_count}")
    else:
        logging.warning("No files were found in the scan directory!")
        
    logging.info("-" * 30)
    logging.info(f"All logs saved to: {log_path}")

if __name__ == "__main__":
    main()
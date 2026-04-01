import os
import re
import json
import hashlib
import logging
import pandas as pd
from datetime import datetime, timezone
from typing import List, Dict

# Standard logging configuration for traceability
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def calculate_checksum(file_path: str) -> str:
    """
    Calculate the SHA256 checksum of a file to ensure data integrity.
    """
    sha256_hash = hashlib.sha256()
    try:
        with open(file_path, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
    except Exception as e:
        logging.error(f"Failed to read file {file_path}: {e}")
        return ""
    
    return sha256_hash.hexdigest()

def collect_metadata(data_dir: str, remove_file: str) -> List[Dict[str, str]]:
    """
    Traverse a directory and collect metadata of all .xlsx files, 
    tagging duplicates and excluded models.
    """
    # 1. Handle the Exclusion List safely
    remove_list = []
    if os.path.exists(remove_file):
        try:
            remove_df = pd.read_excel(remove_file)
            # Standardize column names to lowercase for robust matching
            remove_df.columns = [c.lower().replace(" ", ".") for c in remove_df.columns]
            
            if 'leave.out' in remove_df.columns and 'file.name' in remove_df.columns:
                # Filter for 'yes' and strip any accidental whitespace from filenames
                remove_list = remove_df[
                    remove_df['leave.out'].astype(str).str.lower() == 'yes'
                ]['file.name'].str.strip().tolist()
            else:
                logging.warning(f"Exclusion file {remove_file} missing required columns 'leave.out' or 'file.name'.")
        except Exception as e:
            logging.error(f"Could not parse exclusion file {remove_file}: {e}")
    else:
        logging.info(f"No exclusion file found at {remove_file}. All files will be scanned.")

    seen_cores = set() 
    metadata_list = []

    # 2. Traverse the directory
    for root, _, files in os.walk(data_dir):
        folder_name = os.path.basename(root)
        for file in files:
            # Skip temp Excel files and non-xlsx files
            if file.endswith(".xlsx") and not file.startswith("~$"): 
                file_path = os.path.join(root, file)
                
                # Calculate relative path for the Snakemake manifest
                relative_file_path = os.path.relpath(file_path, data_dir)
                
                # Regex to find the "Core" name (stripping 'Copy of' etc.)
                file_core = re.sub(r"(?i)^x?copy\s+of\s+", "", file).strip()
                
                # Determine Model Type (Parental vs Resistant)
                model_type = "resistant" if "res" in file_core.lower() else "parental"
                
                # 3. Tagging Logic
                tags = []
                if file_core.lower() in seen_cores:
                    tags.append("duplicated_file")
                if "mouse" in file_core.lower():
                    tags.append("mouse")
                if any(file.lower() == r.lower() for r in remove_list):
                    tags.append("in_remove_list")
                
                tag_str = ", ".join(tags)
                seen_cores.add(file_core.lower())
                
                # 4. Extract File Statistics
                try:
                    m_time = os.path.getmtime(file_path)
                    # Modern timezone-aware datetime
                    mod_date = datetime.fromtimestamp(m_time, tz=timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
                    
                    checksum = calculate_checksum(file_path)
                    file_size_kb = round(os.path.getsize(file_path) / 1024, 2)

                    metadata_list.append({
                        "folder_name": folder_name,
                        "file_name": file,
                        "file_name_core": file_core,
                        "file_path": relative_file_path,
                        "model_type": model_type,
                        "file_size_kb": file_size_kb,
                        "last_modified_date": mod_date,
                        "checksum": checksum,
                        "tag": tag_str
                    })
                except Exception as e:
                    logging.error(f"Error processing metadata for {file}: {e}")
    
    return metadata_list

def save_metadata(metadata_list: List[Dict[str, str]], output_json: str, output_excel: str, output_tsv: str) -> None:
    """
    Save curated metadata to JSON for pipeline use and Excel for human review.
    """
    try:
        # Save JSON (The primary input for 'run_process.py')
        with open(output_json, "w") as jf:
            json.dump(metadata_list, jf, indent=4)
        logging.info(f"Successfully saved manifest to JSON: {output_json}")
                
        # Save Excel (For you to check before the PPT)
        df = pd.DataFrame(metadata_list)
        df.to_excel(output_excel, index=False)
        logging.info(f"Successfully saved manifest to Excel: {output_excel}")

        # save tsv
        df.to_csv(output_tsv, sep='\t', index=False)
        logging.info(f"Metadata TSV saved: {output_tsv}")

    except Exception as e:
        logging.error(f"Failed to save metadata files: {e}")
        
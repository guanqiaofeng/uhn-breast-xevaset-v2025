import os
import re
import json
import logging
import pandas as pd
from typing import List, Dict, Tuple, Any

# Configure logging for professional traceability
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def process_file(file_entry: Dict[str, Any], base_dir: str, min_days: int) -> List[pd.DataFrame]:
    """
    Parses a single Excel file to extract PDX drug response data.
    
    Args:
        file_entry: Dictionary containing file metadata (path, checksum, etc.)
        base_dir: Root directory for the drug screen data.
        min_days: Minimum required treatment duration (e.g., 14) to include a replicate.
        
    Returns:
        List of DataFrames, one for each valid replicate found in the sheet.
    """
    processed_dfs = []
    file_path = os.path.join(base_dir, file_entry["file_path"])
    model_type = file_entry.get("model_type", "unknown")

    try:
        # Read raw Excel (no header) to allow coordinate-based extraction for Cescon Lab Template
        df = pd.read_excel(file_path, header=None)
        
        # --- Metadata Extraction ---
        # IDs are in Row 0, Col 0 and Row 0, Col 1
        sampleID = str(df.iloc[0, 0]).replace("\"", "")
        sampleID = re.sub(r"\s+", " ", sampleID).upper()
        
        implantation_date = df.iloc[1, 0]
        drug_start_date = df.iloc[2, 0]
        mouse_strain = re.sub(r"\s+", " ", str(df.iloc[0, 1]))
        implantation_location = df.iloc[1, 1]
        implantation_tissue_source = df.iloc[2, 1]

        # --- Timeframe Extraction ---
        date_list = []
        rel_date_list = []
        # Relative dates start at Column 4 (E) and repeat every 3 columns
        columns_to_check = range(4, df.shape[1], 3)
        
        for col in columns_to_check:
            date_val = df.iloc[1, col]
            rel_day = df.iloc[2, col]
            if pd.isna(date_val):
                break
            date_list.append(date_val.date() if hasattr(date_val, 'date') else date_val)
            rel_date_list.append(rel_day)
            
        if rel_date_list:
            # Normalize: Ensure the first recorded day is the 0-point for the timeline
            day_zero = rel_date_list[0]
            rel_date_list = [d - day_zero for d in rel_date_list]
        
        # --- Treatment Block Processing ---
        # Specific row indices where treatment groups start in the template
        row_indices = [3, 13, 23, 33, 43, 53]

        for r in row_indices:
            # Validate block (Check if a measurement exists in the first replicate row)
            if r + 1 >= df.shape[0] or pd.isna(df.iloc[r + 1, 2]):
                break

            # Standardize Drug and Dose names
            drug_name = re.sub(r"\s+", "", str(df.iloc[r + 1, 0])).upper()
            dose_val = re.sub(r"\s+", "", str(df.iloc[r + 6, 1])).lower()

            if drug_name in ["H20", "H2O"]:
                drug_name = "H2O"
                dose_val = "-"

            # Process up to 5 replicates (mice) per treatment block
            for rep in range(5):
                curr_row = r + rep + 1
                if curr_row >= df.shape[0] or pd.isna(df.iloc[curr_row, 2]):
                    break

                s_vals, l_vals, v_vals = [], [], []
                for p in range(len(rel_date_list)):
                    s_vals.append(df.iloc[curr_row, 2 + 3 * p])     # Width
                    l_vals.append(df.iloc[curr_row, 2 + 3 * p + 1]) # Length
                    v_vals.append(df.iloc[curr_row, 2 + 3 * p + 2]) # Volume

                # --- Quality Filter: Minimum Duration ---
                # Only include replicates that reached the configured 'min_days' threshold
                if rel_date_list and rel_date_list[-1] >= min_days:
                    rep_df = pd.DataFrame({
                        "file.name": file_entry["file_path"],
                        "file.checksum": file_entry["checksum"],
                        "model_type": model_type,
                        "sampleID": sampleID,
                        "drug": drug_name,
                        "dose": dose_val,
                        "date": date_list,
                        "time": rel_date_list,
                        "width": s_vals,
                        "length": l_vals,
                        "volume": v_vals,
                        "implantation.date": implantation_date,
                        "mouse.strain": mouse_strain
                    })

                    # Clean data: Remove negative time and non-numeric measurements
                    rep_df = rep_df[rep_df["time"] >= 0]
                    rep_df["width"] = pd.to_numeric(rep_df["width"], errors="coerce")
                    rep_df["length"] = pd.to_numeric(rep_df["length"], errors="coerce")
                    rep_df["volume"] = pd.to_numeric(rep_df["volume"], errors="coerce")
                    
                    # Drop rows where critical measurements (W/L) are missing
                    rep_df.dropna(subset=["width", "length"], inplace=True)

                    # Impute Volume: V = (W^2 * L) / 2 if missing but W/L are present
                    v_mask = rep_df["volume"].isna()
                    rep_df.loc[v_mask, "volume"] = (
                        (rep_df.loc[v_mask, "width"] ** 2 * rep_df.loc[v_mask, "length"]) / 2
                    ).round(2)

                    processed_dfs.append(rep_df)

    except Exception as e:
        logging.error(f"Critical error parsing {file_path}: {e}")
        # We re-raise to ensure the pipeline stops and the error is logged in run_process.py
        raise e
        
    return processed_dfs

def process_all_files(metadata_file: str, base_dir: str, min_days: int) -> Tuple[List[pd.DataFrame], List[Dict[str, str]]]:
    """
    Orchestrates the processing of the manifest, skipping excluded files.
    """
    with open(metadata_file, "r") as f:
        file_list = json.load(f)
    
    all_data = []
    qc_issues = []
    
    for file_entry in file_list:
        file_path = os.path.join(base_dir, file_entry.get("file_path", ""))
        tag = file_entry.get("tag", "")

        # Skip files tagged for removal or duplicates during the Scan step
        if any(bad in tag for bad in ["in_remove_list", "duplicated_file"]):
            logging.info(f"Skipping {file_entry['file_name']} | Reason: {tag}")
            qc_issues.append({"file": file_path, "issue": f"Excluded: {tag}"})
            continue

        try:
            # Pass the min_days parameter from the config
            dfs = process_file(file_entry, base_dir, min_days=min_days)
            if dfs:
                all_data.extend(dfs)
                qc_issues.append({"file": file_path, "issue": "success"})
            else:
                qc_issues.append({"file": file_path, "issue": "failed_duration_threshold"})
        except Exception as e:
            logging.error(f"Processing failed for {file_path}: {e}")
            qc_issues.append({"file": file_path, "issue": f"Error: {str(e)}"})
    
    return all_data, qc_issues

def save_processed_data(all_data: List[pd.DataFrame], processed_tsv: str) -> None:
    """Combines all processed replicates into a final TSV master file."""
    if all_data:
        master_df = pd.concat(all_data, ignore_index=True)
        master_df.to_csv(processed_tsv, index=False, sep="\t")
        logging.info(f"✅ Master dataset saved to {processed_tsv} ({len(master_df)} rows)")
    else:
        logging.warning("No data blocks met the inclusion criteria. TSV not saved.")

def save_qc_report(qc_issues: List[Dict[str, str]], qc_report_file: str) -> None:
    """Saves a JSON report of the processing results for auditing."""
    with open(qc_report_file, "w") as f:
        json.dump(qc_issues, f, indent=4)
    logging.info(f"📋 QC report generated: {qc_report_file}")
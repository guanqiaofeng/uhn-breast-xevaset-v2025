import pandas as pd
import os
import re
import logging
from typing import List, Dict, Tuple, Any, Optional

# Setup logging to track regex successes/failures
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def update_model_id(samplename: str, filename: str) -> Tuple[str, str, str]:
    """
    Extracts ModelID, SampleID, and Model_Type using UHN legacy regex patterns.
    """
    def _extract_from_parts(parts: List[str], label: str):
        passage_numbers = [p for p in parts if re.match(r"^P\d+$", p)]
        
        if "RES" in parts:
            model_type = "Resistant"
            if len(passage_numbers) > 1:
                first_passage = passage_numbers[0]
                resistant_passage = passage_numbers[-1]
                try:
                    index_second_last = parts.index(passage_numbers[-2])
                    model_part = parts[:index_second_last]
                    resistant_drug = parts[index_second_last + 1]
                    model_id = "_".join(model_part).replace(" ", "").upper() + f"_{first_passage}_{resistant_drug}_RES"
                    sample_id = f"{model_id}_{resistant_passage}"
                except Exception:
                    model_id, sample_id = "UNKNOWN_RES", "UNKNOWN_RES_Pn"
            elif len(passage_numbers) == 1:
                resistant_passage = passage_numbers[0]
                parts.remove(resistant_passage)
                try:
                    index_res = parts.index('RES')
                    resistant_drug = parts[index_res - 1]
                    model_part = parts[:index_res - 1]
                    model_id = "_".join(model_part).replace(" ", "").upper() + f"_{resistant_drug}_RES"
                    sample_id = f"{model_id}_{resistant_passage}"
                except Exception:
                    model_id, sample_id = "UNKNOWN_RES", "UNKNOWN_RES_Pn"
            else:
                model_id = "_".join(parts).replace(" ", "").upper()
                sample_id = model_id + "_Pn"
        else:
            model_type = "Parental"
            if not passage_numbers:
                model_id = "_".join(parts).replace(" ", "").upper()
                sample_id = model_id + "_Pn"
            else:
                model_part = parts[:parts.index(passage_numbers[-1])]
                model_id = "_".join(model_part).replace(" ", "").upper()
                sample_id = f"{model_id}_{passage_numbers[-1]}"

        return model_id, sample_id, model_type

    # --- Standardization Layer ---
    s_name = str(samplename).upper()
    # Apply your lab's specific string cleanup (shortened for readability)
    replacements = {
        r'P[XM]\s*\+?\s*(\d+)': r'P\1', "BPTO.": "BPTO", "BXTO.": "BXTO",
        "RES 945": "945 RES", "ERIBULIN": "ERIB", "TAXOL": "TAX"
    }
    for pat, rep in replacements.items():
        s_name = re.sub(pat, rep, s_name)
    
    # Clean common noise words
    noise = ["#2- COHORT", "TNBC", "AXILLA", "PARENTAL", "NOTCH "]
    for word in noise:
        s_name = s_name.replace(word, "")

    parts = re.split(r"\s+", s_name.strip())
    m_id, s_id, m_type = _extract_from_parts(parts, "samplename")

    # Fallback to filename logic if samplename is generic (e.g., contains timestamp)
    if len(s_name) < 4 or "00:00:00" in s_name:
        f_name = os.path.splitext(filename)[0].upper().replace("COPY OF ", "")
        f_parts = re.split(r"\s+", f_name.strip())
        return _extract_from_parts(f_parts, "filename")
    
    return m_id, s_id, m_type

def dose_number_unit_split(dose_series: pd.Series) -> pd.DataFrame:
    """Splits composite doses (e.g., '10mg/kg+5mg/kg') into numeric and unit parts."""
    def split_logic(dose_str):
        if pd.isna(dose_str): return '', ''
        parts = re.findall(r'([0-9]*\.?[0-9]+)\s*([a-zA-Z/µ%]+)', str(dose_str))
        return '+'.join(p[0] for p in parts), '+'.join(p[1] for p in parts)

    extracted = dose_series.apply(split_logic)
    return pd.DataFrame(extracted.tolist(), columns=['dose', 'dose.unit'])

def assign_mouse_ids(group: pd.DataFrame, batch_value: str) -> pd.DataFrame:
    """Assigns unique mouse IDs (e.g., BATCH.m1) based on Time=0 markers."""
    group = group.copy()
    group["batch"] = batch_value
    counter, current_id, ids = 0, None, []
    for t in group["time"]:
        if t == 0:
            counter += 1
            current_id = f"{batch_value}.m{counter}"
        ids.append(current_id)
    group["model.id"] = ids
    return group

def update_model_data(df: pd.DataFrame) -> pd.DataFrame:
    """Applies regex ID extraction and cleans metadata columns."""
    res = df.apply(lambda r: pd.Series(update_model_id(r["sampleID"], r["file.name"])), axis=1)
    df[["modelID", "sampleID_new", "modelType"]] = res
    df = df.rename(columns={'sampleID': 'sampleID.original', 'sampleID_new': 'sampleID'})
    
    # Standardize strain/location/source
    df['mouse.strain'] = df['mouse.strain'].str.replace("Strain: ", "", regex=False)
    #df['implantation.location'] = df['implantation.location'].str.replace("Location: ", "", regex=False)
    return df

def map_model(df: pd.DataFrame, map_file: str, target_type: str) -> pd.DataFrame:
    """Merges with official Excel mapping to finalize IDs and filter by lineage."""
    map_df = pd.read_excel(map_file).drop(columns=["name", "name_clean"], errors='ignore')
    
    # Filter map by type if requested
    if target_type != 'both':
        map_df = map_df[map_df['modelType_update'].str.lower() == target_type.lower()]

    merged = df.merge(map_df, on="file.name", how="left")
    
    # Update columns with mapped values if they exist
    cols = ['modelID', 'sampleID', 'modelType']
    for col in cols:
        map_col = f"{col}_update"
        if map_col in merged.columns:
            merged[col] = merged[map_col].fillna(merged[col])
    
    # Post-merge filtering for the final set
    if target_type != 'both':
        merged = merged[merged['modelType'].str.lower() == target_type.lower()]
        
    return merged.rename(columns={'modelType': 'model.type'}).drop(columns=[c for c in merged.columns if '_update' in c])

def update_drug_data(df: pd.DataFrame, drug_map_file: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Standardizes drug names and doses, returning both valid and filtered data."""
    drug_map = pd.read_excel(drug_map_file).astype({'drug': 'string', 'dose': 'string'})
    
    merged = drug_map.merge(df, on=['drug', 'dose'], how='right')
    merged = merged.rename(columns={'treatment.standardized': 'drug.id'})
    
    # Keep track of originals
    merged['drug.original'] = merged['drug']
    merged['drug'] = merged['drug.id']
    merged['dose.original'] = merged['dose']
    
    valid = merged[merged['keep'] == "yes"].copy()
    invalid = merged[merged['keep'] != "yes"].copy()
    
    # Split dose numeric/unit
    if not valid.empty:
        valid[['dose', 'dose.unit']] = dose_number_unit_split(valid['dose.standardized'])
    
    # --- DYNAMIC COLUMN DROPPING ---
    # 1. Identify all columns that contain 'comment' (case-insensitive)
    cols_to_drop = [c for c in valid.columns if 'comment' in str(c).lower()]
    
    # 2. Add your specific technical columns to the drop list
    cols_to_drop.extend(['keep', 'check_dup', 'dose.standardized'])
    
    # 3. Drop them from the valid dataframe
    valid = valid.drop(columns=cols_to_drop, errors='ignore')
    # -------------------------------
        
    return valid, invalid

def update_batch_data(df: pd.DataFrame) -> pd.DataFrame:
    """Groups observations into batches and assigns specific mouse identifiers."""
    df["batch"] = df["modelID"].astype(str) + "." + df["drug.id"].astype(str)
    return df.groupby("batch", sort=False).apply(
        lambda g: assign_mouse_ids(g, g.name), include_groups=False
    ).reset_index(drop=True)
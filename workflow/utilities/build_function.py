import pandas as pd
import re
import os
import logging
from pathlib import Path
from typing import Tuple, Dict

def build_model(df: pd.DataFrame) -> pd.DataFrame:
    """Builds the unique biological model metadata for XevaSet."""
    # We want unique biological models, not unique mice. 
    # Usually, 'modelID' or 'model.id' stripped of passage info represents the model.
    cols = ['model.id', "drug",'implantation.location', 'implantation.tissue.source']
    existing_cols = [c for c in cols if c in df.columns]
    
    # We drop duplicates
    model_df = df[existing_cols].drop_duplicates().reset_index(drop=True)
    
    # Rename modelID to model.id for Xeva consistency
    model_df = model_df.rename(columns={'modelID': 'model.id'})
    
    model_df['tissue'] = 'BRCA'
    model_df['tissue.name'] = 'Breast Cancer'
    
    # Extract Patient ID
    model_df['patient.id'] = model_df['model.id'].apply(
        lambda x: re.split(r'_P\d+', str(x).split('.')[0])[0]
    )
    return model_df

def build_drug(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregates drug metadata, handling potential multiple doses per drug ID."""
    cols = [
        'drug.id', 'drugname.standardized', 'drug.alternativename', 'dose', 'dose.unit', 'drug.original', 
        'dose.original', 'frequency' ,'treatment.day','is.single', 'is.control','is.ADC','is.PDC','paired.treatment'
    ]
    existing = [c for c in cols if c in df.columns]
    drug_df = df[existing].copy().drop_duplicates()

    def merge_values(series):
        unique_vals = series.dropna().unique()
        return unique_vals[0] if len(unique_vals) == 1 else ';'.join(sorted(set(map(str, unique_vals))))

    return drug_df.groupby("drug.id", dropna=False).agg(merge_values).reset_index()

def build_experiment(df: pd.DataFrame) -> pd.DataFrame:
    """Creates the core experimental observation table."""
    columns = [
        'model.id', 'drug.id', 'dose', 'time', 'volume', 'width', 'length', 
        'batch', 'sampleID', 'model.type', 'date', 'file.name'
    ]
    existing = [c for c in columns if c in df.columns]
    return df[existing].copy()

def build_expDesign(df: pd.DataFrame) -> pd.DataFrame:
    """Maps treatment batches to control (H2O) batches."""
    design = df[['batch', 'model.id']].drop_duplicates().sort_values(['batch', 'model.id'])
    
    # Group mice into batches
    exp_design = design.groupby('batch')['model.id'].apply(lambda x: ','.join(map(str, x))).reset_index()
    exp_design.rename(columns={'batch': 'batch.name', 'model.id': 'treatment'}, inplace=True)
    
    # Map Control: e.g., 'MC1.CISPLATIN' looks for 'MC1.H2O'
    exp_design['prefix'] = exp_design['batch.name'].str.split('.').str[0]
    exp_design['control_key'] = exp_design['prefix'] + '.H2O'
    
    control_map = exp_design.set_index('batch.name')['treatment'].to_dict()
    exp_design['control'] = exp_design['control_key'].map(control_map)
    
    # Remove the H2O batches from the treatment list (they are now in the 'control' column)
    return exp_design[~exp_design['batch.name'].str.endswith('.H2O')].drop(columns=['prefix', 'control_key'])

def build_modToBiobaseMap_pdata(df: pd.DataFrame, omics_dir: Path) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """Maps Model IDs to Omics Sample IDs while ignoring metadata files."""
    biobase_df = df[['model.id', 'modelID']].drop_duplicates()
    biobase_df.columns = ['model.id', 'biobase.id']
    biobase_df['biobase.id'] = biobase_df['biobase.id'].apply(
        lambda x: f'S{x}' if str(x)[0].isdigit() else str(x)
    )

    RNA_KEYWORDS = ['rna', 'tpm']
    MUT_KEYWORDS = ['mutation', 'mut']
    CNV_KEYWORDS = ['cnv']
    
    # NEW: Files to strictly ignore
    IGNORE_KEYWORDS = ['feature', 'metadata', 'annotation', 'probe']

    omics_map_list = []
    pdata_dict = {}

    if not omics_dir.exists():
        return pd.DataFrame(), {}

    for file_path in omics_dir.iterdir():
        if not file_path.is_file(): continue
        
        file_name = file_path.name.lower()
        
        # --- CRITICAL FILTER: Skip metadata files ---
        if any(ig in file_name for ig in IGNORE_KEYWORDS):
            logging.info(f"⏩ Skipping metadata file: {file_path.name}")
            continue

        # 1. Identify Data Type
        if any(k in file_name for k in RNA_KEYWORDS): dtype = 'RNASeq'
        elif any(k in file_name for k in MUT_KEYWORDS): dtype = 'mutation'
        elif any(k in file_name for k in CNV_KEYWORDS): dtype = 'CNV'
        else: continue

        try:
            sep = ',' if file_path.suffix == '.csv' else '\t'
            # 2. Extract Headers
            sample_ids = pd.read_csv(file_path, sep=sep, nrows=0).columns[1:]
            
            # 3. Validation: Ensure we didn't grab 'ensembl' as a sample ID
            invalid_headers = ['ensembl', 'hugo', 'symbol', 'description', 'gene']
            sample_ids = [s for s in sample_ids if not any(inv in str(s).lower() for inv in invalid_headers)]

            if len(sample_ids) == 0:
                logging.warning(f"⚠️ No sample IDs found in {file_path.name}. Skipping.")
                continue

            # 4. Store mapping
            temp_map = pd.DataFrame({'biobase.id': sample_ids, 'mDataType': dtype})
            omics_map_list.append(temp_map)
            
            # 5. Build pData
            pdata = pd.DataFrame({'biobase.id': sample_ids})
            pdata['patient.id'] = pdata['biobase.id'].astype(str).str.replace(r'^S', '', regex=True)
            pdata['tissue'] = 'BRCA'
            pdata['tissue_name'] = "Breast Cancer"
            pdata_dict[dtype] = pdata
            
            logging.info(f"✔️ Successfully mapped {dtype} from: {file_path.name}")
            
        except Exception as e:
            logging.error(f"❌ Error processing {file_path.name}: {e}")

    if omics_map_list:
        all_omics = pd.concat(omics_map_list, ignore_index=True)
        merged_map = biobase_df.merge(all_omics, on='biobase.id', how='inner')
        return merged_map, pdata_dict
    
    return pd.DataFrame(), {}

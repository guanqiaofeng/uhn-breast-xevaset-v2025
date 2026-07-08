import sys
import logging
import pandas as pd
import numpy as np
from pathlib import Path
import argparse

# 1. SETUP PROJECT PATHS
PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.append(str(PROJECT_ROOT))

try:
    from workflow.utilities.curve_function import calculate_mouse_metrics, call_mrecist
except ImportError:
    # Fallback for local execution
    sys.path.append(str(PROJECT_ROOT / "workflow" / "utilities"))
    from curve_function import calculate_mouse_metrics, call_mrecist

def run_analysis(exp_path, design_path, out_batch, out_model, tgi_win, recist_win):
    """
    Main logic for metric calculation. 
    Includes Control/Treatment breakdown for AUC and Slope.
    """
    logging.info(f"📖 Loading experiment data: {exp_path.name}")
    
    if not exp_path.exists() or not design_path.exists():
        logging.error(f"❌ Critical Error: Input files missing.")
        return

    exp_df = pd.read_csv(exp_path)
    design_df = pd.read_csv(design_path)
    
    model_level_results = []
    batch_level_results = []

    logging.info(f"Processing {len(design_df)} treatment batches...")

    for _, row in design_df.iterrows():
        batch_id = row['batch.name']
        trt_ids = [mid.strip() for mid in str(row['treatment']).split(',') if mid.strip()]
        ctrl_ids = [mid.strip() for mid in str(row['control']).split(',') if mid.strip()]
        
        # --- STEP A: BATCH CONTROL BASELINE (Flexible V0) ---
        ctrl_metrics = []
        for cid in ctrl_ids:
            c_full = exp_df[exp_df['model.id'] == cid].sort_values('time')
            if c_full.empty: continue
            
            # Earliest measurement within first 3 days
            cv0_candidates = c_full[c_full['time'] <= 3]
            if cv0_candidates.empty: continue
            cv0 = cv0_candidates.iloc[0]['volume']
            if pd.isna(cv0) or cv0 <= 0: continue

            res = calculate_mouse_metrics(c_full, cv0, tgi_win)
            if res: ctrl_metrics.append(res)

        # Baseline stats for the batch (Control Group)
        if ctrl_metrics:
            mean_ctrl_auc = np.mean([m['auc'] for m in ctrl_metrics])
            mean_ctrl_slope = np.mean([m['slope'] for m in ctrl_metrics])
            mean_ctrl_last_vol = np.mean([m['last_vol'] for m in ctrl_metrics])
        else:
            mean_ctrl_auc = mean_ctrl_slope = mean_ctrl_last_vol = np.nan

        batch_replicates = []

        # --- STEP B: MODEL LEVEL (Treatment Group) ---
        for mid in trt_ids:
            m_full = exp_df[exp_df['model.id'] == mid].sort_values('time')
            if m_full.empty: continue
            
            # Earliest measurement within first 3 days
            v0_candidates = m_full[m_full['time'] <= 3]
            if v0_candidates.empty: continue
            v0 = v0_candidates.iloc[0]['volume']
            
            if pd.isna(v0) or v0 <= 0: continue

            # Independent windows: 60d for RECIST, 28d for Efficacy/Geometric
            m_res_resp = calculate_mouse_metrics(m_full, v0, recist_win)
            m_res_tgi = calculate_mouse_metrics(m_full, v0, tgi_win)

            if m_res_resp or m_res_tgi:
                # Calculate TGI vs Batch Mean
                tgi_val = np.nan
                if m_res_tgi and not pd.isna(mean_ctrl_last_vol):
                    tgi_val = 1 - (m_res_tgi['last_vol'] / mean_ctrl_last_vol)
                
                mouse_data = {
                    'model.id': mid,
                    'batch.name': batch_id,
                    'best.response': m_res_resp['br'] if m_res_resp else np.nan,
                    'best.average.response': m_res_resp['bar'] if m_res_resp else np.nan,
                    'mRECIST': call_mrecist(m_res_resp['br'], m_res_resp['bar']) if m_res_resp else "NA",
                    'TGI': tgi_val,
                    'AUC': m_res_tgi['auc'] if m_res_tgi else np.nan,
                    'slope': m_res_tgi['slope'] if m_res_tgi else np.nan
                }
                model_level_results.append(mouse_data)
                batch_replicates.append(mouse_data)

        # --- STEP C: BATCH AGGREGATION (Xeva-Aligned Statistics) ---
        if batch_replicates:
            rep_df = pd.DataFrame(batch_replicates)
            
            # Median Response (Drops NaNs for calculation)
            m_br = rep_df['best.response'].dropna().median()
            m_bar = rep_df['best.average.response'].dropna().median()
            
            # Group Means for Treatment
            avg_trt_auc = rep_df['AUC'].dropna().mean()
            avg_trt_slope = rep_df['slope'].dropna().mean()

            batch_level_results.append({
                'batch.name': batch_id,
                'TGI': rep_df['TGI'].dropna().median(),
                
                # Comparative Metrics (Nature Paper Focus)
                'abc': mean_ctrl_auc - avg_trt_auc,
                # slope is already an angle in degrees (see calculate_mouse_metrics),
                # so the batch angle is a plain difference, matching Xeva's angle()
                'angle': mean_ctrl_slope - avg_trt_slope,
                
                # Internal Xeva Slots (Breakdown)
                'auc.control': mean_ctrl_auc,
                'auc.treatment': avg_trt_auc,
                'slope.control': mean_ctrl_slope,
                'slope.treatment': avg_trt_slope,

                # Batch-level Response Calls
                'BR_median': m_br,
                'BAR_median': m_bar,
                'mRECIST': call_mrecist(m_br, m_bar) if not (pd.isna(m_br) or pd.isna(m_bar)) else "NA",
                'n_replicates': len(rep_df)
            })

    # --- SAVE OUTPUTS ---
    out_batch.parent.mkdir(parents=True, exist_ok=True)
    df_model = pd.DataFrame(model_level_results)
    df_batch = pd.DataFrame(batch_level_results)
    
    df_model.to_csv(out_model, index=False)
    df_batch.to_csv(out_batch, index=False)
    
    # --- FINAL INTEGRITY SUMMARY ---
    logging.info("=" * 45)
    logging.info("       METRICS STEP: FINAL SUMMARY")
    logging.info("=" * 45)
    logging.info(f"Mice Processed (Model Level):  {len(df_model)}")
    logging.info(f"Batches Processed (Batch Level): {len(df_batch)}")
    
    if not df_batch.empty:
        counts = df_batch['mRECIST'].value_counts()
        logging.info("-" * 45)
        logging.info("BATCH-LEVEL mRECIST SUMMARY:")
        for call in ['CR', 'PR', 'SD', 'PD']:
            logging.info(f"  - {call}: {counts.get(call, 0)}")
        
        valid_abc = df_batch['abc'].dropna()
        if not valid_abc.empty:
            logging.info("-" * 45)
            logging.info(f"Avg Batch Stats: AUC_Ctrl={df_batch['auc.control'].mean():.2f}, AUC_Trt={df_batch['auc.treatment'].mean():.2f}")
            logging.info(f"Avg Comparative: ABC={valid_abc.mean():.2f}, Angle={df_batch['angle'].mean():.2f}")
    logging.info("=" * 45)

def main():
    PROC_DIR = PROJECT_ROOT / "data" / "procdata" / "drug_screen"
    LOG_DIR = PROJECT_ROOT / "logs" / "drug_screen"
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    parser = argparse.ArgumentParser(description="Calculate PDX Metrics with Xeva Alignment")
    parser.add_argument("--exp", default=PROC_DIR / "build" / "experiment.csv")
    parser.add_argument("--design", default=PROC_DIR / "build" / "expDesign.csv")
    parser.add_argument("--out_batch", default=PROC_DIR / "build" / "batch_metrics.csv")
    parser.add_argument("--out_model", default=PROC_DIR / "build" / "model_metrics.csv")
    parser.add_argument("--tgi_win", type=int, default=28)
    parser.add_argument("--recist_win", type=int, default=60)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(LOG_DIR / "calculate_metrics.log"), logging.StreamHandler(sys.stdout)],
        force=True
    )

    logging.info("📈 Running Calculation (Flexible Baseline + Xeva Breakdown)...")
    run_analysis(
        Path(args.exp), Path(args.design), 
        Path(args.out_batch), Path(args.out_model), 
        args.tgi_win, args.recist_win
    )

if __name__ == "__main__":
    main()
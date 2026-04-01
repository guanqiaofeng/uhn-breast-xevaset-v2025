import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

def calculate_mouse_metrics(m_df, v0, window):
    """
    Unified metric calculation:
    - RECIST: Uses best response and best average (3-point mean).
    - Geometric: Uses Linear Regression for slope and trapezoid for AUC.
    """
    # 1. Filter by window
    m_win = m_df[m_df['time'] <= window].sort_values('time').copy()
    if m_win.empty or len(m_win) < 2:
        return None
        
    vols = m_win['volume'].values
    times = m_win['time'].values
    
    # --- GEOMETRIC MATH ---
    # Slope for Angle calculation
    model = LinearRegression().fit(times.reshape(-1, 1), vols)
    slope = model.coef_[0]
    
    # AUC calculation
    auc_val = np.trapezoid(vols, times)

    # --- RECIST MATH ---
    # We only look at volumes AFTER the start (time > 0) to find response
    # If the study is short, we include all volumes to avoid NaNs
    post_v = m_win.loc[m_win['time'] > 0, 'volume'].dropna().values
    if len(post_v) == 0:
        post_v = vols # Fallback to all volumes if no post-day-0 data exists

    # Best Response (BR): The single deepest drop
    br = ((np.min(post_v) / v0) - 1) * 100
    
    # Best Average Response (BAR): Mean of 3 smallest volumes
    # This prevents one noisy low-reading from creating a fake "CR"
    sorted_v = np.sort(post_v)
    bar = ((np.mean(sorted_v[:3]) / v0) - 1) * 100

    return {
        'slope': slope,
        'auc': auc_val,
        'last_vol': vols[-1],
        'br': br,
        'bar': bar
    }

def call_mrecist(br, bar):
    """
    Applies standard Gao et al. (Nature Medicine 2015) thresholds.
    """
    if pd.isna(br) or pd.isna(bar): 
        return "NA"
    
    # CR: Both deep drop and sustained low average
    if br <= -95 and bar <= -40: 
        return "CR"
    # PR: Significant drop and sustained partial low
    if br <= -50 and bar <= -20: 
        return "PR"
    # SD: Growth is controlled within 35% of baseline
    if br <= 35 and bar <= 30: 
        return "SD"
    
    # PD: Growth exceeds 35% from baseline
    return "PD"
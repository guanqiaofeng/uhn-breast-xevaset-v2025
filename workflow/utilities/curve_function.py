import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

def calculate_mouse_metrics(m_df, v0, window, min_time=10):
    """
    Unified metric calculation:
    - RECIST: Best response (min % change) and best average response
      (min running-average % change), both restricted to time >= min_time.
    - Geometric: Uses Linear Regression for slope and trapezoid for AUC.
    """
    # 1. Filter by window
    m_win = m_df[m_df['time'] <= window].sort_values('time').copy()
    if m_win.empty or len(m_win) < 2:
        return None
        
    vols = m_win['volume'].values
    times = m_win['time'].values

    # --- GEOMETRIC MATH ---
    # Normalize to baseline so AUC/TGI reflect relative growth rather than
    # absolute tumor size, consistent with how BR/BAR are normalized below
    vols_norm = vols / v0

    # Slope/angle: baseline-subtract (not divide) and force the regression
    # through the origin, matching Xeva's slope() — the fitted line passes
    # through the mouse's actual first observation, in absolute mm^3/day
    # units. The angle (atan, in degrees) is the value itself, as in Xeva.
    times_shift = times - times[0]
    vols_shift = vols - vols[0]
    model = LinearRegression(fit_intercept=False).fit(
        times_shift.reshape(-1, 1), vols_shift
    )
    slope = np.degrees(np.arctan(model.coef_[0]))

    # AUC calculation (raw volume — mice are size-matched at treatment
    # start, so baseline normalization isn't needed; matches Xeva's AUC())
    auc_val = np.trapezoid(vols, times)

    # --- RECIST MATH (Gao et al. 2015 Nature Medicine mRECIST definition) ---
    # % change in volume from baseline at every time point
    pct_change = (vols / v0 - 1) * 100
    # Running average of % change from day 0 through each time point
    running_avg = np.cumsum(pct_change) / np.arange(1, len(pct_change) + 1)

    # Both BR and BAR are only considered beyond min_time, to avoid an
    # early transient reading being mistaken for a durable response
    eligible = times >= min_time
    if eligible.any():
        # Best Response (BR): deepest single-point drop after min_time
        br = np.min(pct_change[eligible])
        # Best Average Response (BAR): lowest running average after min_time,
        # rewarding sustained response rather than isolated low readings
        bar = np.min(running_avg[eligible])
    else:
        br = np.nan
        bar = np.nan

    return {
        'slope': slope,
        'auc': auc_val,
        # Relative to v0 so TGI (T_last/T0 vs C_last/C0) isn't biased by
        # baseline volume differences between arms
        'last_vol': vols_norm[-1],
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
    if br < -95 and bar < -40:
        return "CR"
    # PR: Significant drop and sustained partial low
    if br < -50 and bar < -20:
        return "PR"
    # SD: Growth is controlled within 35% of baseline
    if br < 35 and bar < 30:
        return "SD"
    
    # PD: Growth exceeds 35% from baseline
    return "PD"
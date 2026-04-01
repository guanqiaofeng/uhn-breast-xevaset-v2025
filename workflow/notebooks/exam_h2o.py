import pandas as pd
import numpy as np

def find_h2o_survival_cutoff(exp_path, control_keyword="H2O"):
    # 1. Load your experiment data
    df = pd.read_csv(exp_path)
    
    # 2. Filter for H2O (Control) mice only
    # (Using case-insensitive search to catch H2O, h2o, Water, etc.)
    h2o_df = df[df['drug.id'].str.contains(control_keyword, case=False, na=False)].copy()
    
    if h2o_df.empty:
        print(f"❌ Error: No drug names containing '{control_keyword}' found in experiment.csv")
        return

    # 3. Find the maximum duration (last day) for each individual mouse
    mouse_durations = h2o_df.groupby('model.id')['time'].max().reset_index()
    mouse_durations.columns = ['model.id', 'duration']
    
    # 4. Sort from lowest (shortest duration) to highest (longest duration)
    sorted_mice = mouse_durations.sort_values(by='duration').reset_index(drop=True)
    
    # 5. Calculate the 80% survival value
    # To find the day where 80% are still "alive," we exclude the bottom 20% of drop-outs.
    n_mice = len(sorted_mice)
    # The 20th percentile index represents the point where 80% of data remains
    cutoff_index = int(np.floor(n_mice * 0.2)) 
    survival_80_value = sorted_mice.loc[cutoff_index, 'duration']
    
    # 6. Output the results
    print("\n" + "="*50)
    print(f"💧 H2O COHORT DURATION ANALYSIS (N={n_mice})")
    print("="*50)
    print(sorted_mice.to_string(index=False))
    print("-" * 50)
    print(f"✅ 80% SURVIVAL CUTOFF VALUE: {survival_80_value} days")
    print(f"   (80% of H2O mice survived until or past this day)")
    print("="*50)

    return survival_80_value

# Run the analysis
cutoff = find_h2o_survival_cutoff("../../data/procdata/drug_screen/build/experiment.csv")
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_survival_cutoff(exp_path, control_name="H2O"):
    # 1. Load data
    df = pd.read_csv(exp_path)
    
    # 2. Identify unique mice and their last recorded day
    # We group by drug to see if different treatments have different drop-off rates
    survival_data = []
    
    for drug in df['drug.id'].unique():
        drug_df = df[df['drug.id'] == drug]
        total_mice = drug_df['model.id'].nunique()
        
        # Get all unique timepoints in the study
        all_times = sorted(df['time'].unique())
        
        for t in all_times:
            # Count how many mice have data at this time OR later
            # (This handles cases where a mouse might miss one day but return)
            mice_at_time = drug_df[drug_df['time'] >= t]['model.id'].nunique()
            survival_pct = (mice_at_time / total_mice) * 100
            
            survival_data.append({
                'Drug': drug,
                'Time': t,
                'Survival_Pct': survival_pct,
                'Is_Control': (control_name.upper() in drug.upper())
            })
    
    survival_df = pd.DataFrame(survival_data)
    
    # 3. Plotting
    plt.figure(figsize=(12, 7))
    sns.set_style("whitegrid")
    
    # Plot all drugs in light gray
    sns.lineplot(data=survival_df[~survival_df['Is_Control']], 
                 x='Time', y='Survival_Pct', units='Drug', estimator=None, 
                 color='lightgray', alpha=0.5, label='Treatment Groups')
    
    # Highlight the Control Group(s) in Red
    sns.lineplot(data=survival_df[survival_df['Is_Control']], 
                 x='Time', y='Survival_Pct', color='red', linewidth=3, label='Control (H2O/Baseline)')
    
    # 4. Add the 80% Cutoff Line
    plt.axhline(y=80, color='black', linestyle='--', linewidth=2)
    plt.text(0, 82, '80% Survival Threshold', fontsize=12, fontweight='bold')
    
    # 5. Formatting
    plt.title('PDX Cohort Data Presence Over Time', fontsize=16)
    plt.xlabel('Days on Study', fontsize=14)
    plt.ylabel('% Mice with Data Remaining', fontsize=14)
    plt.ylim(0, 105)
    plt.xlim(0, survival_df['Time'].max())
    plt.legend(loc='lower left')
    
    plt.tight_layout()
    plt.savefig("pdx_survival_cutoff_check.png", dpi=300)
    print("📊 Plot saved as 'pdx_survival_cutoff_check.png'")

# Execute
plot_survival_cutoff("../../data/procdata/drug_screen/build/experiment.csv")
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_duration_histogram(exp_path):
    # 1. Load data
    df = pd.read_csv(exp_path)
    
    # 2. Calculate duration (last day) for every individual mouse
    duration_df = df.groupby(['drug.id', 'model.id'])['time'].max().reset_index()
    duration_df.columns = ['Drug', 'Mouse_ID', 'Duration']
    
    # 3. Plotting
    plt.figure(figsize=(14, 8))
    sns.set_style("white")
    
    # Create the histogram
    # We use 'hue' to distinguish drugs, but 'multiple="stack"' or 'element="step"' helps readability
    g = sns.histplot(data=duration_df, x='Duration', hue='Drug', 
                     element="step", kde=True, palette="viridis", alpha=0.5)
    
    # 4. Identify and Draw the Median
    overall_median = duration_df['Duration'].median()
    plt.axvline(overall_median, color='red', linestyle='--', linewidth=2, label=f'Overall Median: {overall_median} days')
    
    # Optional: Per-drug medians (printed to console)
    print("\n" + "="*40)
    print("🕒 MEDIAN TREATMENT DURATION BY DRUG")
    print("="*40)
    medians = duration_df.groupby('Drug')['Duration'].median().sort_values()
    print(medians)
    medians.to_csv("median_duration_per_drug.tsv", sep="\t")
    
    # 5. Formatting
    plt.title('Distribution of Treatment Durations across all Mice & Drugs', fontsize=16)
    plt.xlabel('Duration (Days)', fontsize=14)
    plt.ylabel('Count (Number of Mice)', fontsize=14)
    plt.legend(title="Drug Group", bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig("pdx_duration_histogram.png", dpi=300)
    print("\n📊 Histogram saved as: pdx_duration_histogram.png")

# Execute
plot_duration_histogram("../../data/procdata/drug_screen/build/experiment.csv")
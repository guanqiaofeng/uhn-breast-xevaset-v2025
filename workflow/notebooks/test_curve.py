import pandas as pd
df_batch = pd.read_csv("../../data/procdata/drug_screen/build/batch_metrics.csv")
df_model = pd.read_csv("../../data/procdata/drug_screen/build/model_metrics.csv")

print("--- PD Batch Analysis ---")
pd_batches = df_batch[df_batch['mRECIST_batch_call'] == 'PD']['batch.name'].head(5)
for b in pd_batches:
    models = df_model[df_model['batch.name'] == b]
    print(f"Batch: {b}")
    print(models[['model.id', 'BR_response', 'mRECIST_call']])
    print("-" * 20)
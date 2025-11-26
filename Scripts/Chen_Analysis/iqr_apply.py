import pandas as pd
import os

df = pd.read_csv('Data/Chen_Analysis/chen_data_5.csv')
df_iqr = pd.read_csv('Data/Chen_Analysis/summary_iqr.csv')

files_to_keep = set(df_iqr['file'].tolist())

fasta_dir = "Data/Chen_Analysis/Bins"

for filename in os.listdir(fasta_dir):
    if filename.endswith(".fa.gz"):
        sample_name = filename[:-6] 
        
        
        if sample_name not in files_to_keep:
            file_path = os.path.join(fasta_dir, filename)
            os.remove(file_path)
            print(f"Removido: {file_path}")

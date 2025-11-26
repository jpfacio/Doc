import pandas as pd
import os

df = pd.read_csv('../Data/data_chen_5.csv')
df_iqr = pd.read_csv('../Data/summary_iqr.csv')

files = set(df_iqr['file'].to_list())

df_new = df[df['Sample name'].isin(files)]

fasta_dir = "../Data/Bins"

for i in files:
    file_path = os.path.join(fasta_dir, f"{i}.fa.gz")
    if os.path.exists(file_path):
        os.remove(file_path)
    else:
        print("File not found")

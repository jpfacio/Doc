import pandas as pd
import os
import subprocess
import tempfile
import gzip
import shutil

df = pd.read_csv("../../Data/Chen_Analysis/chen_data_qc.csv")
bins_dir = "../../Data/Chen_Analysis/Bins" 
prokka_dir = "../../Data/Chen_Analysis/Prokka"

os.makedirs(f"{prokka_dir}/Archaea", exist_ok=True)
os.makedirs(f"{prokka_dir}/Bacteria", exist_ok=True)

df_archaea = df[df['species'].str[3] == "A"]
df_bacteria = df[df['species'].str[3] == "B"]

def extract_genus(tax_string):
    if pd.isna(tax_string):
        return ""
    tax_class = tax_string.split(';')
    if len(tax_class) > 5 and tax_class[5].startswith('g__'):
        genus = tax_class[5][3:]
        return genus if genus else ""
    return ""

df_archaea = df_archaea.copy()
df_bacteria = df_bacteria.copy()
df_archaea['genus'] = df_archaea['species'].apply(extract_genus)
df_bacteria['genus'] = df_bacteria['species'].apply(extract_genus)

archaea_w_genus = df_archaea[df_archaea['genus'] != '']
archaea_n_genus = df_archaea[df_archaea['genus'] == '']
bacteria_w_genus = df_bacteria[df_bacteria['genus'] != '']
bacteria_n_genus = df_bacteria[df_bacteria['genus'] == '']

def decompress_to_temp(gz_path, temp_dir, sample_name):
    fa_path = os.path.join(temp_dir, f"{sample_name}.fa")
    try:
        with gzip.open(gz_path, 'rt') as f_in:
            with open(fa_path, 'w') as f_out:
                shutil.copyfileobj(f_in, f_out)
        return fa_path
    except Exception as e:
        print(f"Unzip error {gz_path}: {e}")
        return None

def run_prokka(sample_name, genus, is_archaea, output_dir, fa_path):
    kingdom = "Archaea" if is_archaea else "Bacteria"
    genus_name = genus if genus else "Candidatus"
    
    cmd = [
        "prokka", "--outdir", output_dir,
        "--prefix", sample_name,
        "--kingdom", kingdom,
        "--genus", genus_name,
        "--strain", sample_name,
        "--cpus", "4",
        "--force", fa_path
    ]
    
    print(f"Running Prokka: {sample_name} ({kingdom}, gênero: {genus_name})")
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print(f"Error: {sample_name}")

with tempfile.TemporaryDirectory() as temp_dir:
    print(f"Temp dir: {temp_dir}")
    
    def process_group(df_group, is_archaea):
        for _, row in df_group.iterrows():
            sample = row['sample_name']
            gz_path = os.path.join(bins_dir, f"{sample}.fa.gz")
            
            if not os.path.exists(gz_path):
                print(f"File not found: {gz_path}")
                continue
                
            fa_path = decompress_to_temp(gz_path, temp_dir, sample)
            if not fa_path:
                continue
                
            domain = "Archaea" if is_archaea else "Bacteria"
            output = f"{prokka_dir}/{domain}/{sample}"
    
            genus_val = row['genus'] if 'genus' in row else ""
            run_prokka(sample, genus_val, is_archaea, output, fa_path)

    process_group(archaea_w_genus, True)
    process_group(archaea_n_genus, True)
    process_group(bacteria_w_genus, False)
    process_group(bacteria_n_genus, False)

print("Done")


        



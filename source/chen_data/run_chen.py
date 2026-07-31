from pathlib import Path
import subprocess
import functions as f
import pandas as pd
from time import perf_counter
from datetime import timedelta
import resource

# Control keys

get_json = False
metadata = False

# Defining directories and files

scripts_dir = Path("source/chen_data")
data_dir = Path("Data/Raw/Bins")
tmp=Path("tmp")
csv_file = tmp / "chen_data.csv"
json_file =  tmp / "chen_data.json"
log=Path("log")
log_run=Path("log/run.log")

# Download and organizing Data

if get_json:
    
    print("Downloading and processing the Chen et al. (2022) database")
    print("Downloading metadata json")
    
    fetch_start = perf_counter()
    
    subprocess.run([
        "wget", 
        "https://db.cngb.org/maya/api/get_export_data?dataset_id=MDB0000002&table_name=sample",
        "-O",
       str(json_file)
    ], check=True)
    
    print("Converting the json metadata to csv")
    
    chen_df = f.convert_json_csv(json_file, csv_file)
    
    print("Done!")

    print('Filtering metadata')

    chen_df_filtered = f.filtering(chen_df, csv_file)

    print("Done!")

    links = f.generate_ftp_list(chen_df_filtered)

    print('Downloading MAGs from metadata')
    
    f.fetch_urls(links, data_dir)
    
    fetch_elapsed =perf_counter() - fetch_start
    peak_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    
    space = subprocess.run(
    ["du", "-sh", '.'],
    capture_output=True,
    text=True,
    check=True
    )
    
    project_size = space.stdout.split()[0]
    
    with open(log_run, 'a') as run:
        run.write(
            "#####  FETCH URLS CHECKPOINT  #####\n\n"
            f"Execution time: {timedelta(seconds=round(fetch_elapsed))}\n"
            f"Peak memory: {peak_mem / 1024:.2f} MB\n"
            f"Project size: {project_size}\n\n"
        )
        
else:
    pass
    
 # Organizing the metadata file
 
if metadata:
    chen_data = pd.read_csv("tmp/chen_data.csv")
    metadata = pd.read_csv("tmp/metadata.csv")

    metadata['bin'] = chen_data['sample_name']
    metadata['sample'] = chen_data["assembly_id"]
    metadata['project'] = chen_data['project_id']
    metadata['id_study'] = "https://doi.org/10.1038/s41586-024-07891-2"
    metadata['coord'] = chen_data['latitude_and_longitude']
    metadata['date'] = chen_data['collected_date']

    metadata.to_csv("tmp/metadata.csv", index=False)
else:
    pass








    

from pathlib import Path
import functions as f
import pandas as pd
from time import perf_counter
from datetime import timedelta
import subprocess
import resource
import requests
import time
import json
import os

# Control Keys

qc = False
bakta_key= False
ent_key = False
go = False
pah_key = True

# Path definitions

data_dir=Path("Data/Raw/Bins")
tmp=Path("tmp")
log=Path("log")
log_run=Path("log/run.log")
metadata=Path("tmp/metadata.csv")


###########################################
##########     QC CHECKPOINT     ##########          
###########################################

if qc: 
    
    print("Generating summary statistics")
    
    st_start = perf_counter()
    
    st_summary = f.st.seqkit_summary(data_dir, tmp)
    
    f.st.seq_filter(st_summary)
    
    f.st.seq_remove_500(data_dir)
    
    checkm_data = f.st.run_checkm(data_dir, tmp, log)
    
    f.st.final_processing(st_summary, checkm_data, tmp)
    
    st_elapsed = perf_counter() - st_start
    st_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    
    st_space = subprocess.run(
        ["du", "-sh", '.'],
        capture_output=True,
        text=True,
        check=True
        )
    
    st_size = st_space.stdout.split()[0]
    
    with open(log_run, 'a') as run:
            run.write(
                "#####   QC CHECKPOINT  #####\n\n"
                f"Execution time: {timedelta(seconds=round(st_elapsed))}\n"
                f"Peak memory: {st_mem / 1024:.2f} MB\n"
                f"Project size: {st_size}\n\n"
            )
            
else:
    pass

##############################################
##########     BAKTA CHECKPOINT     ##########
##############################################

if bakta_key:
    
    print("Starting protein annotation (Bakta)")
    
    bakta_start = perf_counter()
    
    bins_list = f.bakta.path_to_list(data_dir)
    f.bakta.fetch_bakta(bins_list, Path('../db-light'), Path('Data/Raw/Processed'))
    
    bakta_elapsed = perf_counter() - bakta_start
    bakta_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    
    bakta_space = subprocess.run(
        ["du", "-sh", '.'],
        capture_output=True,
        text=True,
        check=True
        )
    
    bakta_size = bakta_space.stdout.split()[0]
    
    with open(log_run, 'a') as run:
            run.write(
                "#####   BAKTA CHECKPOINT  #####\n\n"
                f"Execution time: {timedelta(seconds=round(bakta_elapsed))}\n"
                f"Peak memory: {bakta_mem / 1024:.2f} MB\n"
                f"Project size: {bakta_size}\n\n"
            )
else:
    pass

#################################################
##########     ENTITIES CHECKPOINT     ##########
#################################################

if ent_key:
    
    ent_start = perf_counter()
    
    print("Building main entities")
    
    processed_dir = Path("Data/Raw/Processed")
    
    for file in processed_dir.rglob('*'):
        if file.is_file() == False:
            continue
        
        str_file = str(file)
        if '.tsv' not in str_file and '.gbff' not in str_file:
            for _ in range(10):
                file.unlink()
                if file.exists() == False:
                    break
                time.sleep(1)
        elif 'hypotheticals' in str_file:
            for _ in range(10):
                file.unlink()
                if file.exists() == False:
                    break
                time.sleep(1)
                
    
    genes_ent = pd.concat(
        f.ent.create_genes_ent(tsv) for tsv in processed_dir.rglob("*.tsv")
    )
    genes_ent.to_csv("Data/Entities/genes.csv", index=False)
    
    bins_ent = f.ent.create_bins_ent(metadata)
    bins_ent.to_csv("Data/Entities/bins.csv", index=False)
    
    studies_ent = f.ent.create_studies_ent(metadata)
    studies_ent.to_csv("Data/Entities/studies.csv", index=False)
    
    ent_elapsed = perf_counter() - ent_start
    ent_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    
    ent_space = subprocess.run(
            ["du", "-sh", '.'],
            capture_output=True,
            text=True,
            check=True
            )
    
    ent_size = ent_space.stdout.split()[0]
    
    with open(log_run, 'a') as run:
                run.write(
                    "#####   ENTITIES CHECKPOINT  #####\n\n"
                    f"Execution time: {timedelta(seconds=round(ent_elapsed))}\n"
                    f"Peak memory: {ent_mem / 1024:.2f} MB\n"
                    f"Project size: {ent_size}\n\n"
                )
else: 
    pass

#############################################
##########     GO CHECKPOINT     ##########
#############################################

if go:
    
    print("Running GO annotation module")
    
    go_start = perf_counter()
    
    genes_ent = pd.read_csv("Data/Entities/genes.csv")
    
    annotations = f.go.create_database_metadata(genes_ent)
    
    annotations = f.go.fetch_uniref2uniparc(annotations)
    
    annotations = f.go.fetch_uniparc2interpro(annotations)
    
    ipr2go = f.go.parse_interpro2go("support_files/interpro2go.txt")
    annotations = f.go.fetch_go_terms(annotations, ipr2go)
    
    annotations = f.go.resolve_go_names(annotations, "support_files/go-basic.obo")
    
    annotations.to_csv("Data/Entities/annotations.csv", index=False)
    
    go_elapsed = perf_counter() - go_start
    go_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    
    go_space = subprocess.run(
        ["du", "-sh", '.'],
        capture_output=True,
        text=True,
        check=True
        )
    
    go_size = go_space.stdout.split()[0]
    
    with open(log_run, 'a') as run:
            run.write(
                "#####   GO CHECKPOINT  #####\n\n"
                f"Execution time: {timedelta(seconds=round(go_elapsed))}\n"
                f"Peak memory: {go_mem / 1024:.2f} MB\n"
                f"Project size: {go_size}\n\n"
            )
            
else:
    pass

############################################
##########     PAH CHECKPOINT     ##########
############################################

if pah_key:
    
    print("Running PAH-degradation gene screening (DIAMOND)")
    
    pah_start = perf_counter()
    
    candidates, pah_summary = f.pah.run()
    
    pah_elapsed = perf_counter() - pah_start
    pah_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    
    pah_space = subprocess.run(
        ["du", "-sh", '.'],
        capture_output=True,
        text=True,
        check=True
        )
    
    pah_size = pah_space.stdout.split()[0]
    
    with open(log_run, 'a') as run:
            run.write(
                "#####   PAH CHECKPOINT  #####\n\n"
                f"Execution time: {timedelta(seconds=round(pah_elapsed))}\n"
                f"Peak memory: {pah_mem / 1024:.2f} MB\n"
                f"Project size: {pah_size}\n\n"
            )
            
else:
    pass

   









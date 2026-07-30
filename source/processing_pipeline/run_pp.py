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

qc = True
bakta_key= False
uniprotkb_key = False
basic_info_key = False
go_key = False
uniparc_key = False
annot_key = False
go_analysis_key = False
pah_key = False

data_dir=Path("Data/Raw/Bins")
tmp=Path("tmp")
log=Path("log")
log_run=Path("log/run.log")

print("Generating summary statistics")

if qc: 
    
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

if bakta_key:
    
    print("Starting protein annotation (Bakta)")
    
    bins_list = f.annot.path_to_list(data_dir)
    f.annot.fetch_bakta(bins_list, Path('db-light'), Path('Data/Raw/Processed'))
else:
    pass

if annot_key:
    
    print("Processing protein annotation")
    
    data_annot = Path("Data/Raw/Processed")

    tsv_annot = []
    for i in data_annot.glob("*.tsv"):
        name = i.name
        if "hypotheticals" not in name and "inference" not in name:
            tsv_annot.append(i)

    genes_bins = []
    for i in tsv_annot:
        genes_temp = f.ent.create_genes_ent(i)
        if not genes_temp.empty:
            genes_bins.append(genes_temp)
        
    if genes_bins:
        genes_ent = pd.concat(genes_bins, ignore_index=True)
        genes_ent.to_csv("Data/Entities/genes.csv", index=False)
    
        
    metadata = "tmp/metadata.csv"

    f.ent.create_bins_ent(metadata)

    bins = pd.read_csv("Data/Entities/bins.csv")

    doi_list = bins["Study_ID"].unique().tolist()

    doi_dfs = []
    for doi in doi_list:
        df = f.ent.create_studies_ent(doi)
        doi_dfs.append(df)
        df = pd.read_csv("tmp/uniparc_info.csv")

    final_doi_df = pd.concat(doi_dfs, ignore_index=True)

    final_doi_df.to_csv("Data/Entities/studies.csv", index=False)

    f.annot.create_domain_metadata("Data/Entities/genes.csv")
else:
    pass

if uniprotkb_key:
    uniref_to_uniprotkb = f.annot.map_prot("tmp/domain_metadata.csv")
    uniref_to_uniprotkb.to_csv("tmp/uniref_to_uniprotkb.csv")
else:
    pass

if basic_info_key:
    uniprotkb_info = f.annot.uniprot_info("tmp/uniref_to_uniprotkb.csv")
    uniprotkb_info.to_csv("tmp/uniprotkb_info.csv")
else:
    pass

if go_key:
    go_terms_info = f.annot.get_go_terms("tmp/uniprotkb_info.csv")
else:
    pass

if uniparc_key:
    df = pd.read_csv("tmp/uniparc_info.csv")
    f.annot.uniparc_info("tmp/uniref_to_uniprotkb.csv")
else:
    pass

if go_analysis_key:
    f.annot.download_interpro2go(
    "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro2go",
    "tmp"
    )

    ipr2go = f.annot.interpro2go_parser("tmp/interpro2go")
    ipr2go_complete = f.annot.interpro2go_convert("tmp/uniparc_info.csv", "InterPro", ipr2go)
    ipr2go_complete = ipr2go_complete[["UniParcId", "IPR", "GO_ID", "GO_Description"]]

    ipr2go_complete.to_csv("tmp/ipr2go_complete.csv")

    final_merge_go = f.annot.final_merge_go("tmp/ipr2go_complete.csv")
    final_merge_go.to_csv("tmp/go_terms.csv")

    genes = pd.read_csv("Data/Entities/genes.csv")
    genes = genes.merge(final_merge_go[['Gene_Tag', 'go', 'Description']], on='Gene_Tag', how='left')
    genes = genes[["Gene_Tag", "Contig", "Start", "Stop", "Product", "Databases", "go", "Description", "Bin", "Product_Sequence"]]
    genes.to_csv("Data/Entities/genes.csv")

    go_terms = genes[["Gene_Tag", "go"]]

    df_go_info = f.go.analyze_go_terms(go_terms)
    df_go_info.to_csv('tmp/go_info.csv', index=False)
    
    go_info = pd.read_csv('tmp/go_info.csv')

    go_info_full = f.go.get_hierarchy(go_info)
    go_info_full.to_csv('tmp/go_info_full.csv')
else:
    pass

if pah_key:
    f.candidates.against_pah("Data/Reports/candidates.faa", "support_files/pah_db_v2.4/")
else:
    pass

   










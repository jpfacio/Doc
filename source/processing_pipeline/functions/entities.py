import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from habanero import cn
import json
import re


def create_genes_ent(tsv: Path) -> pd.DataFrame:
    
    """This function is a core element of the database building module.
    It receives a .tsv as input from the Bakta output and filters only
    CDS rows, create the gene tags (BIN_XXXX), saves useful information
    and retrieve the protein sequence from Bakta's .gbff file. The resulting
    DataFrame will be the core data source from all the upcoming processes.
    
    Args:
        tsv (Path): A .tsv file from Bakta output

    Returns:
        pd.DataFrame: Entity dataframe
    """
    df = pd.read_csv(tsv, sep="\t", skiprows=5)
    df.columns.values[0] = df.columns[0].lstrip("#")
    
    df_cds = df[df["Type"] == "cds"]
    n_row = len(df_cds["Sequence Id"])
    
    tsv_base = str(tsv).removesuffix(".tsv")
    tsv_name = Path(tsv_base).name
    
    gene_tags = [f"{tsv_name}_{i:04d}" for i in range(1, n_row + 1)]

    contig = list(df_cds["Sequence Id"])
    start = list(df_cds["Start"])
    stop = list(df_cds["Stop"])
    function = list(df_cds["Product"])
    databases = list(df_cds["DbXrefs"])
    
    tsv_path = Path(tsv)
    gbff = tsv_path.with_suffix(".gbff")
    
    locus_tag = list(df_cds["Locus Tag"])
    
    seq_dict = {}

    for record in SeqIO.parse(gbff, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                locus = feature.qualifiers.get("locus_tag", [None])[0]
                if locus is not None and locus in locus_tag:
                    seq = feature.qualifiers.get("translation", [""])[0]
                    seq_dict[locus] = seq

        
    sequences = [seq_dict.get(locus) for locus in locus_tag]
    
    
    genes_ent = pd.DataFrame({"Gene_Tag": gene_tags,
                               "Contig": contig,
                               "Start": start,
                               "Stop": stop,
                               "Product": function,
                               "Databases": databases,
                               "Product_Sequence": sequences,
                               "Bin": tsv_name})
    
    return genes_ent
    
def create_bins_ent(path: Path) -> pd.DataFrame:
    
    """This function defines the Bins entity of the database. It takes the 
    previously saved metadata file from the tmp directory and store its info
    in an entity containing the bins build from the data mining process
    
    Args:
        path (Path): Filepath to be processed

    Returns:
        pd.DataFrame: Bins entity DataFrame
    """
    metadata = pd.read_csv(path)
    
    bins = list(metadata["bin"])
    sample = list(metadata['sample'])
    project = list(metadata['project'])
    date = list(metadata['date'])
    study = list(metadata['id_study'])
    
    bins_ent = pd.DataFrame({"Bin": bins,
                             "Sample": sample,
                             "Project": project,
                             "Date": date,
                             "Study_ID": study})
    
    return bins_ent

def create_studies_ent(doi: str):
        raw_data = cn.content_negotiation(ids=doi, format="citeproc-json")
        data = json.loads(raw_data)
        
        title = data.get("title", [""])[0] if isinstance(data.get("title"), list) else data.get("title", "")
        
        keywords = "; ".join(data.get("subject", [])) if data.get("subject") else ""
        
        authors_str = ""
        authors = data.get("author", [])
        
        for author in authors:
            family = author.get("family", "")
            given = author.get("given", "")
            initial = given[0] if given else ""
            authors_str += f"{family}, {initial}.; "
        authors_str = authors_str[:-2] if authors_str.endswith("; ") else authors_str
        
        container = data.get("container-title", [""]) if data.get("container-title") else ""
        
        year = None
        issued = data.get("issued") or data.get("published")
        if issued and "date-parts" in issued and len(issued["date-parts"][0]) > 0:
            year = issued["date-parts"][0][0]
    
        affil_list = []
        for author in authors:
            affil = author.get("affiliation", [])
            if affil:
                for a in affil:
                    country = a.get("country", "") 
                    if country and country not in affil_list:
                        affil_list.append(country)
        affil_str = "; ".join(affil_list)
        
        return pd.DataFrame([{
            "Study_ID": doi,
            "Title": title,
            "Keywords": keywords,
            "Authors": authors_str,
            "Journal": container,
            "Year": year,
            "Affiliations": affil_str
        }])
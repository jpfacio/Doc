import subprocess
from pathlib import Path
import pandas as pd
from itertools import islice
import requests
import time
import ast
import os
import re
from concurrent.futures import ThreadPoolExecutor
from functools import partial

def path_to_list(bins: Path) -> list:
    
    """Receives a filepath and turns its content into a list

    Args:
        bins (Path): The filepath to be converted
    Returns:
        list: list of filepaths
    """
    
    files = list(bins.glob("*"))
    
    return files
    
def bakta_analysis(bin: Path, db: str, output: Path) -> None:
    
    """Run Bakta analysis
    
    Args:
        bin (string): Filepath containing the sequences
        db (string): Filepath with the bakta database.
        output (string): Filepath of the output
        
    Returns:
        None
        
    """
    outdir = output / bin.stem
    outdir.mkdir(parents=True, exist_ok=True)
    
    log_file = "log/bakta.log"
    
    with open(log_file, "w") as log:
        subprocess.run([
            "bakta", "--db", 
            str(db),
            "--output",
            str(outdir),
            "--meta", "--threads", "2",
            "--force",
            str(bin)
        ], check=True, stdout=log, stderr=subprocess.STDOUT)

def fetch_bakta(bins: list, db: Path, out: Path) -> None:
    
    with ThreadPoolExecutor(max_workers=5) as executor:
        list(executor.map(partial(bakta_analysis, db=db, output=out), bins))
    
    
def create_domain_metadata(genes: str):
    genes = pd.read_csv(genes)
    
    genes_filt = genes[genes['Databases'].notna()]
    domain = genes_filt[['Gene_Tag', 'Databases']].copy()
    
    domain['Databases'] = domain['Databases'].str.split(',')
    
    db_types =[]
    for i in domain['Databases']:
        pref = [s.split(':')[0] for s in i]
        db_types.append(pref)
    
    db_types = list(set(i for sublist in db_types for i in sublist))
    db_types = list(set(s.replace(' ', '') for s in db_types))
    
    for i in db_types:
        domain[i] = None
    
    for idx, db_list in domain['Databases'].items():
        for entry in db_list:
            prefix, value = map(str.strip, entry.split(':', 1))
            if prefix in db_types:
                domain.loc[idx, prefix] = (
                    value if domain.loc[idx, prefix] is None
                    else domain.loc[idx, prefix] + f",{value}"
                )
        
    domain.to_csv("tmp/domain_metadata.csv", index=False)
    
def map_prot(df_path, batch_size=50, fallback_batch_size=15):
    df = pd.read_csv(df_path)

    mapping = (
        df[["UniRef", "Gene_Tag"]]
        .dropna()
        .drop_duplicates(subset="UniRef")
        .set_index("UniRef")["Gene_Tag"]
        .to_dict()
    )

    uniref_list = (
        df["UniRef"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    url = "https://rest.uniprot.org/uniprotkb/search"
    results_list = []

    for i in range(0, len(uniref_list), batch_size):
        batch = uniref_list[i:i + batch_size]

        print(f"[Batch {i} → {i + len(batch)}]")

        batch_results = {}

        # =========================
        # 🔹 QUERY PRINCIPAL (UniRef)
        # =========================
        query = " OR ".join([f"uniref_cluster_50:{u}" for u in batch])

        params = {
            "query": query,
            "format": "json",
            "size": batch_size
        }

        try:
            r = requests.get(url, params=params)
            r.raise_for_status()

            data = r.json()

            for entry in data.get("results", []):
                accession = entry.get("primaryAccession")

                uni_ref = None
                for ref in entry.get("uniProtKBCrossReferences", []):
                    if ref.get("database") == "UniRef":
                        uni_ref = ref.get("id")
                        break

                if uni_ref and accession:
                    batch_results[uni_ref] = accession

        except Exception as e:
            print(f"  → Batch error: {e}")

        # =========================
        # 🔸 FALLBACK (USANDO OR — FIX REAL)
        # =========================
        missing = [u for u in batch if u not in batch_results]

        if missing:
            fallback_ids = [u.replace("UniRef50_", "") for u in missing]

            print(f"  → Fallback batch ({len(missing)})")

            for j in range(0, len(fallback_ids), fallback_batch_size):
                fb_batch = fallback_ids[j:j + fallback_batch_size]

                # 🔥 AQUI É O FIX
                query_fb = " OR ".join([f"accession:{x}" for x in fb_batch])

                params_fb = {
                    "query": query_fb,
                    "format": "json",
                    "size": len(fb_batch)
                }

                try:
                    r = requests.get(url, params=params_fb)
                    r.raise_for_status()

                    data = r.json()

                    for entry in data.get("results", []):
                        acc = entry.get("primaryAccession")
                        if acc:
                            batch_results[f"UniRef50_{acc}"] = acc

                except Exception as e:
                    print(f"  → Fallback sub-batch error: {e}")

                time.sleep(0.2)

        # =========================
        # 🧱 CONSOLIDA
        # =========================
        for u in batch:
            results_list.append({
                "Gene_Tag": mapping.get(u),
                "UniRef": u,
                "UniProtKB": batch_results.get(u)
            })

        time.sleep(0.3)

    df_out = pd.DataFrame(results_list)
    df_out["UniProtKB"] = df_out["UniProtKB"].astype("object")

    return df_out

def uniprot_info(df_path, batch_size=100):
    df = pd.read_csv(df_path)

    uniprot_list = (
        df["UniProtKB"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    if not uniprot_list:
        print("No UniProtKB entries found")
        return pd.DataFrame(columns=[
            "UniProtKB", "Protein_name", "Gene",
            "Organism", "Length", "Reviewed", "GO_terms"
        ])

    url = "https://rest.uniprot.org/uniprotkb/search"
    all_results = []

    for b in range(0, len(uniprot_list), batch_size):
        batch = uniprot_list[b:b + batch_size]

        print(f"\nProcessing batch {b} -> {b + len(batch)}")

        query = " OR ".join([f"accession:{acc}" for acc in batch])

        params = {
            "query": query,
            "fields": ",".join([
                "accession",
                "protein_name",
                "gene_primary",
                "organism_name",
                "length",
                "reviewed",
                "go_id"
            ]),
            "format": "json",
            "size": batch_size
        }

        batch_results = {}

        try:
            r = requests.get(url, params=params)
            r.raise_for_status()
            data = r.json()

            while True:
                entries = data.get("results") or []

                for entry in entries:
                    acc = entry.get("primaryAccession")

                    protein = (
                        entry.get("proteinDescription", {})
                        .get("recommendedName", {})
                        .get("fullName", {})
                        .get("value")
                    )

                    gene = None
                    genes = entry.get("genes") or []
                    if genes:
                        gene = genes[0].get("geneName", {}).get("value")

                    organism = entry.get("organism", {}).get("scientificName")
                    length = entry.get("sequence", {}).get("length")

                    entry_type = entry.get("entryType", "")
                    reviewed = "reviewed" if "reviewed" in entry_type.lower() else "unreviewed"

                    go_terms_list = []

                    for ref in entry.get("uniProtKBCrossReferences", []):
                        if ref.get("database") == "GO":
                            go_terms_list.append(ref.get("id"))

                    go_terms = ";".join(go_terms_list) if go_terms_list else None

                    batch_results[acc] = {
                        "UniProtKB": acc,
                        "Protein_name": protein,
                        "Gene": gene,
                        "Organism": organism,
                        "Length": length,
                        "Reviewed": reviewed,
                        "GO_terms": go_terms
                    }

                next_link = data.get("next")
                if not next_link:
                    break

                r = requests.get(next_link)
                r.raise_for_status()
                data = r.json()

            for acc in batch:
                if acc in batch_results:
                    all_results.append(batch_results[acc])
                else:
                    all_results.append({
                        "UniProtKB": acc,
                        "Protein_name": None,
                        "Gene": None,
                        "Organism": None,
                        "Length": None,
                        "Reviewed": None,
                        "GO_terms": None
                    })

        except Exception as e:
            print(f"Error in batch {b}: {e}")

            for acc in batch:
                all_results.append({
                    "UniProtKB": acc,
                    "Protein_name": None,
                    "Gene": None,
                    "Organism": None,
                    "Length": None,
                    "Reviewed": None,
                    "GO_terms": None
                })

        time.sleep(0.5)

    df_out = pd.DataFrame(all_results)

    expected_cols = [
        "UniProtKB", "Protein_name", "Gene",
        "Organism", "Length", "Reviewed", "GO_terms"
    ]

    df_out = df_out.reindex(columns=expected_cols)

    if not df_out.empty:
        df_out = df_out.astype("object")

    return df_out

def get_go_terms(df_path, batch_size=200):
    df = pd.read_csv(df_path)
    base_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms"

    print("Analysing GO terms")

    # remove linhas sem GO
    df = df.dropna(subset=["GO_terms"])

    # explode GO terms mantendo relação com UniProt
    go_ids_list_marked = (
        df.assign(GO_split=df["GO_terms"].str.split(";"))
        .explode("GO_split")
        .assign(
            combined=lambda x: x["UniProtKB"].str.strip() + "_" + x["GO_split"].str.strip()
        )["combined"]
        .dropna()
        .tolist()
    )

    # extrai apenas GO IDs
    go_ids_list = [i.split("_")[1] for i in go_ids_list_marked]

    # =========================
    # 🔁 LOOP API (somente isso)
    # =========================
    go_dict = {}

    for i in range(0, len(go_ids_list), batch_size):
        batch = go_ids_list[i:i + batch_size]

        print(f"Analyzing GO terms {i} -> {i + len(batch)}")

        ids_str = ",".join(batch)

        try:
            r = requests.get(f"{base_url}/{ids_str}")
            r.raise_for_status()  # FIX crítico

            data = r.json()

            for result in data.get("results", []):
                go_dict[result.get("id")] = result.get("name")

        except Exception as e:
            print(f"Error in batch {i}: {e}")

        time.sleep(0.2)

    # =========================
    # 🧱 DATAFRAMES (fora do loop)
    # =========================
    df_go = pd.DataFrame(
        list(go_dict.items()),
        columns=["go", "Description"]
    )

    df_pairs = pd.DataFrame({"pair": go_ids_list_marked})
    df_pairs[["marker", "go"]] = df_pairs["pair"].str.split("_", expand=True)

    df_merged = df_pairs.merge(df_go, on="go", how="left")

    df_merged = (
        df_merged.groupby("marker", as_index=False)
        .agg({
            "go": list,
            "Description": list
        })
    )

    # =========================
    # 🔗 MERGE com UniProtKB
    # =========================
    target_df = pd.read_csv("tmp/uniref_to_uniprotkb.csv")

    uniprotkb_go = target_df.merge(
        df_merged,
        left_on="UniProtKB",
        right_on="marker",
        how="left"
    )

    if "marker" in uniprotkb_go.columns:
        uniprotkb_go = uniprotkb_go.drop(columns=["marker"])

    # garante colunas (evita KeyError)
    for col in ["go", "Description"]:
        if col not in uniprotkb_go.columns:
            uniprotkb_go[col] = None

    # salva intermediário
    uniprotkb_go.to_csv("tmp/uniref_to_uniprotkb.csv", index=False)

    # =========================
    # 🧬 MERGE com genes
    # =========================
    genes = pd.read_csv("Data/Entities/genes.csv")

    filt_uniprotkb_go = uniprotkb_go[
        uniprotkb_go["go"].notna() | uniprotkb_go["Description"].notna()
    ].copy()

    genes_go = genes.merge(
        filt_uniprotkb_go[["Gene_Tag", "go", "Description"]],
        on="Gene_Tag",
        how="left"
    )

    # garante colunas antes de selecionar
    for col in ["go", "Description"]:
        if col not in genes_go.columns:
            genes_go[col] = None

    genes_go = genes_go[[
        "Gene_Tag", "Contig", "Start", "Stop",
        "Product", "Databases", "go", "Description",
        "Bin", "Product_Sequence"
    ]]

    genes_go.to_csv("Data/Entities/genes.csv", index=False)

    return genes_go
        
def uniparc_info(df_path, batch_size=100):
    df = pd.read_csv(df_path)

    print("[INFO] Filtering entries without GO terms...")
    uniparc_query = df[df["go"].isna()]

    acc_list = uniparc_query["UniProtKB"].dropna().unique().tolist()
    print(f"[INFO] Total accessions: {len(acc_list)}")

    results = []

    # =========================
    # 🔹 UniProtKB → UniParcId
    # =========================
    for i in range(0, len(acc_list), batch_size):
        batch = acc_list[i:i + batch_size]

        print(f"[UniProtKB] Batch {i} → {i + len(batch)}")

        for j, acc in enumerate(batch, 1):
            print(f"  [UniProtKB] ({j}/{len(batch)}) Fetching {acc}")

            url = f"https://rest.uniprot.org/uniprotkb/{acc}"

            try:
                r = requests.get(url)
                r.raise_for_status()

                data = r.json()

                uniParcId = None
                if "extraAttributes" in data:
                    uniParcId = data["extraAttributes"].get("uniParcId")

                results.append({
                    "UniProtKB": acc,
                    "UniParcId": uniParcId,
                    "EntryType": data.get("entryType"),
                    "Raw": data
                })

                print(f"    → OK | UniParcId: {uniParcId}")

            except Exception as e:
                print(f"    → ERROR: {e}")

                results.append({
                    "UniProtKB": acc,
                    "UniParcId": None,
                    "EntryType": None,
                    "Raw": None
                })

            time.sleep(0.2)

    uniparc_info_raw = pd.DataFrame(results)
    uniparc_info_raw.to_csv("tmp/uniparc_info_raw.csv", index=False)

    print(f"[INFO] UniProtKB step completed: {len(uniparc_info_raw)} rows")

    # =========================
    # 🔹 UniParc
    # =========================
    upis = uniparc_info_raw["UniParcId"].dropna().unique().tolist()
    print(f"[INFO] Total UniParc IDs: {len(upis)}")

    results_upi = []

    for i in range(0, len(upis), batch_size):
        batch = upis[i:i + batch_size]

        print(f"[UniParc] Batch {i} → {i + len(batch)}")

        for j, upi in enumerate(batch, 1):
            print(f"  [UniParc] ({j}/{len(batch)}) Fetching {upi}")

            url = f"https://rest.uniprot.org/uniparc/{upi}"

            try:
                r = requests.get(url)
                r.raise_for_status()

                data = r.json()

                results_upi.append({
                    "UniParcId": upi,
                    "Raw": data
                })

                print("    → OK")

            except Exception as e:
                print(f"    → ERROR: {e}")

                results_upi.append({
                    "UniParcId": upi,
                    "Raw": None
                })

            time.sleep(0.2)

    uniparc_info = pd.DataFrame(results_upi)

    print(f"[INFO] UniParc step completed: {len(uniparc_info)} rows")

    # =========================
    # 🔹 Extract InterPro
    # =========================
    print("[INFO] Extracting InterPro IDs...")

    interpro_list = []

    for idx, raw in enumerate(uniparc_info["Raw"]):
        if pd.isna(raw):
            interpro_list.append(None)
            continue

        if isinstance(raw, str):
            raw = ast.literal_eval(raw)

        ids = []

        for feat in raw.get("sequenceFeatures", []):
            interpro = feat.get("interproGroup", {})
            ip_id = interpro.get("id")

            if ip_id:
                ids.append(ip_id)

        ids = ";".join(sorted(set(ids))) if ids else None
        interpro_list.append(ids)

        if idx % 500 == 0:
            print(f"[DEBUG] Processed {idx} entries")

    uniparc_info["InterPro"] = interpro_list
    uniparc_info = uniparc_info[["UniParcId", "InterPro"]]

    uniparc_info.to_csv("tmp/uniparc_info.csv", index=False)

    print("[INFO] Pipeline completed successfully")

    return uniparc_info
        
def download_interpro2go(url, dir):
    os.makedirs(dir, exist_ok=True)
    
    filename = os.path.basename(url)
    filepath = os.path.join(dir, filename)
    
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(filepath, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

def interpro2go_parser(i2go):
    ipr2go = {}

    with open(i2go) as f:
        for line in f:
            line = line.strip()

            if line.startswith("!"):
                continue

            ipr_match = re.search(r"InterPro:(IPR\d+)", line)
            if not ipr_match:
                continue

            ipr = ipr_match.group(1)

            go_matches = re.findall(r"GO:([^;]+)\s*;\s*(GO:\d+)", line)

            if not go_matches:
                continue

            if ipr not in ipr2go:
                ipr2go[ipr] = []

            for desc, go_id in go_matches:
                ipr2go[ipr].append({
                    "GO_ID": go_id,
                    "Description": desc.strip()
                })
                
    return ipr2go

def interpro2go_convert(df_path, ipr_col, ipr2go):
    df = pd.read_csv(df_path)
    
    rows = []
    
    for _, row in df.iterrows():
        raw_iprs = row[ipr_col]
        
        if pd.isna(raw_iprs):
            continue
        
        iprs = [i.strip() for i in str(raw_iprs).split(";")]
        
        for ipr in iprs:
            go_entries = ipr2go.get(ipr, [])
            
            if not go_entries:
                continue
            
            for go in go_entries:
                rows.append({
                    **row.to_dict(),
                    "IPR": ipr,
                    "GO_ID": go.get("GO_ID"),
                    "GO_Description": go.get("Description")
                })
                
    return pd.DataFrame(rows)

def final_merge_go(df_path):
    df = pd.read_csv(df_path)
    
    df.columns = df.columns.str.strip() 
    
    rows = []
    
    for uniparc, group in df.groupby("UniParcId"):
        iprs = sorted(set(group["IPR"].dropna()))

        go_pairs = set(zip(
            group["GO_ID"].dropna(),
            group["GO_Description"].dropna()
        ))

        go_ids = [g[0] for g in go_pairs]
        go_desc = [g[1] for g in go_pairs]

        rows.append({
            "UniParcId": uniparc,
            "IPR": iprs,
            "GO_ID": go_ids,
            "GO_Description": go_desc
        })
        
    df_uniparc = pd.DataFrame(rows)
    df_uniparc.to_csv("tmp/ipr2go_complete.csv")
    
    df_uniprotkb = pd.read_csv("tmp/uniparc_info_raw.csv")
    
    df_merged = df_uniparc.merge(df_uniprotkb, on="UniParcId", how="left")
    df_merged = df_merged[["UniProtKB", "UniParcId", "IPR", "GO_ID", "GO_Description"]]
    df_merged.rename(columns={"GO_ID": "go", "GO_Description": "Description"}, inplace=True)
    df_merged_subset = df_merged[["UniProtKB", "go", "Description"]]
    
    
    merge_target = pd.read_csv("tmp/uniref_to_uniprotkb.csv")
    merge_target = merge_target[["Gene_Tag", "UniProtKB", "go", "Description"]]
    
    merge_target = merge_target.merge(df_merged_subset, on="UniProtKB", how="left", suffixes=('', '_new'))
    merge_target['go'] = merge_target['go_new'].combine_first(merge_target['go'])
    merge_target['Description'] = merge_target['Description_new'].combine_first(merge_target['Description'])
    merge_target = merge_target.drop(['go_new', 'Description_new'], axis=1)
    
    return merge_target
    

import pandas as pd
import requests
import json
from concurrent.futures import ThreadPoolExecutor
import time
from tqdm import tqdm

def create_database_metadata(genes: pd.DataFrame) -> pd.DataFrame:
    
    """This function defines a core element of the GO module, is extracts the dbs
    from the Genes entity and reorganizes the data according to databases origin.
    Then the function stores the ids in respective columns, returning a dataframe
    to be processed later in the pipeline. 
    
    The reason for this is the data format delivered by Bakta, which mix all the ids
    retrieved from databases in an unique string, the data needs to be saved in lists 
    to be processed easily.
    
    Args:
        genes (pd.DataFrame): Genes entity saved in the Entities module

    Returns:
        pd.DataFrame: Database metadata dataframe
    """
    databases = genes[['Gene_Tag', 'Databases']]
    databases = databases.dropna(subset=['Databases'])
    databases['Databases'] = databases['Databases'].str.split(', ')
    databases['Databases'] = databases['Databases'].apply(
        lambda x: [item for item in x if 'SO' not in item])
    
    for idx, row in databases.iterrows():
        for item in row["Databases"]:
            pref, id = map(str.strip, item.split(":", 1))

            if pref not in databases.columns:
                databases[pref] = None

            current = databases.loc[idx, pref]

            if pd.isna(current):
                databases.loc[idx, pref] = id
            else:
                databases.loc[idx, pref] = f"{current},{id}"
            
    databases['UniRef'] = databases['UniRef'].str.split('_', expand=True)[1]
    
    return databases

def uniref2uniparc(uniref: str, session: requests.Session, retries: int = 3):

    url = f"https://rest.uniprot.org/uniparc/search?query={uniref}"

    for attempt in range(retries):

        try:
            r = session.get(url, timeout=30)
            r.raise_for_status()

            data = r.json()

            results = data.get("results", [])

            if results:

                uniparc_ids = [
                    x["uniParcId"]
                    for x in results
                ]

                return uniref, uniparc_ids

            return uniref, None

        except requests.exceptions.RequestException as e:

            if attempt < retries - 1:
                time.sleep(5)
            else:
                print(f"[ERROR] {uniref}: {e}")
                return uniref, None


def fetch_uniref2uniparc(data_meta: pd.DataFrame, max_workers: int = 5):

    query = (
        data_meta["UniRef"]
        .dropna()
        .unique()
        .tolist()
    )

    print(f"Processing {len(query)} unique UniRef IDs...\n")

    session = requests.Session()

    def worker(uniref):
        return uniref2uniparc(uniref, session)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:

        results = list(
            tqdm(
                executor.map(worker, query),
                total=len(query),
                desc="Mapping UniRef → UniParc",
                unit="IDs"
            )
        )

    session.close()

    mapping = dict(results)

    data_meta["UniParc"] = (
        data_meta["UniRef"]
        .map(mapping)
    )

    mapped = data_meta["UniParc"].notna().sum()

    print("\nFinished!")
    print(f"Mapped {mapped:,} / {len(data_meta):,} rows.")

    return data_meta


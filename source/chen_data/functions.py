import pandas as pd
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
from functools import partial
import time

def convert_json_csv(json_file: Path, output: Path) -> pd.DataFrame:
    
    """Converts a json input into a csv output
    
    Args:
        json_file (Path): .json filepath to be converted
        output (Path): Filepath to contain the converted .csv
        
        Returns:
            The DataFrame object processed by the function
    """
    
    with open(json_file, 'r') as f:
        data = json.load(f)

    records = data['data']
    
    df = pd.DataFrame(records)
    
    df.to_csv(output, index=False)
    
    return df

def filtering(df: pd.DataFrame, out: Path) -> pd.DataFrame:
    
    """Filters a dataframe based on a column content, randomize a number of rows and 
    changes the download links
    
    Args:
        df (pd.DataFrame): Dataframe to be processed
        out (pd.DataFrame): Resulting Dataframe

    Returns:
        Dataframe resulted from the filtering
    """
    total = len(df)
    df = df[df['latitude_and_longitude'] != 'missing']
    new = len(df)
    removed = total - new
    
    print(f"{removed} MAGs without geographical coordinates removed")
    
    print("Randomizing MAGs from metadata")
    
    df = df.sample(n=20, random_state=42)
    
    df['download'] = df['download'].str.replace('https://', 'ftp://', regex=False)
    
    df.to_csv(out, index=False)
    
    return df

def generate_ftp_list(df: pd.DataFrame) -> list:
    
    """Generate ftp download links as a list
    
    Args:
        df (pd.DataFrame): Dataframe to be processed
    
    Returns:
        Download links
    """
    links = df.apply(
        lambda x: f"{x['download']}{x['assembly_id']}/{x['sample_name']}.fa.gz",
        axis=1
    ).tolist()
    
    return links

def get_url(url: str, out: Path) -> None:
    
    """Runs a wget command
    
    Args:
        url (str): url to be downloaded
        out (Path): file destination
        
    Returns:
        None
    """
    
    MAX_RETRIES = 10
    
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            subprocess.run([
                'wget', '-c', '--tries=3', '--waitretry=10', '--read-timeout=60',
                url,
                '-P', str(out)
            ], check=True)
            
            break
        
        except subprocess.CalledProcessError: 
            print(f"Error downloading {url}. Attempt {attempt}/{MAX_RETRIES}")
            
            with open('log/chen_download.log', 'a') as log:
                log.write(
                    f"{time.ctime()} | Attempt {attempt}/{MAX_RETRIES} failed | {url}\n"
                )
                
            time.sleep(30)
        
def fetch_urls(urls: list, out: Path) -> None:
    
    """Runs a parallel processing for a list of urls under a download function
    
    Args:
        urls (list): A list containing the urls
        out (Path): Filepath of the output 
        
    Returns:
        None
    """
    with ThreadPoolExecutor(max_workers=4) as executor:
        list(executor.map(partial(get_url, out=out), urls))
        
   
    

    

    
    
    
    


    

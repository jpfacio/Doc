from goatools.obo_parser import GODag
import subprocess
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from pathlib import Path

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
    
    
                    
            
            
            
        
        
    
                    
                
        
        
    
    
    



















                    
            
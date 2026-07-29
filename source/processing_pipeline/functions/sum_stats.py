import subprocess
from pathlib import Path
import pandas as pd

def seqkit_summary(data_dir: str, output: str) -> pd.DataFrame:
    
    """Run seqkit analysis on all .fa.gz files in a directory
    
    Args:
        data_dir (string): The data filepath containing the sequences
        output (string): Summary files
        
    Returns:
        Summary statistics in dataframe
    """
    
    out_file = Path(output) / "summary_stats.tsv"
    
    subprocess.run(
        f"seqkit stats -a -T -o {out_file} -j 12 {data_dir}/*.fa.gz",
        shell=True,
        check=True
    )
    
    df = pd.read_csv(out_file, sep='\t')
    
    return df
    
def seq_filter(df: pd.DataFrame) -> None:
    
    """Filter a dataframe removing rows where 'num_seqs" values are bigger than 1000
    
    Args:
        df (pd.DataFrame): The data frame to be processed
        
    Returns:
        None
    """
    
    to_remove = df.loc[df['num_seqs'] > 1000, 'file'].to_list()
    
    removed = 0
    
    for file in to_remove:
        path = Path(file)
        if path.exists():
            path.unlink()
            removed += 1
            
    print(f"{removed} files removed")
    
def seq_remove_500(files_dir: Path):
    
    """Check fa.gz files and apply seqkit to remove contigs below 500bp length
    
    Args:
        files_dir (Path): Directory containing .fa.gz files
    Returns:
        None
    """
    
    for fa in files_dir.glob("*.fa.gz"):
        
        original_contigs = int(
            subprocess.check_output(
                ['seqkit', 'stats', '-T', str(fa)],
                text=True
            ).strip().split("\n")[1].split("\t")[3]
        )
        
        tmp = fa.with_suffix(".tmp.fa.gz")
        
        with open(tmp, 'wb') as f:
            subprocess.run(
                ['seqkit', 'seq', '-m', '500', str(fa)],
                stdout=f,
                stderr=subprocess.DEVNULL,
                check=True
            )
            
        processed_contigs = int(
            subprocess.check_output(
                ['seqkit', 'stats', '-T', str(tmp)],
                text=True
            ).strip().split("\n")[1].split("\t")[3]
        )
        
        removed = original_contigs - processed_contigs
        
        Path.replace(tmp, fa)
        
        print(f'{fa.name}: {removed} contigs from {original_contigs}.')
    
    


    

    

    
    
    

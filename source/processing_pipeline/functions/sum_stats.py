import subprocess
from pathlib import Path
import pandas as pd
import time

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
    
    """Check fa.gz files and apply seqkit to remove contigs below 500bp length, uncompressing the file
    in the end.
    
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
        
        fa.rename(fa.with_suffix(''))
        
        print(f'{fa.name}: {removed} contigs from {original_contigs}.')

def run_checkm(files_dir: Path, output_dir: Path, log: Path, threads: int = 8) -> pd.DataFrame:
    
    """Runs CheckM main analysis and generates its .tsv summary of the results
    
    Args:
        files_dir (Path): Path containing the files to be processed
        output_dir (Path): Directory where the CheckM output will be sent
        log (Path): Path of the log file
        threads (int): Number of threads to be passed to CheckM command

    Returns:
        pd.DataFrame: Dataframe containing the CheckM results
    """
    output_dir = output_dir / 'checkm'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    log_file = log / 'checkm.log'
    
    print("Starting CheckM Analysis")
    
    with open(log_file, 'w') as log:
        subprocess.run(
            [
                'checkm', 'lineage_wf', '-t',
                str(threads),
                '-x', 'fa',
                str(files_dir),
                str(output_dir)
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
            text=True
        )
        
    final_table = "tmp/checkm_results.tsv"
    
    with open(final_table, 'w') as f:
        subprocess.run(
            [
                "checkm", "qa", "--tab_table", "-o", "2",
                str(output_dir / "lineage.ms"),
                str(output_dir)
            ],
            stdout=f,
            check=True,
            text=True
        )
    
    with open(final_table, "r") as f:
        lines = [line for line in f if not line.startswith("[")]
        
    with open(final_table, "w") as f:
        f.writelines(lines)

    print("Done!")
    
    df = pd.read_csv(final_table, sep="\t")
    
    return df 

def final_processing(seqkit: pd.DataFrame, checkm: pd.DataFrame, tmp_file: Path) -> None:
    
    
    checkm_dir = tmp_file / "checkm"
    
    for _ in range(10):
        subprocess.run(["rm", "-rf", str(checkm_dir)],
                       stdout = subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        
        if not checkm_dir.exists():
            break
        
        time.sleep(1)
    
    subprocess.run(["rm", "-rf", str(checkm_dir)], check=True)
    
    seqkit = seqkit[
        [
            "file",
            "num_seqs",
            "N50", 
            "GC(%)"
        ]
    ]
    
    seqkit.loc[:, 'file'] = seqkit['file'].str.removesuffix('.gz')
    seqkit = seqkit.loc[seqkit['num_seqs'] <= 1000]
    
    checkm = checkm[
        [
            "Completeness",
            "Contamination",
            "Strain heterogeneity"
        ]
    ]
    
    df = pd.concat([seqkit.reset_index(drop=True), 
                    checkm.reset_index(drop=True)], 
                   axis=1)
    
    df['MIMAG'] = 'Low'
    
    df.loc[
        (df["Completeness"] >= 50) &
        (df["Contamination"] <= 10), "MIMAG"
    ] = "Medium"
    
    df.loc[
        (df["Completeness"] >= 90) &
        (df["Contamination"] <= 5), "MIMAG"
    ] = "High"
    
    low = df.loc[df['MIMAG'] == 'Low', 'file'].to_list()
    
    for file in low:
        subprocess.run(['rm', '-f', file],
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        
    df = df.loc[df['MIMAG'] != "Low"]
        
    print(f"{len(low)} low quality MAGs removed.")
    
    table_path = tmp_file / "qc_metrics.tsv"
    
    df.to_csv(table_path, sep="\t", index=False)
    
    subprocess.run(["rm", "tmp/summary_stats.tsv", "tmp/checkm_results.tsv"], check=True)
    


    

    
    
    

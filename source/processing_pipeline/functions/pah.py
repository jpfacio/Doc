import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
from tqdm import tqdm


DEFAULT_OUTFMT = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp",
]


def build_query_faa(genes_csv: Path, out_faa: Path) -> int:

    """Reads the Genes entity and writes the query protein fasta used for the
    DIAMOND screen. Genes without a translation are dropped, so every record
    written holds a real sequence.

    Args:
        genes_csv (Path): Path to the Genes entity CSV
        out_faa (Path): Output fasta path

    Returns:
        int: Number of sequences written
    """
    genes = pd.read_csv(genes_csv)
    sequences = genes[["Gene_Tag", "Product_Sequence"]].dropna()

    with open(out_faa, "w") as fh:
        for tag, seq in sequences.itertuples(index=False):
            fh.write(f">{tag}\n{seq}\n")

    return len(sequences)


def prepare_dbs(base_dir: Path, force: bool = False) -> list[str]:

    """Discovers the PAH reference families and makes sure every one has a
    DIAMOND database. Databases are committed to the repo, so makedb is only
    triggered for missing files or when force=True (e.g. after editing a .faa).

    Args:
        base_dir (Path): Directory holding the .faa/.dmnd reference files
        force (bool): Rebuild every database, even if it already exists

    Returns:
        list: Family names (reference fasta stems)
    """
    base_dir = Path(base_dir)
    families = []

    for path in sorted(base_dir.glob("*.faa")):
        family = path.stem
        families.append(family)
        dmnd = base_dir / f"{family}.dmnd"
        if force or not dmnd.exists():
            subprocess.run(
                ["diamond", "makedb", "--in", str(path),
                 "--db", str(base_dir / family)],
                check=True,
            )

    return families


def screen_one(query_faa: Path, family: str, base_dir: Path, out_dir: Path,
               evalue: float = 1e-5, threads: int = 2) -> Path:

    """Runs diamond blastp for a single PAH family against the query fasta.
    Sensitive mode is used because the curated reference sets are small;
    the best hit per query is kept. Stderr is captured to its own log file.

    Args:
        query_faa (Path): Query protein fasta
        family (str): PAH family name
        base_dir (Path): Directory holding the reference .dmnd files
        out_dir (Path): Directory for the result and log files
        evalue (float): Maximum e-value for reported hits
        threads (int): DIAMOND threads for this process

    Returns:
        Path: Result tsv path
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    out = out_dir / f"results_{family}.tsv"

    result = subprocess.run(
        [
            "diamond", "blastp",
            "--query", str(query_faa),
            "--db", str(Path(base_dir) / f"{family}.dmnd"),
            "--out", str(out),
            "--outfmt", "6", *DEFAULT_OUTFMT,
            "--evalue", str(evalue),
            "--max-target-seqs", "1",
            "--sensitive",
            "--threads", str(threads),
            "--quiet",
        ],
        capture_output=True,
        text=True,
    )

    with open(out_dir / f"stderr_{family}.log", "w") as fh:
        fh.write(result.stderr)

    result.check_returncode()

    return out


def aggregate(out_dir: Path, families: list[str], genes_csv: Path,
              min_identity: float | None = None,
              min_qcov: float | None = None) -> pd.DataFrame:

    """Collects the per-family BLAST results into a single candidates frame
    and joins gene metadata (Bin, Product). Optional post-filters on identity
    and query coverage.

    Args:
        out_dir (Path): Directory holding the results_{family}.tsv files
        families (list): PAH family names
        genes_csv (Path): Path to the Genes entity CSV
        min_identity (float | None): Keep only hits at/above this %% identity
        min_qcov (float | None): Keep only hits at/above this query coverage

    Returns:
        pd.DataFrame: Candidates frame
    """
    out_dir = Path(out_dir)
    frames = []

    for family in families:
        out = out_dir / f"results_{family}.tsv"
        if not out.exists():
            continue
        df = pd.read_csv(out, sep="\t", header=None, names=DEFAULT_OUTFMT)
        df.insert(0, "Family", family)
        frames.append(df)

    if not frames:
        raise FileNotFoundError(
            f"No PAH screening results found in {out_dir}")

    hits = pd.concat(frames, ignore_index=True)
    hits = hits.rename(columns={"qseqid": "Gene_Tag", "sseqid": "Subject"})

    genes = pd.read_csv(genes_csv)[["Gene_Tag", "Bin", "Product"]]
    candidates = hits.merge(genes, on="Gene_Tag", how="left")

    cols = {
        "pident": "Identity",
        "evalue": "Evalue",
        "bitscore": "Bitscore",
        "qcovhsp": "Qcov",
    }
    candidates = candidates.rename(columns=cols)[
        ["Gene_Tag", "Bin", "Product", "Family", "Subject",
         "Identity", "Evalue", "Bitscore", "Qcov"]
    ]

    if min_identity is not None:
        candidates = candidates[candidates["Identity"] >= min_identity]
    if min_qcov is not None:
        candidates = candidates[candidates["Qcov"] >= min_qcov]

    return candidates.reset_index(drop=True)


def summarize(candidates: pd.DataFrame) -> pd.DataFrame:

    """Per-bin candidate counts for each PAH family.

    Args:
        candidates (pd.DataFrame): Candidates frame

    Returns:
        pd.DataFrame: Bin x Family count table
    """
    summary = (
        candidates.groupby(["Bin", "Family"])
        .size()
        .reset_index(name="N")
    )
    return summary.sort_values(["Bin", "Family"]).reset_index(drop=True)


def run(base_dir: str = "support_files/pah_db_v2.4",
        genes_csv: str = "Data/Entities/genes.csv",
        query_faa: str = "tmp/sequences.faa",
        out_dir: str = "log/pah_log",
        out_ent: str = "Data/Entities/pah.csv",
        out_report: str = "Data/Reports/pah_summary.csv",
        evalue: float = 1e-5,
        min_identity: float | None = None,
        min_qcov: float | None = None,
        max_workers: int = 4,
        threads: int = 2,
        force_dbs: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:

    """Full PAH-degradation gene screen: build the query fasta from the Genes
    entity, screen every family with DIAMOND in parallel, then aggregate and
    summarize the hits into the PAH entity and report.

    Args:
        base_dir (str): Directory holding the .faa/.dmnd reference files
        genes_csv (str): Path to the Genes entity CSV
        query_faa (str): Path where the query fasta is written
        out_dir (str): Directory for raw BLAST and stderr logs
        out_ent (str): Path of the PAH entity CSV
        out_report (str): Path of the per-bin summary CSV
        evalue (float): Maximum e-value for reported hits
        min_identity (float | None): Optional identity post-filter
        min_qcov (float | None): Optional query coverage post-filter
        max_workers (int): Families screened in parallel
        threads (int): DIAMOND threads per family process
        force_dbs (bool): Rebuild DIAMOND databases even if present

    Returns:
        tuple: (candidates frame, per-bin summary frame)
    """
    if shutil.which("diamond") is None:
        raise RuntimeError(
            "diamond not found on PATH - activate the MbXDb-env")

    n_seqs = build_query_faa(Path(genes_csv), Path(query_faa))
    print(f"Query fasta: {n_seqs:,} sequences -> {query_faa}")

    families = prepare_dbs(Path(base_dir), force=force_dbs)
    print(f"PAH databases ready: {len(families)} families")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def worker(family):
        screen_one(Path(query_faa), family, Path(base_dir), out_dir,
                   evalue=evalue, threads=threads)
        return family

    failed = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(worker, family): family for family in families}
        for future in tqdm(as_completed(futures), total=len(families),
                           desc="Screening against PAH DBs", unit="DB"):
            family = futures[future]
            try:
                future.result()
            except Exception as e:
                print(f"[ERROR] {family}: {e}")
                failed.append(family)

    if failed:
        raise RuntimeError(
            f"PAH screening failed for families: {', '.join(failed)}")

    candidates = aggregate(out_dir, families, Path(genes_csv),
                           min_identity=min_identity, min_qcov=min_qcov)
    candidates.to_csv(out_ent, index=False)

    summary = summarize(candidates)
    summary.to_csv(out_report, index=False)

    print(f"PAH candidates: {len(candidates):,} hits across "
          f"{summary['Bin'].nunique()} bins")
    print(f"Entity: {out_ent}")
    print(f"Summary: {out_report}")

    return candidates, summary

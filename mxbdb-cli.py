"""mxbdb-cli - query the MbXDb entity CSVs.

Standalone, pandas-only query CLI over Data/Entities/{genes,pah,bins}.csv.
Run from the repo root. Does not import the pipeline `functions` package
(that would pull in the full toolchain).
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

ENTITY_DIR = Path("Data/Entities")
GENES_CSV = ENTITY_DIR / "genes.csv"
PAH_CSV = ENTITY_DIR / "pah.csv"
BINS_CSV = ENTITY_DIR / "bins.csv"

pd.set_option("display.max_colwidth", 60)
pd.set_option("display.width", 140)


def load_genes() -> pd.DataFrame:
    return pd.read_csv(GENES_CSV)


def load_pah() -> pd.DataFrame:
    return pd.read_csv(PAH_CSV)


def load_bins() -> pd.DataFrame:
    return pd.read_csv(BINS_CSV)


def emit(df: pd.DataFrame, args: argparse.Namespace) -> None:
    if args.json:
        sys.stdout.write(df.to_json(orient="records", indent=2) + "\n")
    elif args.csv:
        df.to_csv(sys.stdout, index=False)
    else:
        sys.stdout.write(df.to_string(index=False) + "\n")


def add_io_args(sp: argparse.ArgumentParser) -> None:
    sp.add_argument("--csv", action="store_true",
                    help="print as CSV instead of a pretty table")
    sp.add_argument("--json", action="store_true",
                    help="print as JSON records instead of a pretty table")


def cmd_genes_search(args: argparse.Namespace) -> None:
    genes = load_genes()
    mask = genes["Product"].fillna("").str.contains(
        args.keyword, case=False, regex=False)
    mask |= genes["Gene_Tag"].str.contains(args.keyword, case=False, regex=False)
    if args.bin:
        mask &= genes["Bin"] == args.bin
    result = genes.loc[mask, ["Gene_Tag", "Bin", "Product", "Databases"]]
    if args.limit:
        result = result.head(args.limit)
    emit(result, args)


def cmd_genes_get(args: argparse.Namespace) -> None:
    genes = load_genes()
    row = genes[genes["Gene_Tag"] == args.tag]
    if row.empty:
        sys.exit(f"error: no gene with tag '{args.tag}'")
    emit(row, args)
    hits = load_pah()
    hits = hits[hits["Gene_Tag"] == args.tag]
    print("\nPAH hits:")
    emit(hits, args)


def cmd_pah_top(args: argparse.Namespace) -> None:
    pah = load_pah()
    no_filters = not any([args.family, args.bin,
                          args.min_identity is not None,
                          args.min_qcov is not None])
    if no_filters and args.top is None:
        sys.exit("error: 'pah top' needs at least one filter "
                 "(--family/--bin/--min-identity/--min-qcov) or -n/--top")
    if args.family:
        pah = pah[pah["Family"] == args.family]
    if args.bin:
        pah = pah[pah["Bin"] == args.bin]
    if args.min_identity is not None:
        pah = pah[pah["Identity"] >= args.min_identity]
    if args.min_qcov is not None:
        pah = pah[pah["Qcov"] >= args.min_qcov]
    result = pah.sort_values("Bitscore", ascending=False)
    if args.top:
        result = result.head(args.top)
    emit(result, args)


def cmd_pah_summary(args: argparse.Namespace) -> None:
    pah = load_pah()
    if args.family:
        pah = pah[pah["Family"] == args.family]
    if args.bin:
        pah = pah[pah["Bin"] == args.bin]
    summary = pah.groupby(["Bin", "Family"]).size().reset_index(name="N")
    summary = summary.sort_values(["Bin", "Family"]).reset_index(drop=True)
    emit(summary, args)


def cmd_bins(args: argparse.Namespace) -> None:
    bins = load_bins()
    genes = load_genes()
    pah = load_pah()
    gene_counts = genes.groupby("Bin").size().reset_index(name="Gene_Count")
    pah_hits = pah.groupby("Bin").size().reset_index(name="Pah_Hits")
    pah_fam = pah.groupby("Bin")["Family"].nunique().reset_index(
        name="Pah_Families")
    result = bins.merge(gene_counts, on="Bin", how="left")
    result = result.merge(pah_hits, on="Bin", how="left")
    result = result.merge(pah_fam, on="Bin", how="left")
    for col in ("Gene_Count", "Pah_Hits", "Pah_Families"):
        result[col] = result[col].fillna(0).astype(int)
    emit(result[["Bin", "Sample", "Study_ID",
                 "Gene_Count", "Pah_Hits", "Pah_Families"]], args)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="mxbdb-cli",
        description="Query the MbXDb entity CSVs (genes, pah, bins). "
                    "Run from the repo root.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_genes = sub.add_parser("genes", help="query the genes entity")
    gsub = p_genes.add_subparsers(dest="genes_cmd", required=True)
    g_search = gsub.add_parser("search", help="search genes by product/tag keyword")
    g_search.add_argument("keyword", help="substring to match against Product or Gene_Tag")
    g_search.add_argument("--bin", help="restrict to a single bin")
    g_search.add_argument("--limit", type=int, help="show at most N rows")
    add_io_args(g_search)
    g_search.set_defaults(func=cmd_genes_search)

    g_get = gsub.add_parser("get", help="show a gene's full record and its PAH hits")
    g_get.add_argument("tag", help="Gene_Tag, e.g. GOMC.bin.13102_1547")
    add_io_args(g_get)
    g_get.set_defaults(func=cmd_genes_get)

    p_pah = sub.add_parser("pah", help="query the pah entity")
    psub = p_pah.add_subparsers(dest="pah_cmd", required=True)
    p_top = psub.add_parser("top", help="list PAH hits sorted by bitscore")
    p_top.add_argument("--family", help="PAH family (e.g. nahAc, nidA, phdF)")
    p_top.add_argument("--bin", help="restrict to a single bin")
    p_top.add_argument("--min-identity", type=float, help="minimum %% identity")
    p_top.add_argument("--min-qcov", type=float, help="minimum query coverage")
    p_top.add_argument("-n", "--top", type=int, help="show only the N best hits")
    add_io_args(p_top)
    p_top.set_defaults(func=cmd_pah_top)

    p_sum = psub.add_parser("summary", help="per-bin x family hit counts")
    p_sum.add_argument("--family", help="restrict to a single family")
    p_sum.add_argument("--bin", help="restrict to a single bin")
    add_io_args(p_sum)
    p_sum.set_defaults(func=cmd_pah_summary)

    p_bins = sub.add_parser("bins", help="overview: one row per bin")
    add_io_args(p_bins)
    p_bins.set_defaults(func=cmd_bins)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()

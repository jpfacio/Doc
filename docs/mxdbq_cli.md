# mxbdb-cli.py — Entity Query CLI

> **Module:** `mxbdb-cli.py`
> **Purpose:** a standalone, pandas-only CLI for querying the MbXDb entity CSVs
> (`Data/Entities/{genes,pah,bins}.csv`).
> **Requirements:** only `pandas`. Does **not** import the pipeline `functions`
> package (that would drag in habanero/goatools/biopython).
> **CWD rule:** like every other script in the repo, run it from the repo root —
> all entity paths are CWD-relative.

```
python mxbdb-cli.py <command> <subcommand> [options]
```

All commands accept `--csv` and `--json` to switch from the default pretty
table to machine-readable output.

---

## `genes`

### `genes search <keyword>`

Substring match (case-insensitive) against `Product` **or** `Gene_Tag`.
Columns: `Gene_Tag, Bin, Product, Databases`.

```
python mxbdb-cli.py genes search dioxygenase
python mxbdb-cli.py genes search dioxygenase --bin GOMC.bin.13102 --limit 10
```

| Option | Description |
|---|---|
| `keyword` | substring to match (positional) |
| `--bin BIN` | restrict to a single bin |
| `--limit N` | show at most N rows |
| `--csv` / `--json` | output format |

### `genes get <tag>`

Full record for one `Gene_Tag` (all genes columns, including
`Product_Sequence`), followed by that gene's PAH hits from `pah.csv`.

```
python mxbdb-cli.py genes get GOMC.bin.13102_1547
```

This is the only subcommand that joins two entities. A gene with no PAH hits
prints an empty `PAH hits:` block — not an error. A nonexistent tag exits with
`error: no gene with tag '<tag>'` (exit code 1).

---

## `pah`

### `pah top`

PAH hits (DIAMOND blastp results) filtered and sorted by `Bitscore`
descending. Columns: `Gene_Tag, Bin, Product, Family, Subject, Identity,
Evalue, Bitscore, Qcov`.

```
python mxbdb-cli.py pah top --family nidA
python mxbdb-cli.py pah top --min-identity 35 --min-qcov 80 -n 20
python mxbdb-cli.py pah top --bin GOMC.bin.9704 --family phdF
```

| Option | Description |
|---|---|
| `--family FAMILY` | PAH family, e.g. `nahAc`, `nidA`, `phdF` |
| `--bin BIN` | restrict to a single bin |
| `--min-identity X` | keep only hits with identity >= X (%) |
| `--min-qcov X` | keep only hits with query coverage >= X (%) |
| `-n, --top N` | show only the N best hits |
| `--csv` / `--json` | output format |

> Requires at least one filter or `-n`; running bare exits with a usage error
> (exit code 1). The `--min-identity`/`--min-qcov` flags apply the same
> thresholds that `pah.py` exposes but leaves off by default.

### `pah summary`

Per-bin × family hit counts, computed live from `pah.csv` (mirrors
`summarize()` in `pah.py`). Columns: `Bin, Family, N`. Only bins with hits
appear.

```
python mxbdb-cli.py pah summary
python mxbdb-cli.py pah summary --family phdF
python mxbdb-cli.py pah summary --bin GOMC.bin.13102
```

| Option | Description |
|---|---|
| `--family FAMILY` | restrict to a single family |
| `--bin BIN` | restrict to a single bin |
| `--csv` / `--json` | output format |

---

## `bins`

One row per bin: `Bin, Sample, Study_ID, Gene_Count, Pah_Hits,
Pah_Families`. Joins `bins.csv` with per-bin counts from `genes.csv` and
`pah.csv`.

```
python mxbdb-cli.py bins
python mxbdb-cli.py bins --csv
```

| Option | Description |
|---|---|
| `--csv` / `--json` | output format |

No filters in the first version — the table is only 20 rows.

---

## Output formats

| Flag | Behavior |
|---|---|
| *(none)* | pretty table via `df.to_string(index=False)` |
| `--csv` | CSV via `df.to_csv(sys.stdout, index=False)` |
| `--json` | JSON array of records, 2-space indented |

`--csv` and `--json` are mutually exclusive by convention; if both are given,
`--json` wins (it is checked first).

---

## Notes

- Entity files are loaded lazily per subcommand (genes ~45k rows, pah ~922,
  bins 20), so startup stays fast.
- `genes search` treats the keyword as a literal substring (no regex).
- Verification = execution (no test framework in this repo); smoke-check with
  the examples above from the repo root.

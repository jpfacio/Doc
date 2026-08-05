# MbXDb – Metagenomic biorremediation Xenobiotics Database

v0.0.1

MbXDb is a modular bioinformatics pipeline designed to process metagenome-assembled genomes (MAGs/bins) through quality filtering, functional annotation, and downstream analysis. The pipeline integrates state-of-the-art tools to ensure high-quality, biologically meaningful results.

> **Note**: This project is under active development. Taxonomic and xenobiotic analysis modules are upcoming.

---

## Pipeline Overview

The MXD pipeline executes the following steps sequentially:

1. **Filtering**: Retains only bins with **<1,000 contigs** (as a proxy for complexity/completeness).
2. **Quality Control**: Assesses bin quality using **CheckM2** (completeness, contamination, strain heterogeneity).
3. **Functional Annotation**: Annotates genes and proteins using **Bakta**, providing:
   - COG, Pfam, TIGRFAM annotations
   - rRNA, tRNA, ncRNA detection
   - Signal peptides and transmembrane helices
   - Taxonomic hints via conserved marker genes
4. **(In Progress)** Taxonomic profiling & xenobiotic metabolism potential analysis.

---

## Directory Structure

```
MXD/
├── Data/ # All input and output data
│ ├── Entities/ # Metadata, sample info, entities (reserved)
│ ├── Raw/ # Raw input data
│ │ ├── Bins/ # Input bins
│ │ └── Processed/ # Intermediate filtered bins
│ └── Reports/ # Final analysis reports
│
├── source/ # Analysis code modules
│ ├── chen_data/ # Scripts for Chen et al. (2022) database processing
│ ├── processing_pipeline/ # Core pipeline logic 
│ └── test/ # Test subset and validation scripts
│
├── mxbdb-build.sh # Main entrypoint script 
├── mxbdb-cli.py # Query CLI (Data/Entities CSVs)
└── environment.yml # Conda environment 
```

## Quick Start

### Prerequisites
- Linux
- `conda` (or `mamba`)
- At least 16 GB RAM 

### 1. Clone and set up
```bash
git clone https://github.com/jpfacio/MXD.git
cd MXD
conda env create -f environment.yml
conda activate mxd
```

### 2. Run the pipeline

-Run full pipeline on all bins
`./mxbdb-build.sh --all`

-Or run only on test subset
`./mxbdb-build.sh --t`

-Or process Chen et al. (2022) database + your bins
`./mxbdb-build.sh --chen`

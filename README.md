# Transcriptomic-Based Drug Repurposing Pipeline

A command-line pipeline for identifying drug repurposing candidates from public RNA-seq data. Starting from a GEO accession number and a disease name, the pipeline downloads expression data, identifies differentially expressed genes, builds a disease transcriptomic signature, and queries the LINCS L1000 database for drugs that reverse that signature.

```
GEO accession
     │
     ▼
[Tool 1] Download count matrix + sample metadata (GEOparse)
     │
     ▼
[Tool 1b] Interactive sample group assignment (disease / control / exclude)
     │
     ▼
[Tool 2] Differential gene expression — DESeq2 (PyDESeq2)
     │
     ▼
[Tool 3] Disease transcriptomic signature (MyGeneInfo)
     │
     ▼
[Tool 4] Drug repurposing candidates (L1000FWD / Enrichr fallback)
     │
     ▼
Provenance log  +  Markdown report
```

---
## Warning

This code is my own personal project. Therefore you'll likely find bugs or unexpected behavior that might not be documented here. Please let me kow if you find anything and feel free to contribute to it as 
you wish. I have tested it with a couple of rare disease examples, but since the GEO datasets are not standardized and are very messy, be especially carefil with the sample assignment parts of the pipeline. 
Right now it only does single factor differential gene expression anaysis, control vs disease. 

---
## Output

All outputs for an accession (e.g. `GSE297335`) are written to `data/GSE297335/`:

| File | Description |
|------|-------------|
| `GSE297335_counts.csv` | Raw count matrix (genes × samples) |
| `GSE297335_metadata.csv` | GEO sample metadata |
| `GSE297335_groups.csv` | User-assigned disease/control labels |
| `GSE297335_deseq2_results.csv` | Full DESeq2 results table |
| `GSE297335_signature.json` | Disease signature (up/down gene lists) |
| `GSE297335_signature.csv` | Signature in human-readable CSV form |
| `GSE297335_drug_candidates.json` | Ranked drug candidates |
| `GSE297335_drug_candidates.csv` | Ranked drug candidates (CSV) |

Reports are written to `reports/`:

| File | Description |
|------|-------------|
| `{Disease}_provenance.json` | Full audit trail (parameters, warnings, file paths) |
| `{Disease}_report.md` | Human-readable Markdown report |

---

## Setup

### 1. Python version

Python **3.10 or later** is required (the code uses union-type hints such as `str | None`).

```bash
python3 --version
```

### 2. Create a virtual environment

```
python3 -m venv .venv
source .venv/bin/activate        # macOS / Linux
# .venv\Scripts\activate         # Windows
```
or

```
conda create -n tdre --file requirements.txt
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Configure API credentials

Create a `.env` file in the project root (this file is gitignored and never committed):

```bash
cp .env.example .env   # if an example file exists, otherwise create from scratch
```

Edit `.env` and fill in your credentials:

```dotenv
NCBI_API_KEY=your_ncbi_api_key_here
EMAIL_ADDRESS=your_email@example.com
ORPHANET_API_KEY=your_orphanet_key_here   # optional
```

**Getting an NCBI API key** (recommended — raises rate limits from 3 to 10 requests/second):

1. Create a free account at [https://www.ncbi.nlm.nih.gov/account/](https://www.ncbi.nlm.nih.gov/account/)
2. Go to **Account Settings → API Key Management**
3. Generate a key and paste it into `.env`

The pipeline will work without an NCBI API key, but GEO downloads may be throttled.

---

## Running the pipeline

### Interactive mode (recommended for first use)

```bash
python run_pipeline.py
```

You will be prompted for a GEO accession number and a disease name, then walked through sample group assignment interactively.

### Command-line arguments

```bash
python run_pipeline.py --accession GSE297335 --disease "Barth Syndrome"
# or short flags:
python run_pipeline.py -a GSE297335 -d "Barth Syndrome"
```

### Sample group assignment

This is the most important interactive step. GEO metadata is inconsistent, so the pipeline infers suggested groups from sample characteristics but always asks for confirmation.

Before the per-sample loop you will see:
- A summary of all characteristics fields present in the dataset (field names + unique values + counts), so you can understand the metadata structure at a glance
- A **batch assignment** prompt where you can map `field: value → group` for entire subsets at once (useful for clean datasets — just press Enter to skip to per-sample mode)

For each sample you choose:

| Key | Meaning |
|-----|---------|
| `d` | disease — include as a disease sample |
| `c` | control — include as a healthy/control sample |
| `e` | exclude — omit from analysis (e.g. treated samples, outliers) |
| Enter | accept the suggested assignment |

Group assignments are saved to `data/{accession}/{accession}_groups.csv`. If you re-run Tool 1 for the same accession, you will be offered the option to reuse the saved assignments.

DESeq2 requires **at least 3 samples in each group**.

---

## Running tools individually

Each tool can be run as a standalone script. Pass `--help` to see all options.

**Required for every tool:** `-a / --accession` and `-d / --disease`. All file paths default to the standard `data/{accession}/` layout and can be overridden explicitly.

```bash
# Tool 1 — download GEO data + interactive sample group assignment
python tools/1_geo_download.py -a GSE297335 -d "Barth Syndrome"

# Tool 2 — DESeq2 differential expression
python tools/2_deg_analysis.py -a GSE297335 -d "Barth Syndrome"
# with relaxed thresholds:
python tools/2_deg_analysis.py -a GSE297335 -d "Barth Syndrome" --padj 0.1 --lfc 0.5

# Tool 3 — build disease transcriptomic signature
python tools/3_signature.py -a GSE297335 -d "Barth Syndrome"
# with custom signature size:
python tools/3_signature.py -a GSE297335 -d "Barth Syndrome" --max-genes 200

# Tool 4 — query L1000FWD for drug candidates
python tools/4_drug_repurposing.py -a GSE297335 -d "Barth Syndrome"
```

**Explicit path overrides** (useful if files are in non-standard locations):

```bash
python tools/2_deg_analysis.py -a GSE297335 -d "Barth Syndrome" \
    --counts /path/to/counts.csv \
    --groups /path/to/groups.csv \
    --data-dir /path/to/output/

python tools/4_drug_repurposing.py -a GSE297335 -d "Barth Syndrome" \
    --signature /path/to/signature.json \
    --data-dir /path/to/output/
```

---

## External services

The pipeline calls the following external APIs (all free, no authentication required beyond NCBI):

| Service | Used for | Fallback |
|---------|----------|---------|
| [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/) | Downloading expression data | — |
| [MyGeneInfo](https://mygene.info/) | Gene ID → symbol mapping | — |
| [L1000FWD](https://maayanlab.cloud/l1000fwd/) | Drug connectivity scoring | Enrichr |
| [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | BRD ID → drug name lookup | Raw BRD ID |
| [Enrichr](https://maayanlab.cloud/Enrichr/) | Fallback if L1000FWD is down | — |

---

## Local drug annotation resources

The `resources/` directory contains two files from the [Broad Institute Drug Repurposing Hub](https://www.broadinstitute.org/drug-repurposing-hub) that are bundled with the repository and used to annotate drug candidates with mechanism of action, protein target, and clinical phase:

| File | Contents |
|------|----------|
| `repo-sample-annotation-*.txt` | BRD compound ID → drug name mapping |
| `repo-drug-annotation-*.txt` | Drug name → MoA, target, clinical phase, indication |

---

## Report

At the end of a successful run, a Markdown report is written to `reports/{Disease}_report.md`. It includes:

- Dataset metadata (GEO accession, sample counts, organism, tissue)
- Sample group assignment summary
- DEG analysis statistics (genes tested, significant up/down, thresholds used)
- Disease signature composition (top up/down genes)
- Ranked drug candidate table (score, confidence, clinical phase, MoA, target)
- Link to interactive L1000FWD results (if available)
- Quality control warnings raised during the run
- Methods section

A machine-readable `{Disease}_provenance.json` log captures every parameter, output file path, and warning for full reproducibility.

---

## Project structure

```
tdre/
├── run_pipeline.py          # Main orchestrator — start here
├── requirements.txt
├── .env                     # API credentials (not committed)
│
├── tools/
│   ├── 1_geo_download.py    # GEO download + interactive sample assignment
│   ├── 2_deg_analysis.py    # DESeq2 differential expression
│   ├── 3_signature.py       # Disease signature construction
│   └── 4_drug_repurposing.py # L1000FWD drug connectivity query
│
├── utils/
│   └── provenance.py        # ProvenanceLog class + report generation
│
├── resources/               # Broad Drug Repurposing Hub annotation files
│   ├── repo-sample-annotation-*.txt
│   └── repo-drug-annotation-*.txt
│
├── data/                    # Created at runtime, one subdirectory per accession
│   └── GSE297335/
│       ├── GSE297335_counts.csv
│       └── ...
│
└── reports/                 # Provenance logs and Markdown reports
    ├── Barth_Syndrome_provenance.json
    └── Barth_Syndrome_report.md
```

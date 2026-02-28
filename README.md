# TDRE — Transcriptomic Drug Repurposing Engine

TDRE is a bioinformatics pipeline that generates drug repurposing hypotheses for rare and orphan diseases using publicly available transcriptomic data. Given a disease name, it automatically finds bulk RNA-seq datasets on NCBI GEO, runs differential expression analysis, builds a disease transcriptomic signature, and identifies existing drugs whose gene expression effects reverse that signature.

Output is a ranked list of drug candidates with confidence scores, mechanism of action, and molecular targets — framed as hypotheses for experimental follow-up, not validated targets.

---

## How It Works

```
"Disease Name"
      ↓
1. Disease lookup      — NCBI MeSH + Orphanet + PubMed feasibility check
      ↓
2. GEO dataset search  — Find and score candidate RNA-seq datasets
      ↓
3. Dataset selection   — Human reviews and selects best dataset [interactive]
      ↓
4. Data download       — Download count matrix and sample metadata from GEO
      ↓
5. Sample assignment   — Human assigns samples to disease/control groups [interactive]
      ↓
6. DEG analysis        — Differential expression via DESeq2 (PyDESeq2)
      ↓
7. Signature building  — Filter significant DEGs, map to HGNC gene symbols
      ↓
8. Drug repurposing    — Query L1000FWD for drugs that reverse the signature
      ↓
Ranked drug candidates + Markdown report
```

Two steps require human input in the terminal: selecting which GEO dataset to use, and assigning sample groups (disease vs. control). Both decisions are saved and can be reloaded on subsequent runs.

---

## Requirements

- Python 3.10+
- Miniconda or Anaconda (recommended)
- NCBI API key (free — [register here](https://www.ncbi.nlm.nih.gov/account/))
- Orphanet API key (free — [register here](https://api.orphacode.org))
- Broad Drug Repurposing Hub reference files (see Setup below)

---

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/thomascnelson/TDRE.git
cd TDRE
```

### 2. Create the conda environment

```bash
conda create -n tdre python=3.10
conda activate tdre
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Create your `.env` file

Create a file named `.env` in the project root with your API keys:

```
NCBI_API_KEY=your_ncbi_api_key_here
ORPHANET_API_KEY=your_orphanet_api_key_here
EMAIL_ADDRESS=your_email@example.com
```

> Your email address is required by NCBI's Entrez API. It is only sent to NCBI
> as part of API requests and is never stored or shared otherwise.

### 5. Download the Broad Drug Repurposing Hub reference files

Visit [repo-hub.broadinstitute.org/repurposing](https://repo-hub.broadinstitute.org/repurposing)
and download both:
- **Drug-level annotations** (contains drug names, mechanism of action, targets, clinical phase)
- **Sample-level annotations** (contains BRD compound IDs mapped to drug names)

Place both files in the `resources/` directory. The pipeline will find them automatically
regardless of the date suffix in the filename.

```
resources/
├── repo-drug-annotation-YYYYMMDD.txt
└── repo-sample-annotation-YYYYMMDD.txt
```

---

## Usage

```bash
conda activate tdre
python main.py "Disease Name"
```

**Example:**
```bash
python main.py "Friedreich ataxia"
```

The pipeline will prompt you for input at two points:

1. **Dataset selection** — A scored list of GEO datasets is shown. Review the
   metadata and enter the number of the dataset you want to use.

2. **Sample group assignment** — Each sample in the dataset is shown with its
   metadata. Confirm or correct the suggested group assignment (disease / control / exclude).
   Assignments are saved so you won't be asked again on subsequent runs.

### Output

Results are written to disease-specific subdirectories:

```
data/friedreich_ataxia/
├── GSE226646_counts.csv             # Raw count matrix
├── GSE226646_metadata.csv           # Sample metadata
├── GSE226646_groups.csv             # Saved sample group assignments
├── GSE226646_deseq2_results.csv     # Full DESeq2 results
├── GSE226646_signature.json         # Disease signature (for pipeline)
├── GSE226646_signature.csv          # Disease signature (for investigator)
├── GSE226646_drug_candidates.json   # Drug candidates (for pipeline)
└── GSE226646_drug_candidates.csv    # Drug candidates (for investigator)

reports/friedreich_ataxia/
├── friedreich_ataxia_provenance.json   # Machine-readable run log
└── friedreich_ataxia_report.md         # Human-readable report
```

The Markdown report includes an executive summary, dataset description, DEG
statistics, disease signature genes, ranked drug candidates table, QC warnings,
and a full methods section.

---

## Using Individual Tools

Each tool can be run independently from the project root. This is useful when
you want to explore GEO datasets for a disease, re-run a single step with
different parameters, or inspect intermediate results without running the full
pipeline. All tools read API keys from your `.env` file automatically.

> Run all commands from the project root directory with the `tdre` environment active.

---

### Tool 1 — Disease Lookup (`tools/disease_lookup.py`)

Checks whether a disease is a viable drug repurposing candidate by querying
NCBI MeSH, Orphanet, and PubMed. Returns a feasibility score and flags any
issues (e.g. no MeSH term, not in Orphanet, sparse treatment literature).

**Use this when:** You want to quickly assess whether a disease is well enough
characterised to be worth running through the full pipeline.

```bash
python tools/disease_lookup.py
```

Edit the disease name at the bottom of the file before running:
```python
result = lookup_disease("Your Disease Name", provenance=prov)
```

**Output:** Prints MeSH terms, Orphanet match, PubMed treatment literature count,
feasibility score (0–100), and any warning flags.

---

### Tool 2 — GEO Dataset Search (`tools/geo_search.py`)

Searches NCBI GEO for bulk RNA-seq datasets relevant to a disease and scores
each one for suitability (sample size, organism, presence of a control group).
Presents the top recommended datasets interactively for you to select one.

**Use this when:** You want to survey what transcriptomic data exists for a
disease before committing to the full pipeline, or to find datasets for your
own analysis.

```bash
python tools/geo_search.py
```

Edit the disease name at the bottom of the file before running:
```python
geo_search_results = search_geo_datasets("Your Disease Name", max_results=20, provenance=prov)
```

**Output:** Prints a JSON list of scored datasets sorted by suitability score,
then launches the interactive dataset selector.

---

### Tool 3 — GEO Download (`tools/geo_download.py`)

Downloads a specific GEO dataset by accession number. Extracts the raw count
matrix (from supplementary files if available, otherwise from the GEO matrix
file) and sample metadata, then guides you through assigning samples to
disease/control/exclude groups. Group assignments are saved and reloaded on
subsequent runs.

**Use this when:** You already know which GEO accession you want to work with,
or you want to re-do the sample group assignment for a dataset you have
already downloaded.

```bash
python tools/geo_download.py
```

Edit the accession at the bottom of the file before running:
```python
accession = "GSE226646"
```

**Output:** Downloads files to `data/`, saves `{accession}_metadata.csv`,
`{accession}_counts.csv`, and `{accession}_groups.csv`. Prints a summary of
what was found and the sample group assignment results.

---

### Tool 4 — Differential Expression Analysis (`tools/deg_analysis.py`)

Runs DESeq2 (via PyDESeq2) on a count matrix using saved sample group
assignments. Requires raw integer counts. Returns the full DESeq2 results table
with log2 fold changes and adjusted p-values.

**Use this when:** You want to re-run DEG analysis with different thresholds,
or you have a count matrix and group file from outside the pipeline.

```bash
python tools/deg_analysis.py
```

Edit the file paths at the bottom before running:
```python
COUNT_MATRIX = "data/GSE226646_counts.csv"
GROUPS_FILE  = "data/GSE226646_groups.csv"
```

The significance thresholds can be changed by passing parameters to `run_deseq2()`:
```python
result = run_deseq2(COUNT_MATRIX, sample_groups,
                    padj_threshold=0.1,   # default 0.05
                    lfc_threshold=0.5,    # default 1.0
                    min_count_filter=10)  # default 10
```

**Output:** Saves `data/{prefix}_deseq2_results.csv`. Prints gene counts and
top up/downregulated genes.

---

### Tool 5 — Signature Building (`tools/signature.py`)

Filters DESeq2 results to the most significant differentially expressed genes,
maps gene IDs to HGNC symbols via the mygene API, and saves the disease
signature as both JSON (for the pipeline) and CSV (for manual inspection).

**Use this when:** You want to inspect or export the disease signature, or
re-build it from existing DESeq2 results with different thresholds.

```bash
python tools/signature.py
```

Edit the DESeq2 results path at the bottom before running:
```python
DEG_RESULTS = "data/GSE226646_deseq2_results.csv"
```

Thresholds can be adjusted:
```python
signature = build_signature(DEG_RESULTS,
                            padj_threshold=0.1,  # default 0.05
                            lfc_threshold=0.5,   # default 1.0
                            max_genes=150)        # default 150
```

**Output:** Saves `{prefix}_signature.json` and `{prefix}_signature.csv`.
Prints gene counts, quality assessment, and top up/downregulated gene symbols.

---

### Tool 6 — Drug Repurposing (`tools/drug_repurposing.py`)

Submits a disease signature to the L1000FWD API (Ma'ayan Lab, Mount Sinai) and
retrieves drugs whose transcriptomic effects reverse the disease signature.
Drug names and annotations are resolved from the Broad Drug Repurposing Hub.
Falls back to Enrichr if L1000FWD is unavailable.

**Use this when:** You have a gene signature (from any source) and want to find
drugs that reverse it, or you want to re-run the drug query for an existing
signature.

```bash
python tools/drug_repurposing.py
```

Edit the signature file path at the bottom before running:
```python
SIGNATURE_FILE = "data/GSE226646_signature.json"
```

The signature JSON must contain `"upregulated"` and `"downregulated"` lists of
HGNC gene symbols. You can create one manually if needed:
```json
{
  "upregulated":   ["GENE1", "GENE2", "GENE3"],
  "downregulated": ["GENE4", "GENE5", "GENE6"]
}
```

**Output:** Saves `data/{prefix}_drug_candidates.json` and `.csv`. Prints a
ranked table of drug candidates with scores, confidence levels, and mechanism
of action. Also prints a link to the interactive L1000FWD results viewer.

---

## Project Structure

```
tdre/
├── main.py                     # CLI entry point
├── tools/
│   ├── disease_lookup.py       # Tool 1: disease characterisation
│   ├── geo_search.py           # Tool 2: GEO dataset search and selection
│   ├── geo_download.py         # Tool 3: data download and sample grouping
│   ├── deg_analysis.py         # Tool 4: DESeq2 differential expression
│   ├── signature.py            # Tool 5: transcriptomic signature building
│   └── drug_repurposing.py     # Tool 6: L1000FWD drug connectivity scoring
├── utils/
│   └── provenance.py           # Provenance log and report generation
├── resources/
│   └── repo-*-annotation-*.txt # Broad Repurposing Hub reference files
├── data/                       # Pipeline outputs (per-disease subdirectories)
├── reports/                    # Provenance logs and reports
├── requirements.txt
└── .env                        # API keys (not committed to git)
```

---

## Drug Name Resolution

Drug candidates are returned as BRD compound IDs by the L1000FWD API. TDRE
resolves these to common drug names using a three-tier lookup:

1. **Broad Drug Repurposing Hub** (local files) — primary source; also provides
   mechanism of action, molecular target, and clinical development phase
2. **PubChem API** — fallback for compounds not in the Hub
3. **Raw BRD ID** — last resort if both lookups fail

> Note: CLUE.io (the Broad's compound API) is not used as it requires an
> institutional email address for API access.

---

## Interpreting Results

Drug candidates are ranked by L1000FWD combined score. More negative scores
indicate stronger reversal of the disease transcriptomic signature. Confidence
is based on q-value:

| Confidence | Q-value |
|------------|---------|
| High       | < 0.05  |
| Medium     | < 0.20  |
| Low        | ≥ 0.20  |

**All results should be treated as hypotheses for experimental investigation.**
The pipeline identifies gene expression correlations — it does not validate
drug efficacy, safety, or disease relevance. Experimental validation in
disease-relevant cell or animal models is required before drawing conclusions.

---

## License

Mozilla Public License 2.0 — see [LICENSE](LICENSE) for details.

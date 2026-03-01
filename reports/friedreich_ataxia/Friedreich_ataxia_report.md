# Drug Repurposing Report: Friedreich ataxia

**Generated:** 2026-03-01 08:14  
**Pipeline version:** 0.1.0  

---

## Executive Summary

This analysis identified 48 drug repurposing candidates for Friedreich ataxia using transcriptomic data from GSE226646 (Skeletal muscle transcriptomics dissects the pathogenesis of Friedreich's Ataxia). Differential expression analysis compared 7 disease samples against 6 controls, identifying 560 upregulated and 619 downregulated genes. The disease signature comprised 68 upregulated and 39 downregulated genes after filtering and symbol mapping. High-confidence candidates (q < 0.05) include: PP-1, GlyTouCan:G32632JB, tenidap, BRD-K88378636-001-02-8. 68 quality control warning(s) were raised during the analysis — see the QC Warnings section below.

---

## Disease Overview

**Disease:** Friedreich ataxia  
**MeSH ID(s):** 68005621, 67566594, 67565561, 67564999, 67564657, 67564446, 67538061, 67535393, 2106500  
**PubMed treatment literature:** 260 publications  

## Dataset

**GEO Accession:** [GSE226646](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226646)  
**Title:** Skeletal muscle transcriptomics dissects the pathogenesis of Friedreich's Ataxia  
**Organism:** Homo sapiens  
**Total samples:** 20  
**Disease samples:** 7  
**Control samples:** 6  
**Count data source:** supplementary_counts  

**Datasets considered but rejected:**

- `GSE262763` — Correction of splicing defect in compound heterozygous FRDA patient carrying FXN 165+5G>C point mutation.: Insufficient sample size (1 samples) — need at least 6; No clear mention of control group in summary
- `GSE232826` — Loss of homeostatic functions in microglia from a murine model of Friedreich’s ataxia: Borderline sample size (8 samples) — DESeq2 results may be underpowered; Mouse data — gene ortholog mapping will be required; No clear mention of control group in summary
- `GSE261654` — Frataxin Deficiency Drives a Shift from Mitochondrial Metabolism to Glucose Catabolism, Triggering an Inflammatory Phenotype in Microglia (scRNA-Seq): Insufficient sample size (2 samples) — need at least 6; Mouse data — gene ortholog mapping will be required
- `GSE242192` — Genetic interaction screen for FXN in K562 cells, to identify genes that are synthetic sick/buffering in the absence of robust Fe-S cluster synthesis: Insufficient sample size (4 samples) — need at least 6
- `GSE239980` — Butyrate prevents visceral adipose tissue inflammation and metabolic alterations in a mouse model of Friedreich’s ataxia: Borderline sample size (8 samples) — DESeq2 results may be underpowered; Mouse data — gene ortholog mapping will be required; No clear mention of control group in summary
- `GSE176575` — Disruption of the histone H3 cupric reductase activity prevents global transcriptional rewiring when YFH1 is diminished: No clear mention of control group in summary
- `GSE148212` — CRISPR/Cas9 targeted deletions in frataxin rescue Friedreich’s ataxia pathology in a human stem cell-based peripheral sensory neural circuitry: Insufficient sample size (2 samples) — need at least 6; No clear mention of control group in summary
- `GSE133755` — In vitro restoration of Friedreich's Ataxia defect with novel Cas9-edited hiPSC-derived Dorsal Root Ganglion-like Organoids (DRGOs) and Microfluidic System: Insufficient sample size (3 samples) — need at least 6; No clear mention of control group in summary
- `GSE108200` — Frataxin deficiency in Friedreich?s Ataxia is associated with reduced levels of HAX-1, a regulator of cardiomyocyte death and survival.: Insufficient sample size (2 samples) — need at least 6; No clear mention of control group in summary
- `GSE109510` — The transcriptional regulator CCCTC-binding factor limits oxidative stress in endothelial cells: Borderline sample size (6 samples) — DESeq2 results may be underpowered; Mouse data — gene ortholog mapping will be required; No clear mention of control group in summary
- `GSE15843` — Functional genomic analysis of frataxin deficiency, Agilent data: Mouse data — gene ortholog mapping will be required; No clear mention of control group in summary
- `GSE11425` — Chromosome fragility at the expanded GAA/TTC tracks depends on repeat orientation and requires mismatch repair system.: No clear mention of control group in summary
- `GSM7797185` — PAT12 (LOG): Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM7797184` — PAT11 (LOG): Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM7797183` — PAT10 (LOG): Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM7797182` — PAT9 (LOG): Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM7797181` — PAT8 (IOG): Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM7797180` — PAT7 (IOG): Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM7797179` — PAT6 (IOG): Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM7797178` — PAT5 (IOG): Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM7797177` — PAT4 (EOG): Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM7797176` — PAT3 (EOG): Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary

## Sample Group Assignment

| Group | N |
|-------|---|
| disease | 7 |
| exclude | 7 |
| control | 6 |

## Differential Gene Expression Analysis

**Method:** DESeq2 (PyDESeq2)  
**Genes tested:** 14592  
**Significantly upregulated:** 560  
**Significantly downregulated:** 619  

**Parameters:**

- Adjusted p-value threshold: 0.05
- log₂ fold-change threshold: 1.0
- Minimum count filter: 10

## Disease Transcriptomic Signature

**Quality:** good  
**Upregulated genes:** 68  
**Downregulated genes:** 39  

**Top upregulated:**  
`DENND5B`, `MYOM1`, `NDRG2`, `PIK3R1`, `NUCKS1`, `MIR133A1HG`, `RBM25`, `YTHDF3`, `EIF5B`, `CAVIN1`, `PTMA`, `EZR`, `RERGL`, `ZDHHC21`, `BTBD7`

**Top downregulated:**  
`COX6A2`, `GABARAP`, `HLA-B`, `VPS28`, `MT-CO1`, `MRPS24`, `JTB`, `SLC29A1`, `NPAS2`, `ABHD14A`, `HCFC1`, `EFHD1`, `DBET`, `ATP5MK`, `ATP5F1D`

## Drug Repurposing Candidates

Drugs are ranked by L1000FWD combined score (more negative = stronger reversal of disease signature).  

**Interactive results:** [https://maayanlab.cloud/l1000fwd/#result/69a43bbcd33550002d585915](https://maayanlab.cloud/l1000fwd/#result/69a43bbcd33550002d585915)  

| Rank | Drug | Score | Q-value | Confidence | Cell Line | MoA | Target |
|------|------|-------|---------|------------|-----------|-----|--------|
| 1 | PP-1 | -11.9913 | 0.00383294839283436 | high | PC3 | SRC inhibitor | HCK | RET | SRC |
| 2 | GlyTouCan:G32632JB | -9.0777 | 0.0400531231230799 | high | PC3 |  |  |
| 3 | tenidap | -8.8274 | 0.0414257784104472 | high | VCAP | cyclooxygenase inhibitor | ALOX5 | KCNJ4 | PTGS1 |
| 4 | BRD-K88378636-001-02-8 | -7.7022 | 0.0435657138716559 | high | VCAP |  |  |
| 5 | BRD-K25373946-001-04-6 | -7.6479 | 0.0957349706036802 | medium | PC3 |  |  |
| 6 | GANT-58 | -7.2845 | 0.0957349706036802 | medium | MCF7 | GLI antagonist | DHH | GLI1 | IHH |
| 7 | mafenide | -7.0878 | 0.0957349706036802 | medium | MCF7 | carbonic anhydrase inhibitor | CA12 | CA14 | CA2 | CA4 | CA6 | CA9 |
| 8 | BRD-K66792149-001-01-8 | -7.0669 | 0.0957349706036802 | medium | PC3 |  |  |
| 9 | BRD-K32644160-001-16-4 | -6.9396 | 0.122410114907044 | medium | PC3 |  |  |
| 10 | Ro-28-1675 | -6.9275 | 0.0958041294275837 | medium | HCC15 | glucokinase activator | GCK |

> **Interpretation note:** Candidates should be treated as hypotheses for further investigation, not validated drug targets. Experimental validation in disease-relevant cell or animal models is required before drawing conclusions.

## Quality Control Warnings

The following issues were flagged during the analysis:

- **geo_search** (2026-03-01): GSE262763: Insufficient sample size (1 samples) — need at least 6
- **geo_search** (2026-03-01): GSE262763: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE255810: Borderline sample size (6 samples) — DESeq2 results may be underpowered
- **geo_search** (2026-03-01): GSE255810: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE271798: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE232826: Borderline sample size (8 samples) — DESeq2 results may be underpowered
- **geo_search** (2026-03-01): GSE232826: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE232826: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE261655: Borderline sample size (8 samples) — DESeq2 results may be underpowered
- **geo_search** (2026-03-01): GSE261655: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE261654: Insufficient sample size (2 samples) — need at least 6
- **geo_search** (2026-03-01): GSE261654: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE261653: Borderline sample size (6 samples) — DESeq2 results may be underpowered
- **geo_search** (2026-03-01): GSE261653: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE242192: Insufficient sample size (4 samples) — need at least 6
- **geo_search** (2026-03-01): GSE239980: Borderline sample size (8 samples) — DESeq2 results may be underpowered
- **geo_search** (2026-03-01): GSE239980: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE239980: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE229143: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE200907: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE204808: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE203632: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE176575: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE145115: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE148212: Insufficient sample size (2 samples) — need at least 6
- **geo_search** (2026-03-01): GSE148212: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE133755: Insufficient sample size (3 samples) — need at least 6
- **geo_search** (2026-03-01): GSE133755: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE108200: Insufficient sample size (2 samples) — need at least 6
- **geo_search** (2026-03-01): GSE108200: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE102008: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE109510: Borderline sample size (6 samples) — DESeq2 results may be underpowered
- **geo_search** (2026-03-01): GSE109510: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE109510: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE98790: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE98790: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE105052: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE65399: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE31208: Borderline sample size (8 samples) — DESeq2 results may be underpowered
- **geo_search** (2026-03-01): GSE31208: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE30933: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE22651: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE15848: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE15848: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE15843: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE15843: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE11425: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE10745: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSM7797185: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM7797185: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM7797184: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM7797184: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM7797183: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM7797183: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM7797182: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM7797182: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM7797181: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM7797181: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM7797180: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM7797180: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM7797179: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM7797179: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM7797178: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM7797178: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM7797177: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM7797177: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM7797176: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM7797176: No clear mention of control group in summary

## Output Files

| File | Path |
|------|------|
| disease_lookup_raw | `data/friedreich_ataxia/Friedreich_ataxia_lookup.json` |
| count_matrix | `data/friedreich_ataxia/GSE226646_counts.csv` |
| sample_metadata | `data/friedreich_ataxia/GSE226646_metadata.csv` |
| sample_groups | `data/friedreich_ataxia/GSE226646_groups.csv` |
| deg_results | `data/friedreich_ataxia/GSE226646_deseq2_results.csv` |
| signature_json | `data/friedreich_ataxia/GSE226646_signature.json` |
| signature_csv | `data/friedreich_ataxia/GSE226646_signature.csv` |
| drug_candidates_json | `data/friedreich_ataxia/GSE226646_drug_candidates.json` |
| drug_candidates_csv | `data/friedreich_ataxia/GSE226646_drug_candidates.csv` |

## Methods

**Transcriptomic data:** Retrieved from NCBI Gene Expression Omnibus (GEO). Raw count matrices were downloaded and quality-assessed prior to analysis.  

**Differential expression:** Performed using PyDESeq2, a Python implementation of the DESeq2 negative binomial model. Genes were filtered to a minimum count sum of 10 across all samples. Significance was assessed using Wald tests with Benjamini-Hochberg FDR correction.  

**Drug repurposing:** Disease signature genes were submitted to L1000FWD (Ma'ayan Laboratory, Mount Sinai) for connectivity scoring against the LINCS L1000 chemical perturbagen database. Drugs with negatively concordant signatures (i.e. those that reverse the disease gene expression pattern) were ranked by combined score.  

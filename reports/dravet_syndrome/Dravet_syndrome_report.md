# Drug Repurposing Report: Dravet syndrome

**Generated:** 2026-03-01 07:37  
**Pipeline version:** 0.1.0  

---

## Executive Summary

This analysis identified 48 drug repurposing candidates for Dravet syndrome using transcriptomic data from GSE256142 (Variability vs. phenotype: Multimodal analysis of Dravet syndrome brain organoids powered by deep learning). Differential expression analysis compared 5 disease samples against 8 controls, identifying 1608 upregulated and 1670 downregulated genes. The disease signature comprised 57 upregulated and 49 downregulated genes after filtering and symbol mapping. High-confidence candidates (q < 0.05) include: alvespimycin, tretinoin, BRD-K35716340-001-02-0, BRD-A79465854-001-01-3, BIBR-1532. 51 quality control warning(s) were raised during the analysis — see the QC Warnings section below.

---

## Disease Overview

**Disease:** Dravet syndrome  
**MeSH ID(s):** 68004831  
**PubMed treatment literature:** 770 publications  

## Dataset

**GEO Accession:** [GSE256142](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256142)  
**Title:** Variability vs. phenotype: Multimodal analysis of Dravet syndrome brain organoids powered by deep learning  
**Organism:** Homo sapiens  
**Total samples:** 18  
**Disease samples:** 8  
**Control samples:** 10  
**Count data source:** supplementary_counts  

**Datasets considered but rejected:**

- `GSE301399` — Brainstem metabolic changes in a mouse model of Dravet Syndrome: Insufficient sample size (2 samples) — need at least 6; Mouse data — gene ortholog mapping will be required; No clear mention of control group in summary
- `GSE245113` — Potentiating NaV1.1 in Dravet Syndrome patient iPSC-derived GABAergic neurons increases neuronal firing frequency and decreases network synchrony: Insufficient sample size (1 samples) — need at least 6
- `GSE169485` — Deletion of a non-canonical regulatory sequence causes loss of Scn1a expression and epileptic phenotypes in mice II: Borderline sample size (8 samples) — DESeq2 results may be underpowered; Mouse data — gene ortholog mapping will be required; No clear mention of control group in summary
- `GSE169481` — Deletion of a non-canonical regulatory sequence causes loss of Scn1a expression and epileptic phenotypes in mice: Mouse data — gene ortholog mapping will be required; No clear mention of control group in summary
- `GSE145801` — EPI-00649-2019 New insights into the early mechanisms of epileptogenesis in a zebrafish model of Dravet syndrome: Insufficient sample size (4 samples) — need at least 6; No clear mention of control group in summary
- `GSE111436` — dCas9-based Scn1a gene activation restores inhibitory interneuron excitability and restrains epileptic crises in Dravet syndrome mice: Borderline sample size (6 samples) — DESeq2 results may be underpowered; Mouse data — gene ortholog mapping will be required; No clear mention of control group in summary
- `GSE43312` — Drug screening for Dravet syndrome using Scn1a zebrafish mutants: Borderline sample size (8 samples) — DESeq2 results may be underpowered; No clear mention of control group in summary
- `GSM9082361` — Scn1a knock-in biol rep 1: Insufficient sample size (0 samples) — need at least 6; Mouse data — gene ortholog mapping will be required; No clear mention of control group in summary
- `GSM9082360` — wildtype biol rep 1: Insufficient sample size (0 samples) — need at least 6; Mouse data — gene ortholog mapping will be required; No clear mention of control group in summary
- `GSM8086398` — DS_iPSC_d90_VPA_n_2: Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM8086397` — DS_iPSC_d90_VPA_n_1: Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM8086396` — DS_iPSC_d90_n_2: Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM8086395` — DS_iPSC_d90_n_1: Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM8086394` — DS_iPSC_d35_n_2: Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM8086393` — DS_iPSC_d35_n_1: Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM8086389` — DS_iPSC_d35_n_3: Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM8086387` — DS_iPSC_d90_n_3: Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary
- `GSM7836933` — DIV 21 Dravet Patient iPSC neurons: Insufficient sample size (0 samples) — need at least 6; No clear mention of control group in summary

## Sample Group Assignment

| Group | N |
|-------|---|
| control | 10 |
| disease | 8 |

## Differential Gene Expression Analysis

**Method:** DESeq2 (PyDESeq2)  
**Genes tested:** 26403  
**Significantly upregulated:** 1608  
**Significantly downregulated:** 1670  

**Parameters:**

- Adjusted p-value threshold: 0.05
- log₂ fold-change threshold: 1.0
- Minimum count filter: 10

## Disease Transcriptomic Signature

**Quality:** good  
**Upregulated genes:** 57  
**Downregulated genes:** 49  

**Top upregulated:**  
`CACNB4`, `HS3ST2`, `SYNDIG1`, `TTTY14`, `CCDC112`, `ZNF439`, `PCDHA6`, `SMIM17`, `ETNK1`, `PRKY`, `TMEM163`, `CAMTA1`, `PCDHB18P`, `NLGN4Y`, `BTNL9`

**Top downregulated:**  
`LRRC61`, `ALG10B`, `ARL17B`, `SLC2A1`, `SH3D19`, `ARL17A`, `YBX3`, `POTEF`, `LGR4`, `RAX`, `SIK1`, `WDR72`, `IL17RB`, `MAN2B1`, `ZNF37A`

## Drug Repurposing Candidates

Drugs are ranked by L1000FWD combined score (more negative = stronger reversal of disease signature).  

**Interactive results:** [https://maayanlab.cloud/l1000fwd/#result/69a432efd33550002d585914](https://maayanlab.cloud/l1000fwd/#result/69a432efd33550002d585914)  

| Rank | Drug | Score | Q-value | Confidence | Cell Line | MoA | Target |
|------|------|-------|---------|------------|-----------|-----|--------|
| 1 | alvespimycin | -11.5519 | 0.00114569675489634 | high | PC3 | HSP inhibitor | HSP90AA1 |
| 2 | tretinoin | -10.9953 | 0.00149273478441859 | high | MCF7 | retinoid receptor agonist | retinoid receptor ligand | ALDH1A1 | ALDH1A2 | GPRC5A | NR0B1 | NR2C2 | PPARD | RARA | RARB | RARG | RARRES1 | RORB | RORC | RXRB | RXRG |
| 3 | BRD-K35716340-001-02-0 | -10.7093 | 0.00149273478441859 | high | HEC108 |  |  |
| 4 | BRD-A79465854-001-01-3 | -10.4845 | 0.00149273478441859 | high | HT29 |  |  |
| 5 | BIBR-1532 | -10.3899 | 0.00149273478441859 | high | A673 | telomerase inhibitor | TERT |
| 6 | Amg-9810 | -10.3604 | 0.00149273478441859 | high | VCAP | TRPV antagonist | TRPV1 |
| 7 | motesanib | -10.1806 | 0.00149273478441859 | high | MCF7 | KIT inhibitor | PDGFR tyrosine kinase receptor inhibitor | VEGFR inhibitor | FLT1 | FLT4 | KDR | KIT | PDGFRA | RET |
| 8 | (Z)-3-((1H-pyrrol-2-yl)methylene... | -10.0969 | 0.00332922730570087 | high | VCAP |  |  |
| 9 | BIIB021 | -9.9935 | 0.00149273478441859 | high | HA1E | HSP inhibitor | AOX1 | HSP90AA1 |
| 10 | N'-benzoyl-2-phenyl-4-thiazoleca... | -9.9254 | 0.00149273478441859 | high | HA1E |  |  |

> **Interpretation note:** Candidates should be treated as hypotheses for further investigation, not validated drug targets. Experimental validation in disease-relevant cell or animal models is required before drawing conclusions.

## Quality Control Warnings

The following issues were flagged during the analysis:

- **geo_search** (2026-03-01): GSE301399: Insufficient sample size (2 samples) — need at least 6
- **geo_search** (2026-03-01): GSE301399: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE301399: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE256142: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE280241: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE245113: Insufficient sample size (1 samples) — need at least 6
- **geo_search** (2026-03-01): GSE171191: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE169485: Borderline sample size (8 samples) — DESeq2 results may be underpowered
- **geo_search** (2026-03-01): GSE169485: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE169485: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE169481: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE169481: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE153461: Borderline sample size (8 samples) — DESeq2 results may be underpowered
- **geo_search** (2026-03-01): GSE153461: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE145801: Insufficient sample size (4 samples) — need at least 6
- **geo_search** (2026-03-01): GSE145801: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE136535: Borderline sample size (8 samples) — DESeq2 results may be underpowered
- **geo_search** (2026-03-01): GSE136535: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE128573: Borderline sample size (8 samples) — DESeq2 results may be underpowered
- **geo_search** (2026-03-01): GSE128573: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE111436: Borderline sample size (6 samples) — DESeq2 results may be underpowered
- **geo_search** (2026-03-01): GSE111436: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE111436: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE112627: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSE112627: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSE43312: Borderline sample size (8 samples) — DESeq2 results may be underpowered
- **geo_search** (2026-03-01): GSE43312: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM9082361: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM9082361: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSM9082361: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM9082360: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM9082360: Mouse data — gene ortholog mapping will be required
- **geo_search** (2026-03-01): GSM9082360: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM8086398: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM8086398: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM8086397: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM8086397: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM8086396: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM8086396: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM8086395: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM8086395: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM8086394: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM8086394: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM8086393: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM8086393: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM8086389: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM8086389: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM8086387: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM8086387: No clear mention of control group in summary
- **geo_search** (2026-03-01): GSM7836933: Insufficient sample size (0 samples) — need at least 6
- **geo_search** (2026-03-01): GSM7836933: No clear mention of control group in summary

## Output Files

| File | Path |
|------|------|
| disease_lookup_raw | `data/dravet_syndrome/Dravet_syndrome_lookup.json` |
| count_matrix | `data/dravet_syndrome/GSE256142_counts.csv` |
| sample_metadata | `data/dravet_syndrome/GSE256142_metadata.csv` |
| sample_groups | `data/dravet_syndrome/GSE256142_groups.csv` |
| deg_results | `data/dravet_syndrome/GSE256142_deseq2_results.csv` |
| signature_json | `data/dravet_syndrome/GSE256142_signature.json` |
| signature_csv | `data/dravet_syndrome/GSE256142_signature.csv` |
| drug_candidates_json | `data/dravet_syndrome/GSE256142_drug_candidates.json` |
| drug_candidates_csv | `data/dravet_syndrome/GSE256142_drug_candidates.csv` |

## Methods

**Transcriptomic data:** Retrieved from NCBI Gene Expression Omnibus (GEO). Raw count matrices were downloaded and quality-assessed prior to analysis.  

**Differential expression:** Performed using PyDESeq2, a Python implementation of the DESeq2 negative binomial model. Genes were filtered to a minimum count sum of 10 across all samples. Significance was assessed using Wald tests with Benjamini-Hochberg FDR correction.  

**Drug repurposing:** Disease signature genes were submitted to L1000FWD (Ma'ayan Laboratory, Mount Sinai) for connectivity scoring against the LINCS L1000 chemical perturbagen database. Drugs with negatively concordant signatures (i.e. those that reverse the disease gene expression pattern) were ranked by combined score.  

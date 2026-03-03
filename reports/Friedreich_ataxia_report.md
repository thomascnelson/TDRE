# Drug Repurposing Report: Friedreich ataxia

**Generated:** 2026-03-03 10:11  
**Pipeline version:** 0.1.0  

---

## Disease Overview

**Disease:** Friedreich ataxia  

## Dataset

**GEO Accession:** [GSE226646](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226646)  
**Title:** Skeletal muscle transcriptomics dissects the pathogenesis of Friedreich's Ataxia  
**Organism:**   
**Total samples:** 20  
**Disease samples:** 7  
**Control samples:** 6  
**Count data source:** supplementary_counts  

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
`TNPO1`, `MMS19`, `NUCKS1`, `DENND5B`, `BPTF`, `SMARCC1`, `CEBPZ`, `PRMT1`, `TGFB2`, `NDRG2`, `ZBTB20`, `P2RY1`, `EZR`, `ASPH`, `APPL1`

**Top downregulated:**  
`ADD1`, `PFDN5`, `ATP5F1D`, `FBLN2`, `COX6A2`, `HCFC1`, `ATP5F1A`, `ADH1B`, `EFHD1`, `MRPL36`, `IMPA2`, `MRLN`, `NPAS2`, `SLC29A1`, `CHCHD10`

## Drug Repurposing Candidates

Drugs are ranked by L1000FWD combined score (more negative = stronger reversal of disease signature).  

**Interactive results:** [https://maayanlab.cloud/l1000fwd/#result/69a6fa20d33550002d58592b](https://maayanlab.cloud/l1000fwd/#result/69a6fa20d33550002d58592b)  

| Rank | Drug | Score | Q-value | Confidence | Phase | Cell Line | MoA | Target |
|------|------|-------|---------|------------|-------|-----------|-----|--------|
| 1 | PP-1 | -11.9913 | 0.00383294839283436 | high | Preclinical | PC3 | SRC inhibitor | HCK | RET | SRC |
| 2 | GlyTouCan:G32632JB | -9.0777 | 0.0400531231230799 | high |  | PC3 |  |  |
| 3 | tenidap | -8.8274 | 0.0414257784104472 | high | Phase 3 | VCAP | cyclooxygenase inhibitor | ALOX5 | KCNJ4 | PTGS1 |
| 4 | BRD-K88378636-001-02-8 | -7.7022 | 0.0435657138716559 | high |  | VCAP |  |  |
| 5 | BRD-K25373946-001-04-6 | -7.6479 | 0.0957349706036802 | medium |  | PC3 |  |  |
| 6 | GANT-58 | -7.2845 | 0.0957349706036802 | medium | Preclinical | MCF7 | GLI antagonist | DHH | GLI1 | IHH |
| 7 | mafenide | -7.0878 | 0.0957349706036802 | medium | Launched | MCF7 | carbonic anhydrase inhibitor | CA12 | CA14 | CA2 | CA4 | CA6 | CA9 |
| 8 | BRD-K66792149-001-01-8 | -7.0669 | 0.0957349706036802 | medium |  | PC3 |  |  |
| 9 | BRD-K32644160-001-16-4 | -6.9396 | 0.122410114907044 | medium |  | PC3 |  |  |
| 10 | Ro-28-1675 | -6.9275 | 0.0958041294275837 | medium | Preclinical | HCC15 | glucokinase activator | GCK |
| 11 | (S)-2-((9-([1,1'-biphenyl]-4-ylm... | -6.8179 | 0.143674306134858 | medium |  | HT29 |  |  |
| 12 | Awl-II-38.3 | -6.6564 | 0.0983899927902734 | medium |  | PC3 |  |  |
| 13 | Acetic acid [2-(4-acetamido-6-ph... | -6.5439 | 0.0957349706036802 | medium |  | HT29 |  |  |
| 14 | RG 13022 | -6.3527 | 0.170634907488047 | medium |  | VCAP |  |  |
| 15 | YK-4-279 | -6.2573 | 0.0958041294275837 | medium | Preclinical | VCAP | apoptosis inhibitor | EWSR1 | FLI1 |
| 16 | 10-[(Dimethylamino)methyl]-4-eth... | -6.24 | 0.170634907488047 | medium |  | HCC515 |  |  |
| 17 | N'-(2-aminophenyl)-N-[(2R,3S)-5-... | -6.1963 | 0.102431627040414 | medium |  | PC3 |  |  |
| 18 | BRD-K36153907-001-03-8 | -6.1915 | 0.170634907488047 | medium |  | VCAP |  |  |
| 19 | batimastat | -5.9743 | 0.170634907488047 | medium | Phase 3 | HCC515 | matrix metalloprotease inhibitor | ADAM28 | ADAMTS5 | MMP1 | MMP12 | MMP16 | MMP2 | MMP8 |
| 20 | dovitinib | -5.962 | 0.170634907488047 | medium | Phase 3 | A375 | EGFR inhibitor | FGFR inhibitor | FLT3 inhibitor | PDGFR tyrosine kinase receptor inhibitor | VEGFR inhibitor | CSF1R | EGFR | FGFR1 | FGFR2 | FGFR3 | FLT1 | FLT3 | FLT4 | INSR | KDR | KIT | PDGFRA | PDGFRB |
| 21 | BRD-K10882151-001-06-8 | -5.9584 | 0.170634907488047 | medium |  | A375 |  |  |
| 22 | selumetinib | -5.8764 | 0.170634907488047 | medium | Phase 3 | HCC15 | MEK inhibitor | MAP2K1 | MAP2K2 |
| 23 | BRD-A02481876-001-09-9 | -5.7224 | 0.170634907488047 | medium |  | PC3 |  |  |
| 24 | (4R,8S,9S,11S,12R,13S,19S)-12,19... | -5.7063 | 0.170634907488047 | medium |  | A549 |  |  |
| 25 | 3-(1H-indazol-6-yl)-1-propan-2-y... | -5.6703 | 0.170634907488047 | medium |  | HCC515 |  |  |
| 26 | 7-Ethyl-10-hydroxy camptothecin,... | -5.6036 | 0.170634907488047 | medium |  | RMUGS |  |  |
| 27 | 1-((3,4-Difluorophenyl)methyl)-2... | -5.5877 | 0.170634907488047 | medium |  | HT29 |  |  |
| 28 | lovastatin | -5.5545 | 0.170634907488047 | medium | Launched | HA1E | HMGCR inhibitor | CYP3A5 | HDAC2 | HMGCR | ITGAL | NR1I2 |
| 29 | (7S,9S)-7-[[(4S,5S,6S)-4-amino-5... | -5.5051 | 0.170634907488047 | medium |  | A375 |  |  |
| 30 | Erythrosin(E) | -5.4486 | 0.17314317705545 | medium |  | HA1E |  |  |
| 31 | YC-1 | -5.416 | 0.170634907488047 | medium | Preclinical | HUH7 | guanylyl cyclase activator | GUCY1A2 | GUCY1A3 | GUCY1B3 | HIF1A |
| 32 | Amberlite CG 400 | -5.4145 | 0.170634907488047 | medium |  | MCF7 |  |  |
| 33 | BRD-K56653679-001-01-2 | -5.3409 | 0.170634907488047 | medium |  | HA1E |  |  |
| 34 | ITK inhibitor | -5.3187 | 0.170634907488047 | medium |  | MCF7 |  |  |
| 35 | BRD-A74904029-001-05-7 | -5.2894 | 0.170634907488047 | medium |  | MCF7 |  |  |
| 36 | BRD-K63569039-001-04-0 | -5.2861 | 0.176113814206203 | medium |  | VCAP |  |  |
| 37 | CGP-60474 | -5.2293 | 0.170634907488047 | medium | Preclinical | SKBR3 | CDK inhibitor | CDK1 | CDK2 |
| 38 | (3S,3'S,4'R,5'S)-4'-[fluoro(dime... | -5.1908 | 0.182140433300405 | medium |  | MCF7 |  |  |
| 39 | alvocidib | -5.1894 | 0.170634907488047 | medium | Phase 2 | SKBR3 | CDK inhibitor | BCL2 | BIRC5 | CCNT1 | CDK1 | CDK2 | CDK4 | CDK5 | CDK6 | CDK7 | CDK8 | CDK9 | EGFR | MCL1 | PYGM | XIAP |
| 40 | (1R,9S,12S,15R,16Z,18R,19R,21R,2... | -5.041 | 0.18321523258115 | medium |  | PC3 |  |  |
| 41 | beclomethasone-dipropionate | -4.9384 | 0.176684531266468 | medium | Launched | A375 | glucocorticoid receptor agonist | CYP3A5 | GPR97 | NR3C1 |
| 42 | NSC-95397 | -4.6499 | 0.292575571944643 | low | Preclinical | PC3 | CDC inhibitor | CASP3 | CDC25A | CDC25B |
| 43 | PI-103 | -4.5155 | 0.292575571944643 | low | Preclinical | HCC515 | mTOR inhibitor | PI3K inhibitor | MTOR | PIK3CA | PIK3CB | PIK3CD | PIK3CG | PRKDC |
| 44 | calcipotriol | -4.4439 | 0.292575571944643 | low | Launched | PC3 | vitamin D receptor agonist | VDR |
| 45 | PP-2 | -4.4316 | 0.292575571944643 | low | Preclinical | PC3 | SRC inhibitor | ABL1 | LCK | RIPK2 | SRC |
| 46 | WZ-3146 | -4.227 | 0.292575571944643 | low | Preclinical | ASC | EGFR inhibitor | EGFR |
| 47 | indirubin-3-monoxime | -4.2237 | 0.292575571944643 | low | Preclinical | PC3 | CDK inhibitor | glycogen synthase kinase inhibitor | CCNE1 | CDK1 | CDK2 | CDK4 | CDK5 | CDK5R1 | CDK9 | GSK3B | LCK | LRRK1 | LRRK2 |
| 48 | trichostatin-a | -4.1337 | 0.292575571944643 | low | Phase 1 | MCF7 | HDAC inhibitor | HDAC1 | HDAC10 | HDAC2 | HDAC3 | HDAC4 | HDAC5 | HDAC6 | HDAC7 | HDAC8 | HDAC9 |

> **Interpretation note:** Candidates should be treated as hypotheses for further investigation, not validated drug targets. Experimental validation in disease-relevant cell or animal models is required before drawing conclusions.

## Output Files

| File | Path |
|------|------|
| count_matrix | `data/GSE226646/GSE226646_counts.csv` |
| sample_metadata | `data/GSE226646/GSE226646_metadata.csv` |
| sample_groups | `data/GSE226646/GSE226646_groups.csv` |
| deg_results | `data/GSE226646/GSE226646_deseq2_results.csv` |
| signature_json | `data/GSE226646/GSE226646_signature.json` |
| signature_csv | `data/GSE226646/GSE226646_signature.csv` |
| drug_candidates_json | `data/GSE226646/GSE226646_drug_candidates.json` |
| drug_candidates_csv | `data/GSE226646/GSE226646_drug_candidates.csv` |

## Methods

**Transcriptomic data:** Retrieved from NCBI Gene Expression Omnibus (GEO). Raw count matrices were downloaded and quality-assessed prior to analysis.  

**Differential expression:** Performed using PyDESeq2, a Python implementation of the DESeq2 negative binomial model. Genes were filtered to a minimum count sum of 10 across all samples. Significance was assessed using Wald tests with Benjamini-Hochberg FDR correction.  

**Drug repurposing:** Disease signature genes were submitted to L1000FWD (Ma'ayan Laboratory, Mount Sinai) for connectivity scoring against the LINCS L1000 chemical perturbagen database. Drugs with negatively concordant signatures (i.e. those that reverse the disease gene expression pattern) were ranked by combined score.  

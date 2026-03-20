# Multi-omics pleiotropic association analysis for ovarian aging

This repository contains the analysis workflow and scripts used in the study:

**“Multi-omics pleiotropic association analyses reveal functionally relevant genes and druggable pathways for ovarian aging.”**

---

## Overview

Ovarian aging is characterized by the gradual decline in both the quantity and quality of oocytes and represents a major determinant of female reproductive lifespan.

In this study, we integrated genome-wide association study (GWAS) summary statistics with multiple omics datasets to identify molecular traits associated with ovarian aging.

The analysis includes:

- Proteomics (pQTL)
- Transcriptomics (eQTL)
- Splicing QTL (sQTL)
- Metabolomics (mQTL)

Using Mendelian randomization (SMR) and cross-omics integration (OPERA), we prioritized candidate genes and biological pathways associated with ovarian aging.

---

## Data availability

All datasets used in this study are publicly available:

### GWAS (Age at Natural Menopause, ANM)

- REPROGEN Consortium (Ruth et al.): https://www.reprogen.org/
- UK Biobank replication dataset: https://gwas.mrcieu.ac.uk/datasets/ukb-b-17422/

### pQTL

- deCODE study: https://www.decode.com
- INTERVAL study (replication): https://doi.org/10.1038/s41586-018-0175-2

### eQTL

- eQTLGen consortium: https://www.eqtlgen.org/
- GTEx v8 (ovary tissue): https://www.gtexportal.org

### sQTL

- BrainMeta v2: https://www.nature.com/articles/s41588-022-01154-4
- GTEx v8 (ovary tissue): https://www.gtexportal.org/home/

### mQTL

- CLSA study:https://www.clsa-elcv.ca/
- IEU OpenGWAS datasets:https://gwas.mrcieu.ac.uk/datasets/

### Reference panel

- 1000 Genomes Project Phase 3 (European population)

> Note: Due to file size limitations, raw datasets are not hosted in this repository.  
> Users should download them from the original sources listed above.

---

## ⚙️ Software requirements

- SMR v1.3.1  
- OPERA  
- PLINK v1.9  
- R (≥ 4.2)

Recommended environment: Linux (Ubuntu 20.04 or later)

---

## 🔬 Analysis workflow

The computational pipeline consists of the following steps:

1. GWAS summary statistics preprocessing  
2. QTL data formatting and BESD conversion, SMR analysis and HEIDI test
3. Multi-omics integration using OPERA  
5. Downstream analyses: Functional enrichment analysis, Single-cell analysis, Functional validation.

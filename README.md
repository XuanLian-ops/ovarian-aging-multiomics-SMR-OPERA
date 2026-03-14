# Multi-omics pleiotropic association analysis for ovarian aging

This repository contains the analysis scripts used in the study:

**“Multi-omics pleiotropic association analyses reveal functionally relevant genes and druggable pathways for ovarian aging.”**

## Overview

Ovarian aging is characterized by the gradual decline in both the number and quality of oocytes and represents a major determinant of female reproductive lifespan.

In this study, we integrated genome-wide association study (GWAS) summary statistics with multiple omics datasets to identify molecular traits associated with ovarian aging.

The analysis includes:

- proteomics (pQTL)
- transcriptomics (eQTL)
- splicing QTL (sQTL)
- metabolomics (mQTL)

Using Mendelian randomization and cross-omics association approaches, we prioritized candidate genes and molecular pathways involved in ovarian aging.

## Analysis workflow

The computational pipeline includes:

1. GWAS summary statistics preparation
2. Conversion of QTL summary statistics to BESD format
3. Summary-data-based Mendelian Randomization (SMR)
4. HEIDI test
5. Multi-omics integration using OPERA
6. Functional enrichment and downstream analyses

## Software

The analyses were conducted using:

- SMR v1.3.1
- OPERA
- PLINK v1.9
- R

Reference panel:

- 1000 Genomes Project European population